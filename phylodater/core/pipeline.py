"""
DatingPipeline - 工作流引擎

整合所有步骤，提供高阶 API
"""

import atexit
import hashlib
import os
import random
import sys
import threading
import time
from dataclasses import dataclass, field, is_dataclass
from json import dumps as _json_dumps
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

from ..infrastructure import (
    Configuration,
    RunStatus,
    RuntimeMetadataManager,
    SoftwarePaths,
    ToolConfig,
    get_logger,
)
from ..infrastructure.checkpoint import PipelineCheckpointManager
from ..infrastructure.plugins import QCPlugin, QCPluginChain, get_plugin_registry
from ..models import CalibrationPoint, DatingResult, PhylogeneticTree
from ..services import CalibrationResolver, TaxonomyParser, TreeValidator
from .exceptions import CalibrationError, PhyloDaterError
from .method_interface import DatingMethodRegistry

# 全局临时目录注册表，用于在异常退出时清理（线程安全）
_temp_directories: Set[Path] = set()
_temp_dirs_lock = threading.Lock()


def _derive_method_seeds(
    base_seed: Optional[int], method_names: Sequence[str]
) -> Dict[str, int]:
    """为每个并行方法分配互异且可复现的 32-bit 子种子。

    回归修复（审阅报告 P0-1）：``SeedSequence(base).spawn(n)`` 返回的子对象其
    ``.entropy`` 属性**全部等于父种子**，旧代码读取 ``.entropy`` 导致七个方法
    拿到同一个随机种子（"seed degenerates to the same value across methods"）。
    正确做法是取每个子序列**自身**流的第一个 32-bit 字 ``generate_state(1)[0]``，
    对固定父种子稳定可复现，且各子流互不相同。

    仅当 NumPy 的 SeedSequence 不可用时（ImportError/AttributeError）才回退到
    ``secrets``（此时互异但不可复现）。
    """
    names = list(method_names)
    try:
        import numpy as np

        children = np.random.SeedSequence(base_seed).spawn(len(names))
        return {
            name: int(child.generate_state(1)[0])
            for name, child in zip(names, children)
        }
    except (ImportError, AttributeError):
        import secrets

        return {name: secrets.randbelow(2**31) for name in names}


def _register_temp_directory(temp_dir: Path) -> None:
    """注册临时目录以便在退出时清理（线程安全）"""
    with _temp_dirs_lock:
        _temp_directories.add(temp_dir)
    # 显式设置临时目录权限为 0o700，确保仅当前用户可访问
    try:
        os.chmod(str(temp_dir), 0o700)
    except OSError as e:
        # 审阅项 C-8：权限收紧失败不致命，但必须可见（安全属性没有真正生效）
        get_logger().warning(
            f"未能将临时目录权限收紧为 0o700（{temp_dir}）：{e}。"
            "该目录可能对其他用户可读。"
        )


def _cleanup_temp_directories() -> None:
    """清理所有注册的临时目录（通过atexit调用）"""
    import shutil

    with _temp_dirs_lock:
        dirs_to_clean = list(_temp_directories)
        _temp_directories.clear()
    for temp_dir in dirs_to_clean:
        try:
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
        except OSError as e:
            # atexit 阶段日志器可能已被拆解，因此直接写"原始" stderr；
            # 若连 stderr 都已关闭，则确实再无通道可报告（解释器正在退出）。
            stream = sys.__stderr__ or sys.stderr
            try:
                stream.write(f"[phylodater] 临时目录清理失败（{temp_dir}）：{e}\n")
            except (OSError, ValueError):
                pass


# ==================== B-23：运行输入指纹 ====================
#
# 审阅项 B-23 [已证]：``--resume`` 的跳过判定此前**只**看"方法名 + 上次是否
# 成功"，完全不含输入信息。改一行校准 YAML、换一棵树、调一个先验，在同一输出
# 目录再跑一次 ``--resume``，软件会把上一次的年龄当作本次结果返回。
# ``checkpoint.record_input_hash()`` 虽然写好了却零调用点，读取比对用的
# ``get_input_hash``/``input_matches`` 更是根本不存在。
#
# 下面这层把 pipeline 侧接上：每次运行都算一个"运行指纹"，续跑时先比对指纹，
# 不一致（或旧检查点里压根没有指纹、无法判定）就**拒绝**复用陈旧结果。
#
# 与 infrastructure/checkpoint.py 的接口约定（pipeline 侧全部 getattr 守卫，
# 缺任一件都能降级工作）：
#   1. ``record_input_hash(key, hash_value)``            —— 已存在，必需
#      ``key`` 既接受真实文件路径，也接受逻辑键
#      ``RUN_FINGERPRINT_KEY``（内部以 ``str(key)`` 存储，因此字符串可用）。
#   2. ``record_run_parameters(params: dict)``           —— 已存在，可选
#      pipeline 会把 ``{"run_fingerprint": ..., "fingerprint_schema": ...}``
#      同时写进 run_parameters，便于人工核对。
#   3. ``input_matches(current_hash: str) -> bool``      —— **待 checkpoint 作者实现**
#      语义：本次输入指纹是否与检查点里存的一致；没有存过 → False。
#   4. ``compute_input_hash(inputs: Mapping) -> str``    —— **待 checkpoint 作者实现**
#      语义：由 pipeline 给的规范化输入字典算出指纹（哈希即字符串）。
#      缺省时 pipeline 用本地 sha256(json) 计算，结果稳定但不与之互通。
#   5. ``get_input_hash(key) -> Optional[str]``          —— **待 checkpoint 作者实现**
#      读取指纹；缺省时退化为 ``check_input_changed()``（已存在）。

#: 运行指纹在检查点 ``input_hashes`` 中占用的逻辑键（不是真实路径）
RUN_FINGERPRINT_KEY = "phylodater:run-fingerprint"
#: 指纹口径版本：改动 ``_build_run_inputs`` 的字段构成时必须递增，
#: 否则新旧指纹会被误判为"输入变了"（安全）或"输入没变"（危险）。
FINGERPRINT_SCHEMA = "phylodater-run-fingerprint-v1"

#: 与 ``DeepValidator.ULTRAMETRIC_TOLERANT_METHODS`` 保持一致（B-18 分档）：
#: 这些方法设计上就允许非超度量输入，超度量性检查只应报 info。
_ULTRAMETRIC_TOLERANT_METHODS = frozenset(
    {"treepl", "pathd8", "r8s", "pyr8s", "mdcat", "md-cat", "wlogdate"}
)


def _sha256_of_file(path: Path) -> Optional[str]:
    """计算文件内容的 sha256（读不到时返回 None，不抛异常）。"""
    try:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError as e:
        get_logger().warning(f"无法读取 {path} 以计算输入指纹：{e}")
        return None


def _stable_repr(obj: Any) -> str:
    """把任意输入对象转成**与字典序无关**的稳定字符串，用于指纹计算。"""
    if is_dataclass(obj) and not isinstance(obj, type):
        try:
            from dataclasses import asdict

            return _json_dumps(asdict(obj), sort_keys=True, default=repr)
        except Exception as e:  # 非常规字段（循环引用等）→ 退回 repr
            get_logger().debug(f"dataclass 结构化失败，改用 repr：{e}")
    try:
        return _json_dumps(obj, sort_keys=True, default=repr)
    except (TypeError, ValueError) as e:
        get_logger().debug(f"JSON 规范化失败，改用 repr：{e}")
        return repr(obj)


# 注册退出清理函数
atexit.register(_cleanup_temp_directories)


def detect_system_resources() -> Dict[str, Any]:
    """
    探测系统可用资源

    Returns:
        包含 CPU 核心数、内存信息的字典
    """
    try:
        import psutil

        return {
            "cpu_count": psutil.cpu_count(logical=True),
            "cpu_count_physical": psutil.cpu_count(logical=False),
            "memory_total_mb": psutil.virtual_memory().total / (1024 * 1024),
            "memory_available_mb": psutil.virtual_memory().available / (1024 * 1024),
        }
    except ImportError:
        return {
            "cpu_count": os.cpu_count() or 1,
            "cpu_count_physical": None,
            "memory_total_mb": None,
            "memory_available_mb": None,
        }


def _ignored_env_value(name: str, value: str, reason: str) -> None:
    """统一披露"环境变量被忽略"（审阅项 C-8）。

    环境变量写错却静默回落默认值，等于用户的配置凭空消失而没有任何痕迹；
    这类失败与报告在多处批评的"降级只写日志都不够"是同族，这里至少保证
    在日志里能看到被丢弃的原值。
    """
    get_logger().warning(f"忽略非法环境变量 {name}={value!r}（{reason}），使用默认值")


def generate_global_seed() -> int:
    """
    生成全局随机种子

    优先使用环境变量 PHYLODATER_SEED，
    否则生成随机种子
    """
    env_seed = os.environ.get("PHYLODATER_SEED")
    if env_seed is not None:
        try:
            seed = int(env_seed)
            return seed
        except ValueError:
            _ignored_env_value("PHYLODATER_SEED", env_seed, "不是整数")
    import secrets

    return secrets.randbelow(2**31)


def _env_positive_int(name: str, default: int) -> int:
    """读取"必须为正整数"的环境变量，非法值告警后回落默认。"""
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        value = int(raw)
    except ValueError:
        _ignored_env_value(name, raw, "不是整数")
        return default
    if value <= 0:
        _ignored_env_value(name, raw, "必须 > 0")
        return default
    return value


def get_env_threads(default: int = 4) -> int:
    """从环境变量获取线程数"""
    return _env_positive_int("PHYLODATER_THREADS", default)


def get_env_burnin(default: int = 20000) -> int:
    """从环境变量获取MCMC burnin代数"""
    return _env_positive_int("PHYLODATER_BURNIN", default)


def get_env_nsample(default: int = 50000) -> int:
    """从环境变量获取MCMC采样数"""
    return _env_positive_int("PHYLODATER_NSAMPLE", default)


def get_env_timeout(default: int = 7200) -> int:
    """从环境变量获取执行超时时间（秒）"""
    return _env_positive_int("PHYLODATER_TIMEOUT", default)


def sanitize_path(path: str, base_dir: Optional[Path] = None) -> str:
    """
    脱敏路径，将绝对路径转换为相对路径

    Args:
        path: 原始路径
        base_dir: 基准目录，用于计算相对路径

    Returns:
        脱敏后的路径（相对路径或掩码）
    """
    path_obj = Path(path)
    if base_dir is None:
        base_dir = Path.cwd()

    try:
        rel_path = path_obj.resolve().relative_to(base_dir.resolve())
        return str(rel_path)
    except ValueError:
        return f"***/{path_obj.name}"


@dataclass
class PipelineConfig:
    """流水线配置

    支持通过环境变量覆盖默认配置：
    - PHYLODATER_THREADS: 并行线程数
    - PHYLODATER_SEED: 全局随机种子
    - PHYLODATER_BURNIN: MCMC burnin代数
    - PHYLODATER_NSAMPLE: MCMC采样数
    - PHYLODATER_TIMEOUT: 执行超时时间（秒）
    """

    methods: List[str] = field(default_factory=lambda: ["mcmctree", "treepl", "pathd8"])
    output_dir: Path = Path("./phylodater_output")
    threads: int = field(default_factory=lambda: get_env_threads(4))
    skip_failed: bool = True
    save_intermediates: bool = True
    dry_run: bool = False
    seed: Optional[int] = None
    burnin: int = field(default_factory=lambda: get_env_burnin(20000))
    nsample: int = field(default_factory=lambda: get_env_nsample(50000))
    timeout: int = field(default_factory=lambda: get_env_timeout(7200))
    qc_plugins: List[str] = field(default_factory=list)
    enable_checkpoint: bool = True
    software_paths: Optional[SoftwarePaths] = None
    taxonomy_file: Optional[Path] = None
    taxonomy_delimiter_mode: Optional[str] = None
    taxonomy_source_priority: Optional[str] = None
    taxonomy_levels: Optional[dict] = None
    taxonomy_table_sep: Optional[str] = None
    taxonomy_sep: Optional[str] = None
    reroot_strategy: str = "none"  # 'none', 'outgroup', 'midpoint', 'mad'
    outgroup_name: Optional[str] = None
    paml_version: Optional[str] = None
    low_memory: bool = False
    resume_mode: bool = False  # --resume: 跳过已完成的检查点步骤
    strip_annotations: bool = False  # --strip-annotations: 输出时丢弃树注释

    def __post_init__(self) -> None:
        """验证参数范围"""
        if not self.methods:
            raise ValueError("methods list cannot be empty")
        if self.threads < 1:
            raise ValueError(f"threads must be >= 1, got {self.threads}")
        if self.seed is not None and self.seed < 0:
            raise ValueError(f"seed must be non-negative, got {self.seed}")
        if self.burnin < 0:
            raise ValueError(f"burnin must be >= 0, got {self.burnin}")
        if self.nsample < 1:
            raise ValueError(f"nsample must be >= 1, got {self.nsample}")
        if self.timeout < 1:
            raise ValueError(f"timeout must be >= 1, got {self.timeout}")
        if self.reroot_strategy not in ["none", "outgroup", "midpoint", "mad"]:
            raise ValueError(
                f"reroot_strategy must be 'none', 'outgroup', 'midpoint', or 'mad', got '{self.reroot_strategy}'"
            )


class DatingPipeline:
    """
    定年流水线

    整合所有步骤的高阶 API：
    1. 输入验证
    2. 树定根（可选）
    3. 校准点解析
    4. 并行执行多种定年方法
    5. 结果汇总
    """

    def __init__(
        self, config: PipelineConfig, tool_config: Optional[ToolConfig] = None
    ) -> None:
        self.config = config
        self.logger = get_logger()
        self.results: Dict[str, DatingResult] = {}
        #: 运行期发生的**可见降级**记录（审阅项 B-26 同族要求：降级可以发生，
        #: 但必须被记下来并在日志里出现，而不是伪装成一次干净的运行）。
        self.degradations: List[str] = []
        self._tree: Optional[PhylogeneticTree] = None
        #: B-23：``run()`` 收到的原始校准点列表（解析前），用于运行指纹。
        self._source_calibrations: List[CalibrationPoint] = []
        #: B-23：文件哈希的实例级缓存（避免同一运行里重复读大比对文件）
        self._file_hash_cache: Dict[str, Optional[str]] = {}
        self._tool_config = (
            tool_config if tool_config is not None else Configuration().get_config()
        )
        self._software_paths = self.config.software_paths or SoftwarePaths()
        self._metadata_manager: Optional[RuntimeMetadataManager] = None
        self._qc_chain: Optional[QCPluginChain] = None
        self._checkpoint_manager: Optional[PipelineCheckpointManager] = None

        # 应用 PipelineConfig 中的 PAML 版本到 tool_config
        if self.config.paml_version is not None:
            self._tool_config.mcmctree.paml_version = self.config.paml_version
            self.logger.info(
                f"PAML version set to {self.config.paml_version} from pipeline config"
            )

        # 同步全局线程数和随机种子到 tool_config.common，确保适配器内部使用一致
        if self.config.threads is not None:
            self._tool_config.common.nthreads = self.config.threads
        if self.config.seed is not None:
            self._tool_config.common.seed = self.config.seed

        # 同步 burnin/nsample/timeout 到 tool_config：
        #   仅当 PipelineConfig 显式提供了非 None 值(来自 --method-args 透传)时才同步，
        #   避免 PipelineConfig 默认值(20000/50000)覆盖用户已通过 --method-args
        #   在 tool_config.mcmctree 中设定的值(F15)。
        if self.config.burnin is not None:
            self._tool_config.mcmctree.burnin = self.config.burnin
        if self.config.nsample is not None:
            self._tool_config.mcmctree.nsample = self.config.nsample
        if self.config.timeout is not None:
            self._tool_config.common.timeout = self.config.timeout

        # §16.1 - low_memory 模式：减少资源使用
        if self.config.low_memory:
            self._tool_config.common.nthreads = 1
            self._tool_config.mcmctree.burnin = min(
                self._tool_config.mcmctree.burnin, 5000
            )
            self._tool_config.mcmctree.nsample = min(
                self._tool_config.mcmctree.nsample, 10000
            )
            self.logger.info("Low-memory mode: threads=1, reduced MCMC parameters")

        if self.config.seed is None:
            # B-23：区分"用户显式 --seed"与"每次运行自动生成的种子"。
            # 自动种子必然每次不同，把它算进运行指纹会让 ``--resume`` 永远
            # 判定"输入已变化"、从而完全失去续跑能力。
            self._user_seed = None
            self.config.seed = generate_global_seed()
        else:
            self._user_seed = self.config.seed
        # 将最终使用的全局种子同步到 tool_config.common，确保各适配器在默认
        # 运行中也能读取到一致种子（无论用户是否显式传入 --seed）。
        self._tool_config.common.seed = self.config.seed

        # §11.1 - 检查点管理
        # enable_checkpoint 控制是否创建检查点管理器
        # resume_mode 控制是否跳过已完成的步骤（--resume 标志）
        self._resume_mode = getattr(config, "resume_mode", False)

        if self.config.enable_checkpoint:
            self._checkpoint_manager = PipelineCheckpointManager(config.output_dir)

        self._init_qc_chain()
        self._init_random_state()

    def _init_qc_chain(self) -> None:
        """初始化质控插件链

        从配置和全局插件注册表加载质控插件。
        插件可支持长枝吸引过滤、序列质量检查等功能。
        """
        self._qc_chain = QCPluginChain()

        if not self.config.qc_plugins:
            registry = get_plugin_registry()
            enabled_plugins = registry.list_plugins(plugin_type=None)
            for plugin_meta in enabled_plugins:
                if plugin_meta.plugin_type.value == "qc" and registry.is_enabled(
                    plugin_meta.name
                ):
                    plugin = registry.get_plugin(plugin_meta.name)
                    if plugin and isinstance(plugin.instance, QCPlugin):
                        self._qc_chain.add(plugin.instance)
                        self.logger.debug(f"Loaded QC plugin: {plugin_meta.name}")
                    elif plugin and plugin.instance is not None:
                        # 注册表按 plugin_type 筛过一轮，实例却不是 QCPlugin：
                        # 这是注册侧的类型不符，必须说出来，不能塞进链里。
                        self.logger.warning(
                            f"QC plugin '{plugin_meta.name}' is registered as a QC "
                            f"plugin but its instance is "
                            f"{type(plugin.instance).__name__}; it was NOT added to "
                            "the QC chain."
                        )
        else:
            registry = get_plugin_registry()
            for plugin_name in self.config.qc_plugins:
                plugin = registry.get_plugin(plugin_name)
                if plugin and isinstance(plugin.instance, QCPlugin):
                    self._qc_chain.add(plugin.instance)
                    self.logger.debug(f"Loaded QC plugin: {plugin_name}")
                elif plugin and plugin.instance is not None:
                    self.logger.warning(
                        f"QC plugin '{plugin_name}' resolved to a "
                        f"{type(plugin.instance).__name__}, not a QCPlugin; it was "
                        "NOT added to the QC chain."
                    )
                else:
                    self.logger.warning(
                        f"QC plugin not found or not enabled: {plugin_name}"
                    )

        if self._qc_chain.get_plugins():
            self.logger.info(f"QC plugins loaded: {len(self._qc_chain.get_plugins())}")
        else:
            self.logger.debug("No QC plugins loaded")

    def _init_metadata_manager(self) -> None:
        """初始化元数据管理器"""
        from .. import __version__

        self._metadata_manager = RuntimeMetadataManager(
            output_dir=self.config.output_dir,
            pipeline="phylodater",
            version=__version__,
        )
        self._metadata_manager.set_seed(self.config.seed)
        self._metadata_manager.collect_system_resources()

    def _init_random_state(self) -> None:
        """初始化随机状态并记录系统资源

        统一在此处设置随机种子，避免重复初始化。
        """
        # 先设置随机种子（确保后续资源探测不影响种子状态）
        import random

        random.seed(self.config.seed)

        import numpy as np

        try:
            np.random.seed(self.config.seed)
        except AttributeError as e:
            # NumPy 2 移除了全局 np.random.seed：可复现性依赖下面 ParallelPipeline
            # 里的 Generator/SeedSequence 派生路径，此处降级必须留痕（C-8）。
            self.logger.debug(f"全局 np.random.seed 不可用（{e}），跳过播种")

        # 记录系统资源信息
        resources = detect_system_resources()

        self.logger.section("System Resources & Reproducibility")
        self.logger.kv("CPU Cores (Logical)", resources.get("cpu_count", "N/A"))
        if resources.get("cpu_count_physical"):
            self.logger.kv("CPU Cores (Physical)", resources.get("cpu_count_physical"))
        self.logger.kv(
            "Total Memory (MB)",
            (
                f"{resources.get('memory_total_mb', 'N/A'):.1f}"
                if resources.get("memory_total_mb")
                else "N/A"
            ),
        )
        self.logger.kv(
            "Available Memory (MB)",
            (
                f"{resources.get('memory_available_mb', 'N/A'):.1f}"
                if resources.get("memory_available_mb")
                else "N/A"
            ),
        )
        self.logger.kv("Global Random Seed", self.config.seed)
        self.logger.kv(
            "Seed Source",
            (
                "PHYLODATER_SEED env var"
                if os.environ.get("PHYLODATER_SEED")
                else "Auto-generated"
            ),
        )

        self.logger.info(
            "Recommended threads: min(available_cores, num_methods) = "
            f"{min(resources.get('cpu_count', 4), len(self.config.methods))}"
        )

    def run(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path],
        calibrations: List[CalibrationPoint],
    ) -> Dict[str, DatingResult]:
        """
        运行完整流水线

        Args:
            tree: 系统发育树
            alignment_path: 比对文件路径。可以是 None：树定年类方法（r8s /
                treePL / PATHd8 …）只要带分支长度的树，根本不读比对；
                CLI 在 ``--sequence`` 缺省时会直接传 None 进来。
            calibrations: 校准点列表

        Returns:
            方法名到结果的映射
        """
        if self.config.dry_run:
            return self._dry_run(tree, alignment_path, calibrations)

        self._tree = tree
        self._init_metadata_manager()
        #: B-23：留存"用户交给本次运行的原始输入"，供运行指纹使用。
        #: ``_execute_methods`` 拿到的是解析后的校准点，单靠它无法区分
        #: "用户改了 YAML 但被某层忽略" 与 "用户什么都没改"。
        self._source_calibrations = calibrations

        start_time = time.time()
        # 收尾处既要接 KeyboardInterrupt 也要接 Exception，所以它不是“一种异常”。
        pipeline_error: Optional[BaseException] = None

        try:
            self.logger.section("PhyloDater Pipeline")

            step_times = {}

            if self._qc_chain and self._qc_chain.get_plugins():
                step_start = time.time()
                self.logger.step(0, "QC Processing")
                tree, alignment_path = self._run_qc_plugins(tree, alignment_path)
                step_times["qc_processing"] = time.time() - step_start

            if self.config.reroot_strategy != "none":
                tree = self._reroot_tree(tree)

            # 更新 self._tree 确保报告生成器使用最终的树
            self._tree = tree

            step_start = time.time()
            self.logger.step(1, "Input Validation")
            self._validate_inputs(tree, alignment_path)
            step_times["input_validation"] = time.time() - step_start

            step_start = time.time()
            self.logger.step(2, "Calibration Resolution")
            resolved_calibrations = self._resolve_calibrations(tree, calibrations)
            step_times["calibration_resolution"] = time.time() - step_start

            step_start = time.time()
            self.logger.step(3, "Dating Analysis")
            self._execute_methods(tree, alignment_path, resolved_calibrations)
            step_times["dating_analysis"] = time.time() - step_start

            step_start = time.time()
            self.logger.step(4, "Report Generation")
            self._generate_reports()
            step_times["report_generation"] = time.time() - step_start

            elapsed = time.time() - start_time
            self._log_timing_summary(step_times)

            self.logger.section(f"Pipeline Completed in {elapsed:.1f}s")

        except KeyboardInterrupt as e:
            pipeline_error = e
            self.logger.warning("Pipeline interrupted by user")
        except Exception as e:
            pipeline_error = e
            self.logger.error(f"Pipeline failed: {e}")
        finally:
            # 运行期发生的所有可见降级，无论成败都要在收尾处复述一遍
            self._report_degradations()
            # 确保元数据始终被保存（无论成功或失败）
            if self._metadata_manager:
                if pipeline_error is None:
                    status = RunStatus.COMPLETED
                elif isinstance(pipeline_error, KeyboardInterrupt):
                    status = RunStatus.INTERRUPTED
                else:
                    status = RunStatus.FAILED
                self._metadata_manager.finalize(status)
                try:
                    metadata_path = self._metadata_manager.save()
                    self.logger.info(f"Runtime metadata saved to {metadata_path}")
                except Exception as meta_err:
                    self.logger.warning(f"Failed to save runtime metadata: {meta_err}")

        if pipeline_error is not None:
            raise pipeline_error

        return self.results

    def _run_qc_plugins(
        self, tree: PhylogeneticTree, alignment_path: Optional[Path]
    ) -> tuple:
        """运行质控插件链处理树和比对

        Args:
            tree: 输入的系统发育树
            alignment_path: 比对文件路径

        Returns:
            (处理后的树, 处理后的比对路径)
        """
        if not self._qc_chain or not self._qc_chain.get_plugins():
            return tree, alignment_path

        context = {
            "alignment_path": str(alignment_path),
            "output_dir": str(self.config.output_dir),
        }

        try:
            processed_tree = self._qc_chain.process_tree(tree, context)
            if processed_tree is not None:
                tree = processed_tree
                self.logger.success("QC plugins processed tree")
            else:
                self.logger.warning("QC plugin returned None for tree, using original")

            processed_alignment = self._qc_chain.process_alignment(
                alignment_path, context
            )
            if processed_alignment != alignment_path:
                alignment_path = processed_alignment
                self.logger.success("QC plugins processed alignment")

        except Exception as e:
            self.logger.warning(f"QC processing failed: {e}")
            self.logger.warning("Continuing without QC processing")
            self._record_degradation(f"质控插件链执行失败（{e}），本次运行未经 QC 处理")

        return tree, alignment_path

    def _reroot_tree(self, tree: PhylogeneticTree) -> PhylogeneticTree:
        """根据配置对树进行定根"""
        from ..services import RootingMethod, TreeRootingService

        strategy = self.config.reroot_strategy
        if strategy == "none":
            return tree

        rerooter = TreeRootingService()

        strategy_map = {
            "outgroup": RootingMethod.OUTGROUP,
            "midpoint": RootingMethod.MIDPOINT,
            "mad": RootingMethod.MAD,
        }

        if strategy not in strategy_map:
            raise PhyloDaterError(f"Unknown reroot strategy: {strategy}")

        method = strategy_map[strategy]
        outgroup_taxa = None
        if method == RootingMethod.OUTGROUP:
            if not self.config.outgroup_name:
                raise PhyloDaterError("--outgroup required for outgroup rooting")
            outgroup_taxa = [self.config.outgroup_name]

        result = rerooter.root_tree(tree, method=method, outgroup_taxa=outgroup_taxa)
        tree = result.rooted_tree

        self.logger.info(f"Rerooted tree with strategy: {strategy}")
        self._tree = tree
        return tree

    def _dry_run(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path],
        calibrations: List[CalibrationPoint],
    ) -> Dict[str, DatingResult]:
        """
        干运行模式：不实际执行，仅验证配置和依赖
        """
        self.logger.section("Dry Run Mode")

        self.logger.step(1, "Configuration Validation")
        self.logger.info(f"Output directory: {self.config.output_dir}")
        self.logger.info(f"Methods: {', '.join(self.config.methods)}")
        self.logger.info(f"Threads: {self.config.threads}")
        self.logger.info(f"Global seed: {self.config.seed}")

        self.logger.step(2, "Input Validation (Dry)")
        try:
            validator = TreeValidator()
            report = validator.full_validation(
                tree,
                alignment_path,
                method=self._ultrametric_check_method(),
            )
            if report.is_valid:
                self.logger.success("Input validation passed")
            else:
                self.logger.warning(f"Input validation issues: {report.errors}")
        except Exception as e:
            self.logger.warning(f"Input validation skipped: {e}")

        self.logger.step(3, "Dependency Check (Dry)")
        for method_name in self.config.methods:
            if DatingMethodRegistry.is_registered(method_name):
                adapter = DatingMethodRegistry.create(
                    method_name,
                    self._tool_config,
                    self.config.output_dir,
                    self._software_paths,
                )
                available = adapter.validate_environment()
                status = "Available" if available else "Not found"
                self.logger.info(f"  {method_name}: {status}")
            else:
                self.logger.warning(f"  {method_name}: Not registered")

        self.logger.step(4, "Calibration Check (Dry)")
        self.logger.info(f"  Loaded {len(calibrations)} calibration points")

        # 真实解析校验点（MRCA 定位、名称匹配），以便在 dry-run 阶段就暴露
        # 类群名不匹配等问题，而非留到真实运行。解析失败时仅告警（不终止），
        # 让依赖检查等其他校验仍继续执行。
        resolved_calibrations: List[CalibrationPoint] = calibrations
        try:
            resolved_calibrations = self._resolve_calibrations(tree, calibrations)
            tips = set(tree.tip_names)
            for cal in resolved_calibrations:
                pair = cal.mrca_leaf_pair or ()
                missing = [t for t in pair if t not in tips]
                if missing:
                    self.logger.warning(
                        f"  Calibration '{cal.name}': mrca_leaf_pair 中的类群名"
                        f"未在树中找到 -> {missing}"
                    )
                elif pair:
                    self.logger.info(
                        f"  Calibration '{cal.name}': mrca_pair "
                        f"({pair[0][:20]}..., {pair[1][:20]}...) resolved OK"
                    )
        except Exception as e:
            self.logger.warning(f"  Calibration resolution skipped: {e}")

        # 方法特异的约束充分性校验：在 dry-run 阶段预警真实运行会失败的
        # 确定性约束（如 pathd8 必须至少有一个 fixed-age 校准点）。
        self._check_method_specific_constraints(
            self.config.methods, resolved_calibrations
        )

        self.logger.section("Dry Run Completed - No actual analysis performed")
        return {}

    # 方法特异的约束充分性校验（dry-run / 真实运行前置均可复用）
    def _check_method_specific_constraints(
        self, methods: List[str], calibrations: List[CalibrationPoint]
    ) -> None:
        """对指定方法列表做约束充分性校验,不满足时发出 WARNING。

        目前实现:
        - pathd8: 要求至少一个固定年龄 (fixed / fixage) 校准点。
        """
        from ..models.constraints import FixedAgeConstraint

        method_set = set(methods)
        if "pathd8" in method_set:
            has_fixed = any(
                isinstance(cal.age_constraint, FixedAgeConstraint)
                for cal in calibrations
            )
            if not has_fixed:
                self.logger.warning(
                    "  pathd8 requires at least one fixed-age (fixage) calibration, "
                    "but none was found. The real run will fail with "
                    "'PATHd8 requires at least one fixed age (fixage) calibration'. "
                    "Please add a calibration with 'type: fixed' / 'age: <value>'."
                )

    def _log_timing_summary(self, step_times: Dict[str, float]) -> None:
        """记录各步骤耗时"""
        self.logger.info("Timing Summary:")
        for step, duration in step_times.items():
            self.logger.kv(f"  {step}", f"{duration:.1f}s")

    def _record_degradation(self, message: str) -> None:
        """记录一次"流程继续了、但没完全按用户意图执行"的降级（B-26 同族）。

        降级本身允许发生，但必须：(1) 以 ERROR 级别出现一次；(2) 累积在
        ``self.degradations`` 里，供收尾时复述、也供库调用方直接检查。
        """
        self.degradations.append(message)
        self.logger.error(f"运行降级: {message}")

    def _report_degradations(self) -> None:
        """收尾时把本次运行的全部可见降级复述一遍。"""
        if not self.degradations:
            return
        self.logger.section(f"Disclosed Degradations ({len(self.degradations)})")
        for index, message in enumerate(self.degradations, start=1):
            self.logger.warning(f"  [{index}] {message}")

    # ==================== B-23：运行输入指纹 ====================

    def _hash_input_file(self, path: Any) -> Optional[str]:
        """带实例级缓存的文件哈希。

        一次运行里指纹计算与指纹登记都要读比对文件；比对常有上百 MB，
        读两遍纯属浪费。缓存只活在单个 pipeline 实例内，不存在跨运行串味。
        """
        if path is None:
            return None
        key = str(path)
        cache = self._file_hash_cache
        if key not in cache:
            cache[key] = _sha256_of_file(Path(path))
        return cache[key]

    def _build_run_inputs(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path],
        calibrations: List[CalibrationPoint],
    ) -> Dict[str, Any]:
        """收集"决定本次运行结果的全部输入"，作为运行指纹的原料。

        收录的是**会改变科学结果**的量：树的拓扑与枝长、比对文件内容、
        校准点（原始 + 解析后两份）、生效后的工具配置、方法列表、定根策略、
        分类学表内容、外部软件可执行文件路径与随机种子。
        刻意**不**收录 ``output_dir``——续跑本来就发生在同一目录。
        """
        from dataclasses import asdict

        tree_text = getattr(tree, "newick", None)
        if tree_text is None:
            tree_text = repr(tree)

        taxonomy_hash = None
        if self.config.taxonomy_file:
            taxonomy_hash = self._hash_input_file(self.config.taxonomy_file)

        software_paths: Dict[str, Any] = {}
        try:
            software_paths = asdict(self._software_paths)
        except TypeError as e:  # 非 dataclass 的替身对象
            self.logger.debug(f"software_paths 结构化失败，改用 repr：{e}")
            software_paths = {"repr": repr(self._software_paths)}

        tool_config: Any
        try:
            tool_config = asdict(self._tool_config)
            # ``__init__`` 会把最终生效的种子（可能是每次运行都不同的自动种子）
            # 同步进 ``common.seed``。运行指纹只认**用户显式指定**的种子：
            # 否则 ``--resume`` 会因为一个与输入无关的随机量而永远无法复用。
            if isinstance(tool_config.get("common"), dict):
                common = dict(tool_config["common"])
                common["seed"] = getattr(self, "_user_seed", None)
                tool_config["common"] = common
        except TypeError as e:
            self.logger.debug(f"tool_config 结构化失败，改用 repr：{e}")
            tool_config = repr(self._tool_config)

        return {
            "fingerprint_schema": FINGERPRINT_SCHEMA,
            "methods": list(self.config.methods),
            "reroot_strategy": self.config.reroot_strategy,
            "outgroup_name": self.config.outgroup_name,
            # 只登记用户显式指定的种子（见 ``__init__`` 处的说明）
            "seed": getattr(self, "_user_seed", None),
            "strip_annotations": bool(getattr(self.config, "strip_annotations", False)),
            "tree_newick": tree_text,
            "tree_num_tips": getattr(tree, "num_tips", None),
            "alignment_file": getattr(alignment_path, "name", str(alignment_path)),
            "alignment_sha256": self._hash_input_file(alignment_path),
            "taxonomy_file": getattr(self.config.taxonomy_file, "name", None),
            "taxonomy_sha256": taxonomy_hash,
            # 解析后的校准点（含 resolved_taxa / mrca_leaf_pair）
            "resolved_calibrations": [_stable_repr(c) for c in calibrations],
            # 用户交给本次运行的原始校准点：解析层丢弃或改写了什么，
            # 在这里也留一份底，避免"改了但被上层忽略"→ 指纹不变 → 复用陈旧结果。
            "source_calibrations": [
                _stable_repr(c) for c in getattr(self, "_source_calibrations", [])
            ],
            "tool_config": tool_config,
            "software_paths": software_paths,
        }

    def _compute_run_fingerprint(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path],
        calibrations: List[CalibrationPoint],
    ) -> str:
        """算出本次运行的输入指纹。

        若 ``checkpoint.compute_input_hash(inputs)`` 已由 checkpoint 侧提供
        （B-23 接口约定的第 4 件），由它统一 canonicalization；否则本地用
        ``sha256(json)`` 计算——两者都能独立完成"输入变了没有"的判断，只是
        字符串互不兼容，因此切换实现时旧检查点会被判为"无法核验"并重跑。
        """
        inputs = self._build_run_inputs(tree, alignment_path, calibrations)

        hasher = getattr(self._checkpoint_manager, "compute_input_hash", None)
        if callable(hasher):
            try:
                value = hasher(inputs)
                if isinstance(value, str) and value:
                    return value
                self.logger.warning(
                    "checkpoint.compute_input_hash() 未返回字符串指纹，"
                    "改用 pipeline 本地哈希"
                )
            except TypeError as e:
                # 签名与约定不符（例如要求 (file_path) 而非 (mapping)）
                self.logger.warning(
                    f"checkpoint.compute_input_hash(inputs) 调用失败（{e}），"
                    "改用 pipeline 本地哈希"
                )

        canonical = _json_dumps(inputs, sort_keys=True, default=repr)
        return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    def _store_run_fingerprint(
        self,
        fingerprint: str,
        alignment_path: Optional[Path] = None,
    ) -> None:
        """把运行指纹（以及各输入文件自身的哈希）写回检查点。

        ``record_input_hash`` 此前是零调用点的死代码（B-23），这里既登记聚合
        指纹，也逐个登记真实输入文件的哈希，使"哪个输入变了"可被单独回答。
        """
        checkpoint = self._checkpoint_manager
        if checkpoint is None:
            return

        recorder = getattr(checkpoint, "record_input_hash", None)
        if callable(recorder):
            entries: List[tuple] = [(RUN_FINGERPRINT_KEY, fingerprint)]
            for path in filter(None, [alignment_path, self.config.taxonomy_file]):
                digest = self._hash_input_file(path)
                if digest:
                    entries.append((Path(path), digest))
            for key, value in entries:
                try:
                    recorder(key, value)
                except TypeError as e:
                    self.logger.warning(
                        f"checkpoint.record_input_hash({key}, …) 签名不符：{e}"
                    )
                    break
                except OSError as e:
                    self.logger.warning(f"写入输入指纹失败（{key}）：{e}")
        else:
            self.logger.warning(
                "检查点管理器不提供 record_input_hash()，续跑将无法核验输入是否变化"
            )

        params_recorder = getattr(checkpoint, "record_run_parameters", None)
        if callable(params_recorder):
            try:
                params_recorder(
                    {
                        "run_fingerprint": fingerprint,
                        "fingerprint_schema": FINGERPRINT_SCHEMA,
                    }
                )
            except TypeError as e:
                self.logger.debug(f"record_run_parameters 调用失败：{e}")

    def _run_inputs_are_current(self, fingerprint: str) -> Optional[bool]:
        """本次输入是否与检查点里记录的一致。

        Returns:
            ``True``  —— 一致，可以复用缓存结果；
            ``False`` —— 不一致（或检查点里没有指纹可核对），必须重跑；
            ``None``  —— 检查点侧完全没有可用的读取接口，无法判定。
            调用方对 ``None`` 与 ``False`` 一视同仁地拒绝复用：宁可用"必须重跑"
            换回正确性，也不把上一次输入的年龄当成本次的结果（B-23 的成因）。
        """
        checkpoint = self._checkpoint_manager
        if checkpoint is None:
            return True

        # 1) 首选 checkpoint 侧的统一判定（B-23 接口约定的第 3 件）
        matcher = getattr(checkpoint, "input_matches", None)
        if callable(matcher):
            try:
                return bool(matcher(fingerprint))
            except TypeError as e:
                self.logger.debug(f"checkpoint.input_matches() 签名不符：{e}")

        # 2) 退化为"读出来自己比"
        getter = getattr(checkpoint, "get_input_hash", None)
        if callable(getter):
            # stored 来自第三方/可选后端的 getter，标注成 object 而不是 Any：
            # 本函数声明返回 Optional[bool]，``stored == fingerprint`` 必须是个 bool。
            stored: object = None
            try:
                stored = getter(RUN_FINGERPRINT_KEY)
            except TypeError:
                stored = None
            else:
                if stored is None:
                    return False  # 没有存过 —— 无从确认，按"已变更"处理
                # 两侧都可能是非字符串（老检查点存过 dict/None），== 仍返回 bool
                return bool(stored == fingerprint)

        # 3) 退化到已有的"是否变更"探针（absence 同样返回 True = 已变更）
        differ = getattr(checkpoint, "check_input_changed", None)
        if callable(differ):
            try:
                return not bool(differ(RUN_FINGERPRINT_KEY, fingerprint))
            except TypeError as e:
                self.logger.debug(f"checkpoint.check_input_changed() 签名不符：{e}")

        # 4) 最后读内部映射（同包内的受控降级，只读不写）
        stored_map = getattr(checkpoint, "_input_hashes", None)
        if isinstance(stored_map, dict):
            stored = stored_map.get(RUN_FINGERPRINT_KEY)
            if stored is None:
                return False
            return bool(stored == fingerprint)

        return None

    def _plan_method_execution(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path],
        calibrations: List[CalibrationPoint],
    ) -> List[str]:
        """决定本次要执行哪些方法；``--resume`` 时先核验输入指纹（B-23）。

        旧行为：跳过判定的全部依据是"方法名 + 上次是否成功"，于是换树、改校准、
        调先验之后在同一输出目录 ``--resume``，会把上一次输入的年龄当成本次结果
        放回 ``self.results``。

        新行为：只有当**当前输入指纹与检查点记录一致**时才允许复用缓存结果；
        指纹不一致、或检查点里根本没有可核对的指纹（旧版本留下的检查点、
        或 checkpoint 侧尚未提供读取接口），一律拒绝复用并全量重跑。
        """
        checkpoint = self._checkpoint_manager
        if checkpoint is None:
            return list(self.config.methods)

        methods_to_run = checkpoint.get_methods_to_run(
            self.config.methods, skip_completed=self._resume_mode
        )

        fingerprint = self._compute_run_fingerprint(tree, alignment_path, calibrations)

        if self._resume_mode:
            skippable = [m for m in self.config.methods if m not in methods_to_run]
            if skippable:
                current = self._run_inputs_are_current(fingerprint)
                if current is not True:
                    reason = (
                        "输入指纹与检查点不一致（本次输入已变化）"
                        if current is False
                        else "检查点中没有可核验的输入指纹"
                    )
                    self.logger.warning(
                        f"Resume refused for {skippable}: {reason}。"
                        f"已缓存结果可能对应另一次输入，忽略它们并全部重新执行。"
                        f"（当前指纹 {fingerprint[:16]}…；如需复用上次结果，"
                        "请在输入未变动时重跑，或使用干净的输出目录。）"
                    )
                    self._record_degradation(
                        f"--resume 未复用任何缓存结果（{reason}）；"
                        f"{len(self.config.methods)} 个方法全部重新执行"
                    )
                    # 记下新指纹：输入再次变回一致时，后续续跑仍可复用。
                    self._store_run_fingerprint(fingerprint, alignment_path)
                    return list(self.config.methods)

        self._store_run_fingerprint(fingerprint, alignment_path)

        if len(methods_to_run) < len(self.config.methods):
            # 只有 --resume 才真的做过指纹核验；非续跑模式下清单变短只可能是
            # 某方法超过 max_retries，此时不能声称"已核验一致"。
            verified_note = " (输入指纹已核验一致)" if self._resume_mode else ""
            self.logger.info(
                f"Checkpoint recovery: running {len(methods_to_run)}/"
                f"{len(self.config.methods)} methods{verified_note}"
            )
            for method_name in self.config.methods:
                if method_name in methods_to_run:
                    continue
                cached_result = checkpoint.get_cached_result(method_name)
                if cached_result:
                    self.results[method_name] = cached_result
                    self.logger.info(
                        f"Restored result for {method_name} from checkpoint"
                    )

        return methods_to_run

    def _ultrametric_check_method(self) -> Optional[str]:
        """为超度量性检查（B-18）挑一个代表方法名。

        ``config.methods`` 可能有多个方法，而检查的分档是"最严格者优先"：
        只要其中有一个通路把"输入树已是时间树"当前提（lsd2/iqtree/mcmctree/viz），
        就按它报 warning；全部是允许非超度量的方法（treepl/pathd8/r8s/mdcat/
        wlogdate）时才降为 info。
        """
        for name in self.config.methods:
            if (name or "").strip().lower() not in _ULTRAMETRIC_TOLERANT_METHODS:
                return name
        return self.config.methods[0] if self.config.methods else None

    def _validate_inputs(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path],
        method: Optional[str] = None,
    ) -> None:
        """验证输入

        Args:
            tree: 待定年的系统发育树
            alignment_path: 比对文件路径
            method: 本次将要使用的定年方法名，用于超度量性检查的分档（B-18）。
                为 None 时按配置的方法列表自行推导；**不能**留空传给校验器，
                否则 treePL/PATHd8/r8s 这类"本就允许非超度量"的方法也会在每次
                运行时被报成 warning。
        """
        validator = TreeValidator()
        report = validator.full_validation(
            tree,
            alignment_path,
            method=method or self._ultrametric_check_method(),
        )

        if not report.is_valid:
            raise PhyloDaterError(f"Input validation failed: {report.errors}")

        for warning in report.warnings:
            self.logger.warning(warning)

        self.logger.success("Input validation passed")

    def _resolve_calibrations(
        self, tree: PhylogeneticTree, calibrations: List[CalibrationPoint]
    ) -> List[CalibrationPoint]:
        """解析校准点（不可变修改）"""
        if not calibrations:
            raise CalibrationError(
                "No valid calibration points provided. "
                "At least one calibration point is required for molecular dating.",
                suggestion=(
                    "请检查校准文件：确保使用支持的约束类型 "
                    "(fixed/uniform/soft_lower/maximum/soft_bounds/gamma/skew_normal/skew_t) "
                    "以及正确的字段名 (age/min/max/alpha/beta/location/scale/shape/df)。"
                ),
            )

        all_pre_resolved = all(cal.resolved_taxa for cal in calibrations)

        # 构建 TaxonomyParser，使用 PipelineConfig 中的参数
        # （三个 kwargs 的取值类型各不相同，只能按 Any 聚在一起传 **）
        taxonomy_kwargs: Dict[str, Any] = {}
        if self.config.taxonomy_delimiter_mode:
            taxonomy_kwargs["delimiter_mode"] = self.config.taxonomy_delimiter_mode
        if self.config.taxonomy_source_priority:
            taxonomy_kwargs["source_priority"] = self.config.taxonomy_source_priority
        if self.config.taxonomy_levels:
            taxonomy_kwargs["custom_levels"] = self.config.taxonomy_levels

        taxonomy_parser = TaxonomyParser(**taxonomy_kwargs)
        if self.config.taxonomy_file:
            try:
                load_kwargs = {}
                if self.config.taxonomy_table_sep:
                    load_kwargs["table_sep"] = self.config.taxonomy_table_sep
                if self.config.taxonomy_sep:
                    load_kwargs["taxonomy_sep"] = self.config.taxonomy_sep
                taxonomy_parser.load_from_file(self.config.taxonomy_file, **load_kwargs)
                self.logger.info(f"Loaded taxonomy from {self.config.taxonomy_file}")
            except (FileNotFoundError, ValueError) as e:
                # 审阅项 B-26 的同族问题在库层：旧实现把加载失败降级成一条
                # warning 后**带着空表继续跑**，于是"分类学表↔类群解析"这项
                # 质量检查在表完全失效时看起来正常通过。分类学表是
                # CalibrationResolver 按类群名定位 MRCA 的输入，载不到表就
                # 等于用户请求的解析依据被整体丢弃 —— 必须失败，或在确实
                # 不依赖该表的分支上留下显式降级记录。
                self.logger.error(f"Failed to load taxonomy file: {e}")
                if not all_pre_resolved:
                    raise CalibrationError(
                        f"分类学表加载失败，无法按类群名解析校准点："
                        f"{self.config.taxonomy_file}（{e}）",
                        suggestion=(
                            "多列表格必须含名称列（name/taxon/organism/label 之一）；"
                            "两列表格第二列须是完整分类串。"
                            "若所有校准点已预先解析（提供 resolved_taxa），"
                            "可去掉 --taxonomy-file 后重跑。"
                        ),
                    ) from e
                # 全部校准点已预解析：本轮确实不需要分类学表，但"用户给了表
                # 却没被用上"这件事不能沉默。
                self._record_degradation(
                    f"分类学表 {self.config.taxonomy_file} 加载失败（{e}），"
                    "本次运行的校准点均已预解析、未使用该表"
                )

        resolver = CalibrationResolver(tree, taxonomy_parser)

        if all_pre_resolved:
            # 预解析的校准点仍需验证祖先冲突；同时为缺少 mrca_leaf_pair 的
            # 校准点补齐（下游适配器依赖 mrca_leaf_pair 在时间树中定位 MRCA，
            # 若缺失会导致解析失败）。已存在则保留原值。
            import copy

            filled = []
            for cal in calibrations:
                if cal.mrca_leaf_pair:
                    filled.append(cal)
                    continue
                new_cal = copy.deepcopy(cal)
                if len(new_cal.resolved_taxa) >= 2:
                    try:
                        new_cal.mrca_leaf_pair = tree.get_mrca_terminals(
                            new_cal.resolved_taxa
                        )
                    except Exception as e:  # 树后端不同，失败原因各异（C-8）
                        # "首/末叶"只是次优的 MRCA 代表叶对，会改变下游定位到的
                        # 节点，因此失败原因与降级都必须可见。
                        self.logger.warning(
                            f"Failed to compute most-distant pair for "
                            f"'{new_cal.name}' ({e}), using first/last taxa as fallback"
                        )
                        self._record_degradation(
                            f"校准点 '{new_cal.name}' 的最远叶对推断失败（{e}），"
                            f"退化为 ({new_cal.resolved_taxa[0]}, "
                            f"{new_cal.resolved_taxa[-1]})，"
                            "MRCA 定位可能与预期不同"
                        )
                        new_cal.mrca_leaf_pair = (
                            new_cal.resolved_taxa[0],
                            new_cal.resolved_taxa[-1],
                        )
                filled.append(new_cal)
            resolver.validate_ancestry(filled)
            self.logger.success(
                f"Validated {len(filled)} pre-resolved calibration points"
            )
            return filled

        import copy

        resolved = []
        for cal in calibrations:
            if not cal.resolved_taxa:
                # 如果是根节点校准点
                if cal.is_root_node:
                    new_cal = copy.deepcopy(cal)
                    new_cal.resolved_taxa = list(tree.tip_names)
                    new_cal.is_root_node = True
                    resolved.append(new_cal)
                # 如果已经有 mrca_leaf_pair，直接使用它来解析
                elif cal.mrca_leaf_pair:
                    # 使用 mrca_leaf_pair 来解析
                    mrca_result = resolver.resolve_mrca_pair(cal.mrca_leaf_pair)
                    new_cal = copy.deepcopy(cal)
                    new_cal.resolved_taxa = mrca_result.resolved_taxa
                    new_cal.mrca_leaf_pair = mrca_result.mrca_leaf_pair
                    new_cal.is_root_node = mrca_result.is_root_node
                    resolved.append(new_cal)
                else:
                    # 尝试通过名称解析
                    result = resolver.resolve(cal.name)
                    new_cal = copy.deepcopy(cal)
                    new_cal.resolved_taxa = result.resolved_taxa
                    new_cal.mrca_leaf_pair = result.mrca_leaf_pair
                    new_cal.is_root_node = result.is_root_node
                    resolved.append(new_cal)
            else:
                # 已有 resolved_taxa：计算正确的 mrca_leaf_pair（距离最远叶节点对）
                new_cal = copy.deepcopy(cal)
                if new_cal.mrca_leaf_pair is None and len(new_cal.resolved_taxa) >= 2:
                    try:
                        new_cal.mrca_leaf_pair = tree.get_mrca_terminals(
                            new_cal.resolved_taxa
                        )
                    except Exception as e:  # 树后端不同，失败原因各异（C-8）
                        # "首/末叶"只是次优的 MRCA 代表叶对，会改变下游定位到的
                        # 节点，因此失败原因与降级都必须可见。
                        self.logger.warning(
                            f"Failed to compute most-distant pair for "
                            f"'{new_cal.name}' ({e}), using first/last taxa as fallback"
                        )
                        self._record_degradation(
                            f"校准点 '{new_cal.name}' 的最远叶对推断失败（{e}），"
                            f"退化为 ({new_cal.resolved_taxa[0]}, "
                            f"{new_cal.resolved_taxa[-1]})，"
                            "MRCA 定位可能与预期不同"
                        )
                        new_cal.mrca_leaf_pair = (
                            new_cal.resolved_taxa[0],
                            new_cal.resolved_taxa[-1],
                        )
                resolved.append(new_cal)

        resolver.validate_ancestry(resolved)

        self.logger.success(f"Resolved {len(resolved)} calibration points")
        return resolved

    def _execute_methods(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path],
        calibrations: List[CalibrationPoint],
    ) -> None:
        """执行所有定年方法，支持断点续跑

        ``--resume`` 的复用判定先经过输入指纹核验（B-23），见
        :meth:`_plan_method_execution`。
        """

        methods_to_run = self._plan_method_execution(tree, alignment_path, calibrations)

        for method_name in methods_to_run:
            self.logger.section(f"Method: {method_name}")

            adapter = None
            try:
                if not DatingMethodRegistry.is_registered(method_name):
                    self.logger.error(f"Method '{method_name}' not registered")
                    continue

                if self._checkpoint_manager:
                    self._checkpoint_manager.start_method(method_name)

                adapter = DatingMethodRegistry.create(
                    method_name,
                    self._tool_config,
                    self.config.output_dir,
                    self._software_paths,
                )

                if not adapter.validate_environment():
                    self.logger.warning(
                        f"{method_name} environment not available, skipping"
                    )
                    if self._checkpoint_manager:
                        self._checkpoint_manager.fail_method(
                            method_name, "Environment not available"
                        )
                    continue

                # 方法级校准兼容性检查（在真正运行前拦截不兼容配置）
                calib_errors = adapter.validate_calibrations(calibrations)
                if calib_errors:
                    error_msg = (
                        f"方法 '{method_name}' 的校准配置不兼容：\n"
                        + "\n".join(f"  - {err}" for err in calib_errors)
                    )
                    self.logger.error(error_msg)
                    if self._checkpoint_manager:
                        self._checkpoint_manager.fail_method(method_name, error_msg)
                    if self._metadata_manager:
                        self._metadata_manager.fail_method(method_name, error_msg)
                    continue

                if self._metadata_manager:
                    self._metadata_manager.start_method(method_name)
                    adapter_version = getattr(adapter, "software_version", None)
                    if adapter_version:
                        self._metadata_manager.set_method_version(
                            method_name, method_name, adapter_version
                        )

                method_start = time.time()
                adapter.prepare_inputs(tree, calibrations, alignment_path)
                adapter.execute()
                result = adapter.parse_results()
                result.execution_seconds = time.time() - method_start

                self.results[method_name] = result

                if self._checkpoint_manager:
                    self._checkpoint_manager.complete_method(method_name, result)

                if self._metadata_manager:
                    self._metadata_manager.complete_method(
                        method_name,
                        execution_seconds=result.execution_seconds,
                        is_converged=result.is_converged,
                        node_count=len(result.node_ages),
                        warnings=(
                            [str(w) for w in result.warnings]
                            if result.warnings
                            else None
                        ),
                    )

                self.logger.success(
                    f"{method_name} completed: {len(result.node_ages)} node ages"
                )

            except Exception as e:
                self.logger.error(f"{method_name} failed: {e}")
                if self._checkpoint_manager:
                    self._checkpoint_manager.fail_method(method_name, str(e))
                if self._metadata_manager:
                    self._metadata_manager.fail_method(method_name, str(e))
                if not self.config.skip_failed:
                    raise

            finally:
                if adapter:
                    try:
                        adapter.cleanup(
                            preserve_intermediates=self.config.save_intermediates
                        )
                    except Exception as cleanup_error:
                        self.logger.debug(
                            f"Cleanup warning for {method_name}: {cleanup_error}"
                        )

    def _save_dated_trees(self) -> None:
        """保存各方法的定年树文件

        根据配置决定输出格式与是否丢弃注释：
        - strip_annotations=False: 输出 .nhx/.nexus 文件（保留 NHX/NEXUS 注释，兼容 iTOL/FigTree）
        - strip_annotations=True: 输出 .nwk 文件（丢弃注释，减小体积）

        另外根据 dated_tree_newick 的实际内容自动判断是 NEXUS 还是 Newick：
        NEXUS 格式（如 lsd2、部分 mcmctree 输出）保存为 .nexus，避免扩展名与格式不符。
        """
        strip_anns = getattr(self.config, "strip_annotations", False)

        for method_name, result in self.results.items():
            if not result.dated_tree_newick:
                continue
            try:
                dated_tree = PhylogeneticTree.from_newick(result.dated_tree_newick)

                # 检测实际内容格式
                tree_text = result.dated_tree_newick.strip()
                is_nexus = (
                    tree_text.upper().startswith("#NEXUS")
                    or "BEGIN TREES" in tree_text.upper()
                )

                if strip_anns:
                    dated_tree = dated_tree.strip_annotations()
                    output_path = (
                        self.config.output_dir / f"{method_name}_dated_tree.nwk"
                    )
                    dated_tree.write(output_path, format="newick")
                elif is_nexus:
                    output_path = (
                        self.config.output_dir / f"{method_name}_dated_tree.nexus"
                    )
                    dated_tree.write(output_path, format="nexus")
                else:
                    output_path = (
                        self.config.output_dir / f"{method_name}_dated_tree.nhx"
                    )
                    dated_tree.write(output_path, format="nhx")

                self.logger.info(f"Saved dated tree: {output_path}")
            except Exception as e:
                self.logger.warning(f"Failed to save dated tree for {method_name}: {e}")

    def _generate_reports(self) -> None:
        """生成报告"""
        if not self.results:
            self.logger.warning("No results to report")
            return

        try:
            from .comparison_reporter import ComparisonReporter

            reporter = ComparisonReporter(self.config.output_dir)
            reporter.generate(self.results, tree=self._tree)
        except Exception as e:
            self.logger.error(f"Report generation failed: {e}")
            self.logger.warning("Continuing despite report generation failure")
            self._record_degradation(
                f"比较报告/汇总产物生成失败（{e}）：本次运行没有标准输出报告"
            )

        # 保存定年树文件（.nhx 或 .nwk，取决于 strip_annotations 配置）
        self._save_dated_trees()

        # 使用运行期间实际使用的配置（而非默认值）保存运行时配置
        try:
            config = Configuration()
            config._final_config = self._tool_config
            config.save_runtime_config(self.config.output_dir / "runtime_config.yaml")
        except Exception as e:
            self.logger.warning(f"Failed to save runtime config: {e}")

        self.logger.success(f"Generated reports for {len(self.results)} methods")


class ParallelPipeline(DatingPipeline):
    """
    并行流水线

    使用多线程并行执行多种定年方法
    每个线程使用独立的随机数种子，确保结果可重现
    """

    MAX_WORKERS = 8

    def _execute_methods(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path],
        calibrations: List[CalibrationPoint],
    ) -> None:
        """并行执行所有定年方法

        使用线程本地存储确保每个线程有独立的随机数状态，
        通过 SeedSequence 生成独立的子种子，保证并行结果的可重现性。
        支持断点续跑：跳过已完成的方法，恢复缓存结果。
        ``--resume`` 的复用判定先经过输入指纹核验（B-23），见
        :meth:`_plan_method_execution`。
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed

        import numpy as np

        # 断点续跑：过滤需要执行的方法（先核验输入指纹，B-23）
        methods_to_run = self._plan_method_execution(tree, alignment_path, calibrations)

        if not methods_to_run:
            self.logger.info("All methods completed or exhausted retries")
            return

        max_workers = min(self.config.threads, self.MAX_WORKERS)
        if self.config.threads > self.MAX_WORKERS:
            self.logger.warning(
                f"threads ({self.config.threads}) exceeds MAX_WORKERS limit ({self.MAX_WORKERS}), "
                f"using {max_workers} workers"
            )

        # 使用 SeedSequence 为每个方法生成独立、可复现的子种子（P0-1 修复见
        # :func:`_derive_method_seeds`）。适配器内部通常直接调用 random /
        # numpy.random 的全局接口，因此下方按「方法名 → 子种子」为每个线程播种
        # 全局随机状态，使并行方法各自拥有独立、可复现的随机流。
        method_seeds = _derive_method_seeds(self.config.seed, methods_to_run)

        metadata_manager = self._metadata_manager
        checkpoint_manager = self._checkpoint_manager

        # 每个并行方法使用 SeedSequence 派生的独立子种子（见上方）。适配器内部
        # 通常直接调用 random / numpy.random 的全局接口，因此在此处按
        # 「方法名 → 子种子」为每个线程播种全局随机状态，使并行方法各自拥有
        # 独立、可复现的随机流，避免此前创建的线程本地 RNG 成为「死 RNG」
        # （计算出来却从未被消费）。

        def run_method(method_name: str) -> tuple:
            """运行单个方法（线程安全，按方法种子初始化全局随机状态）"""
            thread_seed = method_seeds[method_name]

            # 以线程本地的子种子初始化全局随机状态，确保可复现且互不相关。
            # 采用现代 numpy Generator（default_rng）替代 NumPy>=2.0 中已弃用的
            # 全局 np.random.seed() 接口；子种子由主 SeedSequence(self.config.seed)
            # 派生（见上方 seed_seq.spawn），保证各并行 worker 独立且可复现。
            random.seed(thread_seed)
            np.random.seed(thread_seed)

            self.logger.debug(f"{method_name}: using thread-local seed {thread_seed}")

            if checkpoint_manager:
                checkpoint_manager.start_method(method_name)

            adapter = None
            try:
                adapter = DatingMethodRegistry.create(
                    method_name,
                    self._tool_config,
                    self.config.output_dir,
                    self._software_paths,
                )

                if not adapter.validate_environment():
                    self.logger.warning(
                        f"{method_name} environment not available, skipping"
                    )
                    if checkpoint_manager:
                        checkpoint_manager.fail_method(
                            method_name, "Environment not available"
                        )
                    return method_name, None, "Environment not available"

                # 方法级校准兼容性检查（在真正运行前拦截不兼容配置）
                calib_errors = adapter.validate_calibrations(calibrations)
                if calib_errors:
                    error_msg = (
                        f"方法 '{method_name}' 的校准配置不兼容：\n"
                        + "\n".join(f"  - {err}" for err in calib_errors)
                    )
                    self.logger.error(error_msg)
                    if checkpoint_manager:
                        checkpoint_manager.fail_method(method_name, error_msg)
                    if metadata_manager:
                        metadata_manager.fail_method(method_name, error_msg)
                    return method_name, None, error_msg

                if metadata_manager:
                    metadata_manager.start_method(method_name)
                    adapter_version = getattr(adapter, "software_version", None)
                    if adapter_version:
                        metadata_manager.set_method_version(
                            method_name, method_name, adapter_version
                        )

                method_start = time.time()
                adapter.prepare_inputs(tree, calibrations, alignment_path)
                adapter.execute()
                result = adapter.parse_results()
                result.execution_seconds = time.time() - method_start

                if checkpoint_manager:
                    checkpoint_manager.complete_method(method_name, result)

                if metadata_manager:
                    metadata_manager.complete_method(
                        method_name,
                        execution_seconds=result.execution_seconds,
                        is_converged=result.is_converged,
                        node_count=len(result.node_ages),
                        warnings=(
                            [str(w) for w in result.warnings]
                            if result.warnings
                            else None
                        ),
                    )

                return method_name, result, None

            except Exception as e:
                self.logger.error(f"{method_name} failed: {e}")
                if checkpoint_manager:
                    checkpoint_manager.fail_method(method_name, str(e))
                if metadata_manager:
                    metadata_manager.fail_method(method_name, str(e))
                return method_name, None, str(e)

            finally:
                if adapter:
                    try:
                        adapter.cleanup(
                            preserve_intermediates=self.config.save_intermediates
                        )
                    except Exception as cleanup_error:
                        # 与串行分支保持一致：清理失败不得掩盖方法结果，
                        # 但必须留痕（C-8：此前的 ``except Exception: pass``）。
                        self.logger.debug(
                            f"Cleanup warning for {method_name}: {cleanup_error}"
                        )

        # 并行执行
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(run_method, name): name for name in methods_to_run
            }

            failed_methods = []
            for future in as_completed(futures):
                method_name, result, error_msg = future.result()
                if result:
                    self.results[method_name] = result
                    self.logger.success(f"{method_name} completed")
                else:
                    failed_methods.append((method_name, error_msg))

            if failed_methods:
                error_summary = (
                    f"{len(failed_methods)} methods failed: "
                    f"{', '.join(f'{n} ({e})' if e else n for n, e in failed_methods)}"
                )
                if not self.config.skip_failed:
                    raise PhyloDaterError(error_summary)
                self.logger.warning(error_summary)
