"""
DatingMethod 接口和 Registry

定义所有定年适配器必须实现的接口
"""

from __future__ import annotations

import numbers
import os
import uuid
from abc import ABC, abstractmethod
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Generic,
    List,
    Optional,
    Type,
    TypeVar,
    Union,
)

if TYPE_CHECKING:
    import dendropy

from ..infrastructure import get_logger
from ..infrastructure.configuration import CommonConfig, SoftwarePaths, ToolConfig
from ..models import CalibrationPoint, DatingResult, NodeAgeEstimate, PhylogeneticTree
from .exceptions import ExecutionError

#: 适配器**实际读取**的那一份配置类型（``ToolConfig`` 的子配置，例如
#: ``TreePLConfig``）。基类用它参数化 ``self.config``，使每个适配器都能声明
#: "我看到的 config 就是这一种"，而不必把 ``self.config`` 降级成 ``Any``。
#: 设计取舍见 ``docs/ARCHITECTURE.md`` §Key Design Decisions 第 5 节
#: （中文版：``docs/ARCHITECTURE_CN.md`` §关键设计决策 5「适配器配置类型」）。
ConfigT = TypeVar("ConfigT")

#: ``prepare_inputs()`` 返回的“输入文件信息字典”里允许出现的值。它既装文件
#: 路径，也装适配器随调用带出的**执行统计量**（位点数、耗时、对数似然）——
#: r8s 与 MD-Cat 就是两者共用一份字典。因此基类 ``cleanup()`` 归档时必须认出
#: 并跳过非路径条目：旧实现无条件 ``Path(value)``，而 ``Path(1000)`` 抛的是
#: ``TypeError``，并不在 cleanup 自己的 (OSError, shutil.Error, ValueError)
#: 里，一开 preserve_intermediates 就会直接砸穿归档。
InputFileInfo = Union[Path, int, float]


class DatingMethod(ABC, Generic[ConfigT]):
    """
    定年方法抽象基类

    所有定年适配器必须实现此接口

    类型参数 ``ConfigT`` 是本适配器读取的配置形态。适配器入口保留双形态
    （``ToolConfig`` 或单个方法的子配置），但两者在 ``__init__`` 里就被归一成
    ``ConfigT``，因此后续代码里的 ``self.config.<参数>`` 在静态视角下始终成立。
    """

    def __init__(
        self,
        config: ConfigT,
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
    ) -> None:
        self.config = config
        # 通用配置在此归一为非 Optional：``common_config`` 缺失时用一份默认值，
        # 而不是把 ``None`` 沿用到读取 nthreads/timeout/seed 的各处——那些读取
        # 没有"跳过"语义，线程数与超时必须有值。
        self.common_config: CommonConfig = (
            common_config if common_config is not None else CommonConfig()
        )
        # 解析为绝对路径，避免 work_dir 为相对路径时与 ProcessRunner(cwd=work_dir)
        #   + 含 work_dir 前缀的相对路径参数组合导致路径重复拼接（FileNotFoundError）。
        self.output_dir = Path(output_dir).resolve()
        self.work_dir = self._create_unique_work_dir(self.output_dir / self.method_name)
        self._software_paths = software_paths or SoftwarePaths()
        # 后端软件版本横幅/版本号，由 ``validate_environment()`` 在探测成功时写入，
        # 并被 pipeline 记录到 runtime_metadata.json 的 ``software_versions``。
        # 定年结果只有在知道是哪个版本的引擎产出时才可复现，因此这里不留空。
        self.software_version: Optional[str] = None

    @staticmethod
    def _package_version(package: str) -> Optional[str]:
        """返回 Python 发行包（pyr8s / wLogDate / MD-Cat）的元数据版本号。"""
        try:
            from importlib import metadata

            return metadata.version(package)
        except Exception:  # pragma: no cover - 包未通过 pip/conda 安装
            return None

    def record_software_version(
        self,
        executable: str,
        version_flag: str = "--version",
        package: Optional[str] = None,
    ) -> Optional[str]:
        """探测并缓存后端软件版本，返回记录的版本字符串。

        优先执行 ``<executable> <version_flag>`` 取第一行非空输出；二进制不可执行
        而软件是以 Python 包形式提供时（pyr8s / wLogDate / MD-Cat），回退到该发行
        包的元数据版本。两条路都失败时返回 ``None``（调用方保持原状，不写假值）。
        """
        version: Optional[str] = None
        try:
            from ..infrastructure import ProcessRunner

            runner = ProcessRunner(timeout=30)
            output = runner.get_version(executable, version_flag)
            if output:
                version = str(output).strip().splitlines()[0].strip()
        except Exception as exc:  # pragma: no cover - 版本探测不应阻断运行
            get_logger().debug(f"Version probe for '{executable}' failed: {exc}")

        if not version and package:
            version = self._package_version(package)

        if version:
            self.software_version = version
        return version

    def _create_unique_work_dir(self, base_dir: Path) -> Path:
        """
        创建具有唯一标识的工作目录

        使用时间戳和随机后缀确保并行任务不会互相覆盖
        格式: {base_dir}_{timestamp}_{random_suffix}

        工作目录会被注册到全局清理列表，确保异常退出时自动清理。
        """
        import time

        from .pipeline import _register_temp_directory

        timestamp = time.strftime("%Y%m%d_%H%M%S")
        unique_suffix = uuid.uuid4().hex[:8]
        unique_dir = Path(str(base_dir) + f"_{timestamp}_{unique_suffix}")
        unique_dir.mkdir(parents=True, exist_ok=True)

        # 注册到全局清理列表，确保异常退出时清理
        _register_temp_directory(unique_dir)

        # 显式设置临时目录权限为 0o700，确保仅当前用户可访问
        try:
            os.chmod(str(unique_dir), 0o700)
        except OSError as exc:
            # C-8：只兜住"权限位设不上"这一类（NFS/特殊文件系统上确实会发生），
            # 并留下痕迹；不再用裸 except Exception 吞掉编码/路径类真实缺陷。
            get_logger().debug(
                f"Could not chmod 0o700 the working directory {unique_dir}: {exc}"
            )

        return unique_dir

    def get_software_path(self, tool_name: str) -> Optional[str]:
        """
        获取软件路径，支持自定义路径配置

        Args:
            tool_name: 软件路径属性名（如 'mcmctree_bin', 'iqtree_bin'）

        Returns:
            自定义路径或 None
        """
        return getattr(self._software_paths, tool_name, None)

    # ------------------------------------------------------------------
    # ``_input_files`` 的类型安全读取
    #
    # 这份字典同时装文件路径与执行统计量（见 InputFileInfo），因此每个
    # 读取处都必须显式收窄。三个 helper 把“收窄不成就是内部状态被写坏”
    # 统一报成 ExecutionError：比拿它去 open() 后收到一个语义不明的
    # FileNotFoundError、或拿到 ``None`` 后一路走到 TypeError 更可诊断。
    # ------------------------------------------------------------------

    def _input_files_map(self) -> Dict[str, InputFileInfo]:
        """返回子类的 ``_input_files``（未准备过时为空字典）。"""
        return getattr(self, "_input_files", {})

    def _require_input_path(self, key: str) -> Path:
        """取 ``_input_files`` 中的文件条目并保证它是个路径。"""
        entry = self._input_files_map()[key]
        if not isinstance(entry, (str, os.PathLike)):
            raise ExecutionError(
                f'{self.method_name}: internal state _input_files["{key}"] must be '
                f"a path, got {type(entry).__name__}: {entry!r}"
            )
        return Path(entry)

    def _require_input_number(self, key: str, default: float = 0.0) -> float:
        """取 ``_input_files`` 中的数值型统计量（耗时、对数似然等）。

        用 :class:`numbers.Real` 而不是 ``float`` 做守卫：上游可能递回 numpy
        标量，它不是 ``float`` 子类但完全可转。
        """
        entry = self._input_files_map().get(key, default)
        if isinstance(entry, bool) or not isinstance(entry, numbers.Real):
            raise ExecutionError(
                f'{self.method_name}: internal state _input_files["{key}"] must be '
                f"a number, got {type(entry).__name__}: {entry!r}"
            )
        return float(entry)

    def _require_input_int(self, key: str) -> int:
        """取 ``_input_files`` 中的整型条目（如位点数）。

        必须保持 ``int``：它常被写进引擎的 ``-s`` 一类参数，变成 float 就与
        上游习惯不一致（"1000" 变 "1000.0"）。
        """
        entry = self._input_files_map()[key]
        if isinstance(entry, bool) or not isinstance(entry, numbers.Integral):
            raise ExecutionError(
                f'{self.method_name}: internal state _input_files["{key}"] must be '
                f"an int, got {type(entry).__name__}: {entry!r}"
            )
        return int(entry)

    def _extract_chronogram_ages(
        self,
        tree: "dendropy.Tree",
        calibrations: List["CalibrationPoint"],
        unit_factor: float = 1.0,
        key_by_cal_name: bool = True,
    ) -> Dict["str", "NodeAgeEstimate"]:
        """从 chronogram（年龄编码在分支长度中）按 MRCA 拓扑解码节点年龄。

        部分定年软件（r8s/lsd2 等）输出的时间树是 chronogram：根节点年龄 = 0，
        向叶节点方向时间递增，年龄编码在分支长度的累加（root-to-node 距离）中，
        节点标签内并不写入 ``age`` 属性或 ``[&age=...]`` 注解。

        本方法通过每个校准点的 ``mrca_leaf_pair`` 在时间树中定位 MRCA 节点，
        取其 ``distance_from_root()`` 作为节点年龄。用作各适配器 parse_results
        的回退（fallback）：当标签级解析未提取到任何年龄时，避免结果被丢弃。

        Args:
            tree: 定年后的时间树（dendropy.Tree，taxon_namespace 含叶节点）。
            calibrations: 校准点列表。
            unit_factor: 距离到 Ma 的换算系数（r8s/lsd2 为 1，mcmctree 的 FigTree
                输出为 Ga，故传入 1000 转 Ma）。
            key_by_cal_name: 为 True 时以校准点名作为 node_ages 的键；回退到节点
                标签（label）作为键仅在前者缺失时。

        Returns:
            {节点键: NodeAgeEstimate}。仅包含成功定位到 MRCA 的校准点。
        """
        from ..models.results import CIType, NodeAgeEstimate

        node_ages: Dict[str, "NodeAgeEstimate"] = {}
        if tree is None or not calibrations:
            return node_ages

        for cal in calibrations:
            # 根节点校准：直接取 seed_node
            if cal.is_root_node and cal.name:
                # 根节点到自身的距离为 0，因此 distance_from_root() 恒为 0；
                # 应使用根到最远叶节点的距离（即树的根年龄）作为校准年龄。
                age = float(tree.max_distance_from_root()) * unit_factor
                key = (
                    cal.name
                    if key_by_cal_name
                    else (tree.seed_node.label or f"node_{id(tree.seed_node)}")
                )
                if key and key not in node_ages:
                    node_ages[key] = NodeAgeEstimate(mean_age=age, ci_type=CIType.NONE)
                continue

            pair = cal.mrca_leaf_pair
            if not pair or len(pair) < 2:
                continue

            # 在时间树的 taxon_namespace 中按名称定位叶节点
            taxa = []
            for tip_name in pair:
                taxon = tree.taxon_namespace.get_taxon(tip_name)
                if taxon is None:
                    get_logger().debug(
                        f"Chronogram fallback: tip '{tip_name}' of calibration "
                        f"'{cal.name}' is not in the dated tree's taxon namespace; "
                        "this calibration will not appear under its name"
                    )
                    break
                taxa.append(taxon)
            if len(taxa) != len(pair):
                continue

            try:
                mrca = tree.mrca(taxa=taxa)
            except (ValueError, KeyError, AttributeError) as exc:
                # C-8：dendropy 的 mrca() 在 taxa 不在同一棵树/命名空间时抛这类异常。
                # 该节点因此不会以校准点名入表——必须留下痕迹，否则"报告里看不到
                # 校准名"就无从追查（与 mcmctree 的同类分支保持一致）。
                get_logger().debug(
                    f"Chronogram fallback: cannot locate MRCA for calibration "
                    f"'{cal.name}' from tips {pair}: {exc}"
                )
                continue
            if mrca is None:
                get_logger().debug(
                    f"Chronogram fallback: MRCA of tips {pair} (calibration "
                    f"'{cal.name}') is undefined in the dated tree; node skipped"
                )
                continue

            try:
                age = float(mrca.distance_from_root()) * unit_factor
            except (TypeError, ValueError) as exc:
                # 分支长度缺失/非数值：只影响这一个节点，但要说明是哪个
                get_logger().debug(
                    f"Chronogram fallback: no usable root distance for calibration "
                    f"'{cal.name}' node: {exc}"
                )
                continue
            key = (
                cal.name
                if (key_by_cal_name and cal.name)
                else (mrca.label or f"node_{id(mrca)}")
            )
            if key and key not in node_ages:
                node_ages[key] = NodeAgeEstimate(mean_age=age, ci_type=CIType.NONE)

        return node_ages

    @property
    @abstractmethod
    def method_name(self) -> str:
        """方法标识符"""
        pass

    @abstractmethod
    def validate_environment(self) -> bool:
        """
        检测外部程序是否可用并验证版本

        Returns:
            True 如果环境有效，False 如果不可用
        """
        pass

    @abstractmethod
    def prepare_inputs(
        self,
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        alignment_path: Optional[Path] = None,
    ) -> Dict:
        """
        生成特定软件所需的输入文件和配置文件

        Args:
            tree: 系统发育树
            calibrations: 校准点列表
            alignment_path: 比对文件路径（可选）

        Returns:
            输入文件信息字典
        """
        pass

    @abstractmethod
    def execute(self) -> bool:
        """
        执行定年分析

        Returns:
            True 如果执行成功
        """
        pass

    @abstractmethod
    def parse_results(self) -> DatingResult:
        """
        解析输出为标准 DatingResult 对象

        Returns:
            定年结果
        """
        pass

    def validate_calibrations(self, calibrations: List[CalibrationPoint]) -> List[str]:
        """校验校准点是否满足本方法的特殊需求。

        各适配器可重写此方法，在真正运行外部软件前检查校准配置是否兼容。
        返回错误信息列表；空列表表示通过。发现错误时应给出可操作的修改建议。

        Args:
            calibrations: 已解析的校准点列表

        Returns:
            错误信息字符串列表
        """
        return []

    def cleanup(self, preserve_intermediates: bool = False) -> Dict[str, Path]:
        """
        归档中间文件，清理临时目录

        Args:
            preserve_intermediates: 是否保留中间文件

        Returns:
            保留的中间文件路径字典 {description: path}
        """
        if not preserve_intermediates:
            return {}

        # 通用归档：将 work_dir 中的中间文件复制到 output_dir/intermediates/<method>/
        import shutil

        intermediates_dir = self.output_dir / "intermediates" / self.method_name
        intermediates_dir.mkdir(parents=True, exist_ok=True)

        archived: Dict[str, Path] = {}
        logger = get_logger()
        input_files: Dict[str, InputFileInfo] = getattr(self, "_input_files", {})
        for desc, src_path in input_files.items():
            if not isinstance(src_path, (str, os.PathLike)):
                # 执行统计量不是文件，不参与归档；work_dir 里的真实产物由下面
                # 的“逐文件扫描”那一段兼顾，不会因此丢失。
                continue
            try:
                src = Path(src_path)
                if src.exists() and src.is_file():
                    dst = intermediates_dir / src.name
                    shutil.copy2(src, dst)
                    archived[desc] = dst
            except (OSError, shutil.Error, ValueError) as exc:
                # C-8：归档失败只影响"中间文件留不留得住"，但静默吞掉会让用户
                # 以为归档成功。逐条记 warning，不再 pass。
                logger.warning(
                    f"Failed to archive intermediate '{desc}' ({src_path}): {exc}"
                )

        # 同时归档 work_dir 下的其他文件（.ctl, .nex, .log 等）
        try:
            work_dir = Path(self.work_dir)
            if work_dir.exists():
                for f in work_dir.iterdir():
                    if f.is_file() and f.name not in [
                        p.name for p in archived.values()
                    ]:
                        try:
                            dst = intermediates_dir / f.name
                            shutil.copy2(f, dst)
                            archived[f.name] = dst
                        except (OSError, shutil.Error) as exc:
                            logger.warning(
                                f"Failed to archive work-dir file '{f}': {exc}"
                            )
        except OSError as exc:
            logger.warning(
                f"Failed to enumerate work directory {getattr(self, 'work_dir', None)} "
                f"for archiving: {exc}"
            )

        return archived


class DatingMethodRegistry:
    """
    定年方法注册表

    工厂注册机制，遵守开放封闭原则
    """

    _registry: Dict[str, Type[DatingMethod]] = {}

    @classmethod
    def register(cls, name: str, method_class: Type[DatingMethod]) -> None:
        """
        注册定年方法

        Args:
            name: 方法名称
            method_class: 方法类
        """
        cls._registry[name.lower()] = method_class

    @classmethod
    def create(
        cls,
        name: str,
        config: ToolConfig,
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[Any] = None,
    ) -> DatingMethod:
        """
        创建定年方法实例

        Args:
            name: 方法名称
            config: 工具配置
            output_dir: 输出目录
            software_paths: 外部软件路径配置
            common_config: 通用配置（线程数、随机种子、verbose等）

        Returns:
            DatingMethod 实例

        Raises:
            UnknownMethodError: 如果方法未注册
        """
        from .exceptions import UnknownMethodError

        name_lower = name.lower()
        if name_lower not in cls._registry:
            raise UnknownMethodError(f"Unknown dating method: {name}")

        method_class = cls._registry[name_lower]

        # 获取通用配置：优先使用传入的 common_config，否则从 ToolConfig 中获取
        if common_config is None:
            common_config = getattr(config, "common", None)

        # 从 ToolConfig 中提取对应方法的子配置
        method_config_map = {
            "mcmctree": "mcmctree",
            "pathd8": "pathd8",
            "lsd2": "lsd2",
            "wlogdate": "wlogdate",
            "mdcat": "mdcat",
            "r8s": "r8s",
            "pyr8s": "r8s",
            "treepl": "treepl",
        }

        config_attr = method_config_map.get(name_lower, name_lower)
        method_config = getattr(config, config_attr, config)

        # 检查适配器是否接受 common_config 参数
        import inspect

        sig = inspect.signature(method_class.__init__)
        if "common_config" in sig.parameters:
            return method_class(
                method_config, output_dir, software_paths, common_config=common_config
            )
        else:
            return method_class(method_config, output_dir, software_paths)

    @classmethod
    def list_methods(cls) -> List[str]:
        """列出所有注册的方法"""
        return list(cls._registry.keys())

    @classmethod
    def is_registered(cls, name: str) -> bool:
        """检查方法是否已注册"""
        return name.lower() in cls._registry
