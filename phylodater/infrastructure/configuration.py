"""
Configuration - 三层配置合并机制

支持：内置默认值 < YAML配置文件 < 命令行参数
"""

import copy
import math
import typing
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Type, TypeVar, Union

import yaml

from ..core.exceptions import UnknownParameterWarning
from .logging import get_logger
from .safe_io import safe_writer

T = TypeVar("T")


def _as_float(value: Any, name: str) -> float:
    """把配置值转成**有限**浮点数（审阅项 B-19 同族）。

    ``--method-args`` 与 YAML 都能把 ``"nan"`` / ``".nan"`` 变成
    ``float('nan')``，而所有 ``<=`` / ``>=`` 阈值比较对 NaN 恒为 False，
    于是"PSRF <= 阈值"会一直判未收敛（或一直判收敛，取决于比较方向），
    用户却看不到任何异常。这里在配置层就把它拦下。
    """
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a number, got bool {value!r}")
    try:
        numeric = float(value)
    except (TypeError, ValueError) as e:
        raise ValueError(f"{name} must be a number, got {value!r}") from e
    if not math.isfinite(numeric):
        raise ValueError(
            f"{name} must be a finite number, got {value!r}: NaN and infinity "
            "silently break every threshold comparison in the convergence checks"
        )
    return numeric


@dataclass
class CommonConfig:
    """所有软件共享的通用配置"""

    nthreads: int = 4  # 统一线程数参数名
    seed: Optional[int] = None  # 随机种子
    verbose: bool = False  # 详细输出
    timeout: int = 7200  # 执行超时时间（秒）

    def __post_init__(self) -> None:
        """验证参数范围"""
        if self.nthreads < 1:
            raise ValueError(f"nthreads must be >= 1, got {self.nthreads}")
        if self.seed is not None and self.seed < 0:
            raise ValueError(f"seed must be non-negative, got {self.seed}")
        if self.timeout < 1:
            raise ValueError(f"timeout must be >= 1, got {self.timeout}")


@dataclass
class MCMCTreeConfig:
    """MCMCTree 配置"""

    clock: int = 2  # 1=严格钟, 2=独立松弛钟, 3=自相关松弛钟
    num_runs: int = 2
    burnin: int = 20000
    nsample: int = 50000
    sampfreq: int = 5
    ndata: int = 1
    root_age: Optional[float] = None
    seqtype: str = "auto"  # auto, dna, protein (注意: MCMCtree 不支持 codon)
    model: str = "auto"
    clean_data: int = 0  # 0=保留所有位点, 1=删除含gap位点
    rate_alpha: float = 2.0
    rate_scale: float = 1.0
    # 阶段 0 速率估算开关：默认 False=运行 baseml/codeml (clock=1) 估算总体替换
    # 速率并用其均值设置 rgene_gamma；True=跳过估算，直接用默认先验 rgene_gamma。
    skip_rate_estimation: bool = False
    time_unit: str = "Ma"  # Ma 或 Ga，影响 rgene_gamma 计算
    paml_version: str = "4.10.8"  # PAML 版本，影响控制文件格式

    # PAML 4.10.8+ 新增参数
    bd_construction: str = "c"  # "c"=conditional, "m"=multiplicative (4.10.8+)
    checkpoint_enabled: bool = False
    checkpoint_mode: int = 1  # 1=save, 2=resume
    checkpoint_prob: float = 0.01
    checkpoint_file: str = "mcmctree.ckpt1"
    duplication: int = 0  # 0=no, 1=yes (for paralogous genes)

    # ---- 收敛判定阈值（审阅项 B-4；适配器经 getattr 读取同名属性）----------
    # ``adapters/mcmctree_method.py`` 以
    # ``_positive_config(self.config, "psrf_threshold", DEFAULT_PSRF_THRESHOLD)``
    # / ``"ess_threshold"`` 读取这两项，此前配置层根本没有这两个字段，
    # 于是永远退回适配器默认值、用户既无法在 YAML 里写、也无法用
    # ``--method-args`` 传（白名单会拒绝）。默认值与适配器常量保持一致：
    # PSRF<=1.05、ESS>=100（MCMC 定年文献的宽松惯例；严格实践取 1.01）。
    psrf_threshold: float = 1.05  # Gelman–Rubin 潜在尺度缩减因子上限
    ess_threshold: float = 100.0  # 最小有效样本量（ESS）下限

    def __post_init__(self) -> None:
        """验证参数范围"""
        # 字符串形式（如 CLI --method-args 直接构造）规范化为 bool，
        # 语义与 Configuration._coerce_value 保持一致
        if isinstance(self.skip_rate_estimation, str):
            self.skip_rate_estimation = self.skip_rate_estimation.strip().lower() in (
                "true",
                "1",
                "yes",
                "on",
            )
        # 阈值同理：--method-args 透传进来的永远是字符串
        self.psrf_threshold = _as_float(self.psrf_threshold, "psrf_threshold")
        self.ess_threshold = _as_float(self.ess_threshold, "ess_threshold")
        if self.clock not in [1, 2, 3]:
            raise ValueError(f"clock must be 1, 2, or 3, got {self.clock}")
        if self.num_runs < 1:
            raise ValueError(f"num_runs must be >= 1, got {self.num_runs}")
        if self.burnin < 0:
            raise ValueError(f"burnin must be >= 0, got {self.burnin}")
        if self.nsample < 1:
            raise ValueError(f"nsample must be >= 1, got {self.nsample}")
        if self.sampfreq < 1:
            raise ValueError(f"sampfreq must be >= 1, got {self.sampfreq}")
        if self.ndata < 1:
            raise ValueError(f"ndata must be >= 1, got {self.ndata}")
        if self.clean_data not in [0, 1]:
            raise ValueError(f"clean_data must be 0 or 1, got {self.clean_data}")
        if self.rate_alpha <= 0:
            raise ValueError(f"rate_alpha must be > 0, got {self.rate_alpha}")
        if self.rate_scale <= 0:
            raise ValueError(f"rate_scale must be > 0, got {self.rate_scale}")
        if self.time_unit not in ["Ma", "Ga"]:
            raise ValueError(f"time_unit must be 'Ma' or 'Ga', got {self.time_unit}")
        if self.checkpoint_mode not in [0, 1, 2]:
            raise ValueError(
                f"checkpoint_mode must be 0, 1, or 2, got {self.checkpoint_mode}"
            )
        if not 0 <= self.checkpoint_prob <= 1:
            raise ValueError(
                f"checkpoint_prob must be between 0 and 1, got {self.checkpoint_prob}"
            )
        # PSRF 恒 >= 1，阈值取 <=1 意味着"任何链都判为未收敛"，
        # 属于会静默改变结论的无意义配置，因此当场拒绝而不是放行。
        if self.psrf_threshold <= 1.0:
            raise ValueError(
                f"psrf_threshold must be > 1.0 (PSRF is >= 1 by construction, so a "
                f"threshold of {self.psrf_threshold} can never be met), got "
                f"{self.psrf_threshold}"
            )
        if self.ess_threshold < 1:
            raise ValueError(f"ess_threshold must be >= 1, got {self.ess_threshold}")


@dataclass
class R8sConfig:
    """r8s/pyr8s 配置"""

    # 默认 NPRS：pyr8s 后端仅支持 NPRS 算法，若默认 PL 会被静默替换。
    # 三个取值对应 r8s 1.7 手册 `divtime method=LF|NPRS|PL algorithm=...`
    # （审阅项 A-6 修复建议 3：原版 r8s 的 LF = penalized-likelihood 变体
    # 此前根本不在允许值里，用户按手册写 method=LF 会被配置层直接拒绝）。
    method: str = "NPRS"  # LF / NPRS / PL（LF 与 PL 仅原生 r8s 子进程通路有效）
    algorithm: str = "TN"
    r8s_bin: str = "r8s"
    smoothing: Optional[float] = None

    def __post_init__(self) -> None:
        """验证参数范围"""
        # 大小写不敏感：手册与既往 YAML 里出现过 LF/NPRS/pl 等多种写法
        if isinstance(self.method, str):
            self.method = self.method.strip().upper()
        if self.method not in ["LF", "PL", "NPRS"]:
            raise ValueError(
                f"method must be 'LF', 'PL' or 'NPRS', got '{self.method}'"
            )
        if self.smoothing is not None and self.smoothing <= 0:
            raise ValueError(f"smoothing must be > 0, got {self.smoothing}")


@dataclass
class TreePLConfig:
    """treePL 配置

    参考: treePL原始代码 main.cpp
    """

    treepl_bin: str = "treepl"
    initial_smooth: float = 100.0

    # 树文件相关参数
    numsites: int = 0  # 位点数，用于处理零分支长度 (默认0=不处理)
    scale: float = 1.0  # 缩放因子
    collapse: bool = False  # 折叠零分支长度
    set1: bool = False  # 将分支长度设为1

    # CV分析参数
    cv: bool = True  # 是否进行交叉验证
    randomcv: bool = True  # 使用随机子采样CV (推荐)
    cvstart: float = 1000.0  # CV起始平滑值
    cvstop: float = 0.1  # CV终止平滑值
    cvmultstep: float = 0.1  # CV步长乘数
    cviter: int = 10  # CV迭代次数 (randomcviter)
    cvsimaniter: int = 0  # CV模拟退火迭代次数
    randomcvsamp: float = 0.1  # 随机CV采样比例
    cvoutfile: str = "cv.out"  # CV输出文件
    cv_replicates: int = 3  # CV重复运行次数，用于检查平滑值稳定性

    # Prime优化参数
    run_prime: bool = True  # 是否运行prime步骤
    prime_replicates: int = 3  # Prime重复运行次数，选择最低opt/optad
    opt: Optional[int] = None  # 优化参数opt，None表示由prime自动确定
    optad: Optional[int] = None  # 优化参数optad，None表示由prime自动确定
    optcvad: Optional[int] = None  # CV优化参数optcvad

    # 自动调整参数
    auto_adjust_opt: bool = True  # 自动调整opt/optad以消除警告
    max_opt_adjust: int = 10  # opt/optad最大调整次数

    # 优化器详细参数
    moredetail: bool = False  # LF更详细的优化
    moredetailad: bool = False  # PL更详细的优化
    moredetailcvad: bool = False  # CV更详细的优化
    calcgrad: bool = False  # 使用梯度而非自动微分

    # 模拟退火参数
    lftemp: Optional[float] = None  # LF起始温度
    pltemp: Optional[float] = None  # PL起始温度
    lfcool: Optional[float] = None  # LF冷却速率
    plcool: Optional[float] = None  # PL冷却速率
    lfstoptemp: Optional[float] = None  # LF停止温度
    plstoptemp: Optional[float] = None  # PL停止温度
    lfrtstep: Optional[float] = None  # LF速率步长
    lfdtstep: Optional[float] = None  # LF日期步长
    plrtstep: Optional[float] = None  # PL速率步长
    pldtstep: Optional[float] = None  # PL日期步长

    # 迭代次数参数
    lfiter: Optional[int] = None  # LF完整迭代次数
    pliter: Optional[int] = None  # PL完整迭代次数
    lfsimaniter: Optional[int] = None  # LF模拟退火迭代次数
    plsimaniter: Optional[int] = None  # PL模拟退火迭代次数

    # 收敛容差
    ftol: Optional[float] = None  # 函数容差
    xtol: Optional[float] = None  # 参数容差

    # 其他参数
    thorough: bool = False  # 详细优化 (MAY TAKE A WHILE)
    paramverbose: bool = False  # 参数详细输出
    checkconstraints: bool = False  # 检查约束有效性
    mapspace: bool = False  # 生成映射空间
    log_pen: bool = False  # 使用对数惩罚
    sample: float = 1.0  # 采样比例

    # 输入文件参数
    ind8s: Optional[str] = None  # r8s日期输入文件
    inr8s: Optional[str] = None  # r8s速率输入文件
    outfile: str = "out_dates.tre"  # 输出文件名

    # Tiny branch处理 (phylodater扩展)
    check_tiny_branches: bool = True  # 检查并处理tiny branch
    min_branch_length_threshold: float = 1e-5  # 最小分支长度阈值
    branch_length_multiplier: float = 1000.0  # 分支长度乘数

    # CV自适应参数 (phylodater扩展)
    adaptive_cv: bool = True  # 自适应调整CV参数
    cvstop_min: float = 1e-36  # cvstop最小值

    def __post_init__(self) -> None:
        """验证参数范围"""
        if self.initial_smooth <= 0:
            raise ValueError(f"initial_smooth must be > 0, got {self.initial_smooth}")
        if self.numsites < 0:
            raise ValueError(f"numsites must be >= 0, got {self.numsites}")
        if self.scale <= 0:
            raise ValueError(f"scale must be > 0, got {self.scale}")
        if self.cvstart <= 0:
            raise ValueError(f"cvstart must be > 0, got {self.cvstart}")
        if self.cvstop <= 0:
            raise ValueError(f"cvstop must be > 0, got {self.cvstop}")
        if self.cvstop >= self.cvstart:
            raise ValueError(
                f"cvstop ({self.cvstop}) must be < cvstart ({self.cvstart})"
            )
        if self.cvmultstep <= 0 or self.cvmultstep >= 1:
            raise ValueError(
                f"cvmultstep must be between 0 and 1, got {self.cvmultstep}"
            )
        if self.cviter < 1:
            raise ValueError(f"cviter must be >= 1, got {self.cviter}")
        if self.cv_replicates < 1:
            raise ValueError(f"cv_replicates must be >= 1, got {self.cv_replicates}")
        if self.prime_replicates < 1:
            raise ValueError(
                f"prime_replicates must be >= 1, got {self.prime_replicates}"
            )
        if self.max_opt_adjust < 0:
            raise ValueError(f"max_opt_adjust must be >= 0, got {self.max_opt_adjust}")
        if self.sample <= 0 or self.sample > 1:
            raise ValueError(f"sample must be between 0 and 1, got {self.sample}")
        if self.min_branch_length_threshold <= 0:
            raise ValueError(
                f"min_branch_length_threshold must be > 0, got {self.min_branch_length_threshold}"
            )
        if self.branch_length_multiplier <= 0:
            raise ValueError(
                f"branch_length_multiplier must be > 0, got {self.branch_length_multiplier}"
            )


@dataclass
class PATHd8Config:
    """PATHd8 配置"""

    pathd8_bin: str = "PATHd8"

    def __post_init__(self) -> None:
        """验证参数范围"""
        if not self.pathd8_bin or not self.pathd8_bin.strip():
            raise ValueError("pathd8_bin cannot be empty")


@dataclass
class LSD2Config:
    """LSD2 配置"""

    # 替换模型。None = 按比对序列类型自动选择（核酸 -> GTR+G，蛋白 -> LG+G）。
    # IQ-TREE2 的 --date 流程必须显式给出 -m（否则跑 modelFinder），但模型与
    # 数据类型不匹配时 IQ-TREE2 会把模型名当成文件路径并报
    # "ERROR: File not found LG"，因此默认不写死氨基酸模型。
    model: Optional[str] = None
    iqtree_bin: str = "iqtree2"
    date_ci: int = 100  # 置信区间重采样次数
    clock_sd: float = 0.2
    date_options: Optional[str] = None
    date_outlier: Optional[str] = None

    def __post_init__(self) -> None:
        """验证参数范围"""
        if self.date_ci < 0:
            raise ValueError(f"date_ci must be >= 0, got {self.date_ci}")
        if self.clock_sd <= 0:
            raise ValueError(f"clock_sd must be > 0, got {self.clock_sd}")


@dataclass
class WLogDateConfig:
    """wLogDate 配置"""

    wlogdate_bin: str = "launch_wLogDate.py"
    backward_time: bool = True
    num_replicates: int = 1
    root_time: Optional[float] = None
    leaf_time: float = 0.0
    sequence_length: Optional[int] = None
    max_iter: int = 50000
    pseudocount: float = 0.01
    zero_branch_len: float = 1e-10

    def __post_init__(self) -> None:
        """验证参数范围"""
        if self.num_replicates < 1:
            raise ValueError(f"num_replicates must be >= 1, got {self.num_replicates}")
        if self.sequence_length is not None and self.sequence_length <= 0:
            raise ValueError(f"sequence_length must be > 0, got {self.sequence_length}")
        if self.max_iter < 1:
            raise ValueError(f"max_iter must be >= 1, got {self.max_iter}")
        if self.pseudocount < 0:
            raise ValueError(f"pseudocount must be >= 0, got {self.pseudocount}")
        if self.zero_branch_len <= 0:
            raise ValueError(f"zero_branch_len must be > 0, got {self.zero_branch_len}")


@dataclass
class MDCatConfig:
    """MD-Cat 配置"""

    mdcat_bin: str = "mdcat"
    ncat: int = 50  # 速率类别数
    nrep: int = 100  # 随机初始化次数
    max_iter: int = 100  # 最大 EM 迭代次数
    backward_time: bool = True  # 反向时间（化石定年）
    as_date: bool = False  # 日期格式
    pseudocount: float = 0.01  # 分支长度伪计数（与直接导入路径及 wLogDate 保持一致）
    use_direct_import: bool = (
        False  # 默认使用 CLI；direct import 在大数据/蛋白比对上不稳定且不可中断
    )
    # CI bootstrap 次数。默认 0（不算置信区间），与上游 md_cat.py 的默认行为一致：
    # 实测 6 分类单元 / 500 位点的比对上，**一次** bootstrap 就要几分钟，
    # 旧默认 100 相当于一个数小时的隐式开销。需要 CI 时显式设 ci_nboots>0。
    ci_nboots: int = 0
    ci_plower: float = 0.025  # CI 下限百分位
    ci_pupper: float = 0.975  # CI 上限百分位
    annotate_level: int = 2  # 注释级别 1-3

    def __post_init__(self) -> None:
        """验证参数范围"""
        if self.ncat < 1:
            raise ValueError(f"ncat must be >= 1, got {self.ncat}")
        if self.nrep < 1:
            raise ValueError(f"nrep must be >= 1, got {self.nrep}")
        if self.max_iter < 1:
            raise ValueError(f"max_iter must be >= 1, got {self.max_iter}")
        if self.ci_nboots < 0:
            raise ValueError(f"ci_nboots must be >= 0, got {self.ci_nboots}")
        if not 0 <= self.ci_plower <= 1:
            raise ValueError(f"ci_plower must be between 0 and 1, got {self.ci_plower}")
        if not 0 <= self.ci_pupper <= 1:
            raise ValueError(f"ci_pupper must be between 0 and 1, got {self.ci_pupper}")
        if self.ci_plower >= self.ci_pupper:
            raise ValueError(
                f"ci_plower ({self.ci_plower}) must be < ci_pupper ({self.ci_pupper})"
            )
        if self.annotate_level not in [1, 2, 3]:
            raise ValueError(
                f"annotate_level must be 1, 2, or 3, got {self.annotate_level}"
            )


#: ``ToolConfig`` 里一层配置节点的并集。``set_cli_override`` 需要把
#: "common" 与七个方法子配置当同一个东西拿字段、做校验，所以它真实的
#: 类型就是这个并集（而不是 Any，也不是只能指一个的 ``object``）。
#: 必须放在八个 dataclass 都定义完之后。
ConfigSection = Union[
    CommonConfig,
    MCMCTreeConfig,
    R8sConfig,
    TreePLConfig,
    PATHd8Config,
    LSD2Config,
    WLogDateConfig,
    MDCatConfig,
]


@dataclass
class ToolConfig:
    """所有工具的配置容器"""

    common: CommonConfig = field(default_factory=CommonConfig)
    mcmctree: MCMCTreeConfig = field(default_factory=MCMCTreeConfig)
    r8s: R8sConfig = field(default_factory=R8sConfig)
    treepl: TreePLConfig = field(default_factory=TreePLConfig)
    pathd8: PATHd8Config = field(default_factory=PATHd8Config)
    lsd2: LSD2Config = field(default_factory=LSD2Config)
    wlogdate: WLogDateConfig = field(default_factory=WLogDateConfig)
    mdcat: MDCatConfig = field(default_factory=MDCatConfig)


@dataclass
class SoftwarePaths:
    """外部软件路径配置"""

    paml_path: Optional[str] = None
    mcmctree_bin: Optional[str] = None
    codeml_bin: Optional[str] = None
    baseml_bin: Optional[str] = None
    iqtree_bin: Optional[str] = None
    r8s_bin: Optional[str] = None
    pyr8s_bin: Optional[str] = None
    treepl_bin: Optional[str] = None
    pathd8_bin: Optional[str] = None
    mdcat_bin: Optional[str] = None
    wlogdate_bin: Optional[str] = None


class Configuration:
    """
    配置管理器

    实现三层配置合并机制：
    1. 内置默认值（最低优先级）
    2. YAML 配置文件（中等优先级）
    3. 命令行参数（最高优先级）
    """

    def __init__(self) -> None:
        self.logger = get_logger()
        self._defaults = ToolConfig()
        self._software_paths = SoftwarePaths()
        self._config_file: Optional[Path] = None
        self._cli_overrides: Dict[str, Any] = {}
        self._file_overrides: Dict[str, Dict[str, Any]] = (
            {}
        )  # tool_name -> {param: value}
        self._final_config: Optional[ToolConfig] = None

    @property
    def software_paths(self) -> SoftwarePaths:
        """获取外部软件路径配置"""
        return self._software_paths

    def set_software_path(self, tool: str, path: str) -> None:
        """设置外部软件路径"""
        if hasattr(self._software_paths, tool):
            setattr(self._software_paths, tool, path)
            self._final_config = None  # 使缓存失效
            self.logger.debug(f"Software path set: {tool} = {path}")

    def load_from_file(self, config_path: Path) -> "Configuration":
        """从 YAML 文件加载配置"""
        self._config_file = config_path

        if not config_path.exists():
            self.logger.warning(f"Config file not found: {config_path}")
            return self

        try:
            with open(config_path, "r") as f:
                data = yaml.safe_load(f)

            if not data:
                return self

            # 解析 software_paths 块
            software_paths = data.get("software_paths", {})
            for path_key, path_value in software_paths.items():
                if path_value and hasattr(self._software_paths, path_key):
                    setattr(self._software_paths, path_key, path_value)

            # 解析 software 块
            software = data.get("software", {})

            # 处理 common 配置
            common_config = software.get("common", {})
            if common_config:
                self._validate_and_merge("common", common_config)

            # 合并各工具配置
            for tool_name, tool_config in software.items():
                if tool_name == "common":
                    continue  # 已经处理过
                if hasattr(self._defaults, tool_name):
                    self._validate_and_merge(tool_name, tool_config)
                else:
                    self.logger.warning(
                        f"Unknown tool in config: {tool_name}",
                        category=UnknownParameterWarning,
                    )

            self.logger.info(f"Loaded configuration from {config_path}")

        except Exception as e:
            self.logger.error(f"Failed to load config file: {e}")
            raise

        return self

    def set_cli_override(self, tool: str, param: str, value: Any) -> None:
        """设置命令行参数覆盖"""
        import difflib

        # 验证工具名和参数名
        if tool == "common":
            key = f"common.{param}"
            tool_obj: ConfigSection = self._defaults.common
        else:
            candidate = getattr(self._defaults, tool, None)
            if candidate is None:
                raise ValueError(f"Unknown method '{tool}'")
            tool_obj = candidate
            key = f"{tool}.{param}"

        # 白名单校验
        valid_fields = getattr(tool_obj, "__dataclass_fields__", {})
        if param not in valid_fields:
            candidates = difflib.get_close_matches(param, list(valid_fields), n=3)
            hint = f" 是否想输入: {', '.join(candidates)}" if candidates else ""
            raise ValueError(f"Unknown parameter '{param}' for method '{tool}'.{hint}")

        self._cli_overrides[key] = value
        self._final_config = None  # 使缓存失效

    def get_config(self) -> ToolConfig:
        """获取最终合并后的配置"""
        if self._final_config is None:
            self._final_config = self._merge_configs()
        return self._final_config

    def _validate_and_merge(self, tool_name: str, config_dict: Dict[str, Any]) -> None:
        """验证文件配置并存储覆盖值（不修改 _defaults）

        验证字段名和类型后，将有效参数存储到 _file_overrides 中，
        在 _merge_configs 中统一应用。
        """
        tool_obj = getattr(self._defaults, tool_name)
        valid_fields = tool_obj.__dataclass_fields__

        overrides = {}
        for key, value in config_dict.items():
            if key in valid_fields:
                field_def = valid_fields[key]
                try:
                    coerced = self._coerce_value(value, field_def.type)
                    if coerced is not None:
                        # 通过构造新实例触发 __post_init__ 校验
                        current = {
                            f: getattr(tool_obj, f)
                            for f in tool_obj.__dataclass_fields__
                        }
                        current[key] = coerced
                        tool_obj.__class__(**current)
                        overrides[key] = value
                except (ValueError, TypeError) as e:
                    self.logger.warning(
                        f"Invalid value for {tool_name}.{key}: {value} ({e})"
                    )
            else:
                self.logger.warning(
                    f"Unknown parameter: {tool_name}.{key}",
                    category=UnknownParameterWarning,
                )

        if overrides:
            self._file_overrides[tool_name] = overrides

    @staticmethod
    def _coerce_value(value: Any, target_type: Type) -> Any:
        """根据目标类型转换配置值，支持 Optional[T] 与基本类型。"""
        if value is None:
            return None

        origin = getattr(target_type, "__origin__", None)
        if origin is typing.Union:
            args = target_type.__args__
            non_none_args = [a for a in args if a is not type(None)]
            if len(non_none_args) == 1:
                if value in (None, "None", "none", ""):
                    return None
                return Configuration._coerce_value(value, non_none_args[0])

        if target_type is bool:
            if isinstance(value, str):
                return value.lower() in ("true", "1", "yes", "on")
            return bool(value)
        if target_type is int:
            return int(value)
        if target_type is float:
            return float(value)
        if target_type is str:
            return str(value)
        return value

    def _set_typed_attr(
        self, tool_obj: Any, param_name: str, value: Any, source: str
    ) -> None:
        """类型安全地设置 dataclass 字段，触发 __post_init__ 校验。"""
        field_def = tool_obj.__dataclass_fields__.get(param_name)
        if field_def is None:
            setattr(tool_obj, param_name, value)
            return

        try:
            coerced = self._coerce_value(value, field_def.type)
        except (ValueError, TypeError) as e:
            self.logger.warning(
                f"Cannot convert {source}.{param_name} value {value!r} to {field_def.type}: {e}",
                category=UnknownParameterWarning,
            )
            return

        # 通过构造新实例触发 __post_init__ 校验
        try:
            current = {f: getattr(tool_obj, f) for f in tool_obj.__dataclass_fields__}
            current[param_name] = coerced
            validated = tool_obj.__class__(**current)
            # 取回 __post_init__ 规范化后的值（例如 R8sConfig.method 会把
            # "pl" 归一成 "PL"）：直接写回未规范化的原始值会让校验通过的
            # "形式"与实际生效的"取值"不一致。
            setattr(tool_obj, param_name, getattr(validated, param_name, coerced))
            self.logger.debug(
                f"{source}: {tool_obj.__class__.__name__}.{param_name} = "
                f"{getattr(tool_obj, param_name)!r}"
            )
        except (ValueError, TypeError) as e:
            self.logger.warning(
                f"Invalid value for {source}.{param_name}: {value!r} ({e})",
                category=UnknownParameterWarning,
            )

    def _merge_configs(self) -> ToolConfig:
        """合并三层配置：默认值 < 文件覆盖 < CLI 覆盖"""
        # 从默认值的深拷贝开始（防止任何覆盖污染 _defaults）
        final = copy.deepcopy(self._defaults)

        # 应用文件覆盖（第二层）
        for tool_name, overrides in self._file_overrides.items():
            if hasattr(final, tool_name):
                tool_obj = getattr(final, tool_name)
                for param_name, value in overrides.items():
                    if hasattr(tool_obj, param_name):
                        self._set_typed_attr(
                            tool_obj, param_name, value, f"File override {tool_name}"
                        )

        # 应用命令行覆盖（第三层，最高优先级）
        for key, value in self._cli_overrides.items():
            if value is None:
                continue

            parts = key.split(".")
            if len(parts) == 2:
                tool_name, param_name = parts
                if hasattr(final, tool_name):
                    tool_obj = getattr(final, tool_name)
                    if hasattr(tool_obj, param_name):
                        self._set_typed_attr(
                            tool_obj, param_name, value, f"CLI override {tool_name}"
                        )

        return final

    def save_runtime_config(self, output_path: Path) -> None:
        """保存运行时配置到文件"""
        config = self.get_config()

        data = {
            "software": {
                "common": asdict(config.common),
                "mcmctree": asdict(config.mcmctree),
                "r8s": asdict(config.r8s),
                "treepl": asdict(config.treepl),
                "pathd8": asdict(config.pathd8),
                "lsd2": asdict(config.lsd2),
                "wlogdate": asdict(config.wlogdate),
                "mdcat": asdict(config.mdcat),
            },
            "software_paths": {
                k: v for k, v in asdict(self._software_paths).items() if v is not None
            },
        }

        output_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = output_path.with_suffix(".tmp")
        with safe_writer(tmp_path) as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        import os

        os.chmod(tmp_path, 0o600)
        tmp_path.replace(output_path)

        self.logger.info(f"Saved runtime configuration to {output_path}")

    def save_to_file(self, output_path: Path) -> None:
        """保存默认配置到文件 (--generate-config 使用)"""
        self.create_default_config_file(output_path)

    @staticmethod
    def create_default_config_file(output_path: Path) -> None:
        """创建默认配置文件"""
        config = ToolConfig()
        software_paths = SoftwarePaths()

        data = {
            "software": {
                "common": asdict(config.common),
                "mcmctree": asdict(config.mcmctree),
                "r8s": asdict(config.r8s),
                "treepl": asdict(config.treepl),
                "pathd8": asdict(config.pathd8),
                "lsd2": asdict(config.lsd2),
                "wlogdate": asdict(config.wlogdate),
                "mdcat": asdict(config.mdcat),
            },
            "software_paths": {k: v for k, v in asdict(software_paths).items()},
        }

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with safe_writer(output_path) as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
