"""
PhyloDater: 多软件并行系统发育定年平台

PhyloDater 是一个统一的系统发育定年平台，整合多种定年软件，
提供一致的接口和标准化的结果输出。

主要特性：
- 支持 MCMCTree、r8s、treePL、PATHd8、LSD2、wLogDate、MD-Cat 等多种定年软件
- 统一的年龄约束类型体系
- 自动化的 MRCA 节点定位
- 标准化的结果比较和报告生成
- 可重复性元数据追踪

示例用法：
    >>> from phylodater import PhylogeneticTree, DatingMethodRegistry
    >>> tree = PhylogeneticTree.from_file("tree.nwk")
    >>> adapter = DatingMethodRegistry.create("pathd8", config, output_dir)
    >>> adapter.prepare_inputs(tree, calibrations)
    >>> adapter.execute()
    >>> result = adapter.parse_results()
"""

from pathlib import Path
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:  # 只供标注使用；services 在模块尾部才导入，避免启动期成环
    from .services.tree_validator import ValidationReport

# 版本信息：从 _version.py 读取（静态版本号，无 git 派生的迭代串）。
try:
    from ._version import __version__, __version_tuple__
except ImportError:  # pragma: no cover - 未构建/未随包安装时的兜底
    __version__ = "0.1.0"
    __version_tuple__ = (0, 1, 0)

__git_hash__ = None

__author__ = "Zichao Zeng (曾子超)"
__license__ = "MIT License"
__update_date__ = "2026-09-17"

MOSAIC_BANNER = r"""
==========================================================================================

   ██████╗ ██╗  ██╗██╗   ██╗██╗       ██████╗  ██████╗  █████╗ ████████╗███████╗██████╗
   ██╔══██╗██║  ██║╚██╗ ██╔╝██║      ██╔═══██╗ ██╔══██╗██╔══██╗╚══██╔══╝██╔════╝██╔══██╗
   ██████╔╝███████║ ╚████╔╝ ██║      ██║   ██║ ██║  ██║███████║   ██║   █████╗  ██████╔╝
   ██╔═══╝ ██╔══██║  ╚██╔╝  ██║      ██║   ██║ ██║  ██║██╔══██║   ██║   ██╔══╝  ██╔══██╗
   ██║     ██║  ██║   ██║   ███████╗ ╚██████╔╝ ██████╔╝██║  ██║   ██║   ███████╗██║  ██║
   ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝  ╚═╝

                                  PhyloDater  v{__version__}
                                  {__license__}
                                  Last updated: {__update_date__}

==========================================================================================

  Dependencies & Licenses:
    BioPython  BSD-3-Clause    DendroPy   BSD-3-Clause    NumPy      BSD-3-Clause
    ETE3       GPL-3.0         PyYAML     MIT License     Pandas     BSD-3-Clause
    Matplotlib PSF-based       SciPy      BSD-3-Clause    psutil     BSD-3-Clause

==========================================================================================
"""


def get_banner() -> str:
    """获取ASCII马赛克图标"""
    return (
        MOSAIC_BANNER.replace("{__version__}", __version__)
        .replace("{__license__}", __license__)
        .replace("{__update_date__}", __update_date__)
    )


def show_banner() -> None:
    """显示ASCII马赛克图标"""
    print(get_banner())


from .adapters import MDCatMethod  # noqa: E402
from .adapters import (  # noqa: E402
    LSD2Method,
    MCMCTreeMethod,
    PATHd8Method,
    R8sPyr8sMethod,
    TreePLMethod,
    WLogDateMethod,
)
from .core import ComparisonReporter  # noqa: E402
from .core import (  # noqa: E402
    CalibrationError,
    DatingMethod,
    DatingMethodRegistry,
    DatingPipeline,
    ExecutionError,
    InputSizeLimitError,
    MonophylyError,
    NegativeBranchLengthError,
    ParallelPipeline,
    PhyloDaterError,
    PhyloFormatError,
    PipelineConfig,
    SemanticDegradationWarning,
    TaxonomyConflictError,
    TreeValidationError,
    ValidationError,
)
from .infrastructure import LogLevel  # noqa: E402
from .infrastructure import (  # noqa: E402
    Configuration,
    PhyloDaterLogger,
    ToolConfig,
    get_logger,
    setup_logging,
)
from .models import AgeConstraint  # noqa: E402
from .models import (  # noqa: E402
    ADAPTER_CAPABILITIES,
    CalibrationPoint,
    CIType,
    DatingResult,
    FixedAgeConstraint,
    FossilConfidence,
    FossilMetadata,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    NodeAgeEstimate,
    NodeDefinition,
    PhylogeneticTree,
    SkewNormalConstraint,
    SkewTConstraint,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)
from .services import CalibrationLoader  # noqa: E402
from .services import (  # noqa: E402
    SHORT_NAME_LENGTH_LIMITS,
    AlignmentFormat,
    CalibrationResolver,
    DeepValidator,
    InputValidator,
    NameCollisionError,
    NameMappingManager,
    NameShortenMode,
    RootingMethod,
    RootingResult,
    TaxonomyConsistencyReport,
    TaxonomyLoadError,
    TaxonomyLoadReport,
    TaxonomyParser,
    TreeFormat,
    TreeRootingService,
    TreeValidator,
    create_default_calibration_config,
    resolve_short_name,
    sanitize_names_uniquely,
)
from .viz import DatingFigure, GeoPlotter, NodeCoordinate, TreePlotter  # noqa: E402


def parse_taxonomy(label: str, mode: str = "reverse") -> dict:
    """解析分类学名称，返回标准字典。

    Args:
        label: 末端节点标签或分类学字符串
        mode: 解析模式（reverse/greedy/segment）

    Returns:
        包含 domain/phylum/class/order/family/genus/species 的字典
    """
    from .services.taxonomy_parser import TaxonomyParser

    parser = TaxonomyParser(delimiter_mode=mode)
    result = parser.parse(label)
    return result or {}


def is_monophyletic(
    tree: Union[PhylogeneticTree, str, Path],
    taxon_label: str,
    rooted: bool = True,
) -> bool:
    """判断给定分类标签在树中是否构成单系群。

    Args:
        tree: 树文件路径或 PhylogeneticTree 对象
        taxon_label: 分类标签（如 'Cyanobacteriota' 或 'LUCA'）
        rooted: 是否已定根

    Returns:
        True 如果单系，False 否则

    Raises:
        ValueError: 如果树未定根或标签无法解析到任何叶节点
    """

    from .models.tree import PhylogeneticTree
    from .services.calibration_resolver import CalibrationResolver
    from .services.taxonomy_parser import TaxonomyParser

    if isinstance(tree, (str, Path)):
        tree = PhylogeneticTree.from_file(Path(tree))
    elif not isinstance(tree, PhylogeneticTree):
        raise TypeError("tree must be a path or PhylogeneticTree")

    if not rooted:
        raise ValueError(
            "is_monophyletic requires a rooted tree; "
            "call tree.get_rooting_info() first."
        )

    resolver = CalibrationResolver(tree, TaxonomyParser())
    result = resolver.resolve(taxon_label)
    if result is None or not result.resolved_taxa:
        raise ValueError(
            f"taxon_label '{taxon_label}' could not be resolved to any tip"
        )
    return bool(tree.is_monophyletic(result.resolved_taxa))


def load_tree(path: Union[str, Path], validate: bool = True) -> PhylogeneticTree:
    """加载树文件。

    Args:
        path: 树文件路径
        validate: 是否验证

    Returns:
        PhylogeneticTree 对象
    """

    from .models.tree import PhylogeneticTree

    if isinstance(path, str):
        path = Path(path)
    tree = PhylogeneticTree.from_file(path)
    if validate:
        # 执行基本验证
        if tree.num_tips == 0:
            raise ValueError("Tree contains no terminal nodes")
    return tree


def cross_validate(
    tree_path: Union[str, Path],
    seq_path: Union[str, Path],
    strict: bool = True,
) -> "ValidationReport":
    """交叉验证树与序列文件的标签一致性。

    Args:
        tree_path: 树文件路径
        seq_path: 序列文件路径
        strict: True 则不匹配时报错，False 则仅警告

    Returns:
        ValidationReport 对象
    """

    from .services.tree_validator import TreeValidator

    if isinstance(tree_path, str):
        tree_path = Path(tree_path)
    if isinstance(seq_path, str):
        seq_path = Path(seq_path)
    tree_obj = PhylogeneticTree.from_file(tree_path)
    validator = TreeValidator()
    report = validator.validate_sequence_consistency(seq_path, tree_obj)
    if strict and not report.is_valid:
        raise TreeValidationError(f"Cross-validation failed: {report.errors}")
    return report


__all__ = [
    # 版本信息
    "__version__",
    "__author__",
    # Banner
    "MOSAIC_BANNER",
    "show_banner",
    "get_banner",
    # 约束类型
    "AgeConstraint",
    "FixedAgeConstraint",
    "UniformAgeConstraint",
    "SoftLowerBoundConstraint",
    "MaximumAgeConstraint",
    "SoftBoundsConstraint",
    "GammaPriorConstraint",
    "SkewNormalConstraint",
    "SkewTConstraint",
    # 校准相关
    "FossilMetadata",
    "FossilConfidence",
    "CalibrationPoint",
    "NodeDefinition",
    # 结果相关
    "NodeAgeEstimate",
    "DatingResult",
    "CIType",
    "ADAPTER_CAPABILITIES",
    # 树
    "PhylogeneticTree",
    # 异常
    "PhyloDaterError",
    "ValidationError",
    "TreeValidationError",
    "CalibrationError",
    "ExecutionError",
    "SemanticDegradationWarning",
    "PhyloFormatError",
    "TaxonomyConflictError",
    "MonophylyError",
    "NegativeBranchLengthError",
    "InputSizeLimitError",
    # 接口
    "DatingMethod",
    "DatingMethodRegistry",
    # 流水线
    "DatingPipeline",
    "ParallelPipeline",
    "PipelineConfig",
    # 报告
    "ComparisonReporter",
    # 基础设施
    "get_logger",
    "PhyloDaterLogger",
    "LogLevel",
    "setup_logging",
    "Configuration",
    "ToolConfig",
    # 服务
    "NameMappingManager",
    "NameShortenMode",
    # 短名唯一性（审阅项 B-21）
    "NameCollisionError",
    "resolve_short_name",
    "sanitize_names_uniquely",
    "SHORT_NAME_LENGTH_LIMITS",
    "TaxonomyParser",
    # 分类表加载/一致性核对（审阅项 B-26）
    "TaxonomyLoadError",
    "TaxonomyLoadReport",
    "TaxonomyConsistencyReport",
    "TreeValidator",
    "TreeRootingService",
    "RootingMethod",
    "RootingResult",
    "CalibrationResolver",
    "CalibrationLoader",
    "create_default_calibration_config",
    # 输入验证
    "InputValidator",
    "TreeFormat",
    "AlignmentFormat",
    # 深度验证
    "DeepValidator",
    # 适配器
    "LSD2Method",
    "MCMCTreeMethod",
    "PATHd8Method",
    "TreePLMethod",
    "R8sPyr8sMethod",
    "WLogDateMethod",
    "MDCatMethod",
    # 可视化
    "DatingFigure",
    "TreePlotter",
    "GeoPlotter",
    "NodeCoordinate",
    # 库模式API
    "parse_taxonomy",
    "is_monophyletic",
    "load_tree",
    "cross_validate",
]
