"""
PhyloDater 服务层

包含名称映射、分类学解析、树验证、校准解析等服务
"""

from .auto_calibrator import (
    AutoCalibration,
    AutoCalibrator,
    CalibratedTree,
    parse_auto_calibrate_args,
)
from .calibration_loader import CalibrationLoader, create_default_calibration_config
from .calibration_resolver import (
    CalibrationResolver,
    ResolutionResult,
    resolve_special_nodes,
)
from .deep_validator import (
    DeepValidator,
    SequenceAlphabet,
    SequenceFormat,
    SequenceValidationDetail,
    TreeValidationDetail,
)
from .input_validator import (
    AlignmentFormat,
    InputValidator,
    TreeFormat,
    ValidationResult,
)
from .name_mapping import (
    SHORT_NAME_LENGTH_LIMITS,
    NameCollisionError,
    NameMappingEntry,
    NameMappingManager,
    NameShortenMode,
    resolve_short_name,
    sanitize_names_uniquely,
)
from .taxonomy_parser import (
    TaxonomyConsistencyReport,
    TaxonomyLoadError,
    TaxonomyLoadReport,
    TaxonomyParser,
)
from .tree_rooting import RootingMethod, RootingResult, TreeRootingService
from .tree_validator import TreeValidator, ValidationReport

__all__ = [
    # 名称映射
    "NameMappingManager",
    "NameMappingEntry",
    "NameShortenMode",
    # 短名唯一性（审阅项 B-21：共享登录号的两个序列不得共用一个短名）
    "NameCollisionError",
    "resolve_short_name",
    "sanitize_names_uniquely",
    "SHORT_NAME_LENGTH_LIMITS",
    # 分类学解析
    "TaxonomyParser",
    # 分类表加载/一致性核对（审阅项 B-26：整表加载失败不得伪装成"0 条但已核对"）
    "TaxonomyLoadError",
    "TaxonomyLoadReport",
    "TaxonomyConsistencyReport",
    # 树验证
    "TreeValidator",
    "ValidationReport",
    # 校准解析
    "CalibrationResolver",
    "ResolutionResult",
    "resolve_special_nodes",
    # 校准加载
    "CalibrationLoader",
    "create_default_calibration_config",
    # 树定根
    "TreeRootingService",
    "RootingMethod",
    "RootingResult",
    # 输入验证
    "InputValidator",
    "TreeFormat",
    "AlignmentFormat",
    "ValidationResult",
    # 深度验证
    "DeepValidator",
    "TreeValidationDetail",
    "SequenceValidationDetail",
    "SequenceAlphabet",
    "SequenceFormat",
    # 自动校准
    "AutoCalibrator",
    "AutoCalibration",
    "CalibratedTree",
    "parse_auto_calibrate_args",
]
