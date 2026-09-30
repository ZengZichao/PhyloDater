"""
PhyloDater 核心层

包含异常层次结构和 DatingMethod 接口
"""

from .comparison_reporter import ComparisonReporter
from .exceptions import (
    AlignmentValidationError,
    CalibrationConflictError,
    CalibrationError,
    CalibrationResolutionError,
    ConfigurationError,
    ConsistencyError,
    ConvergenceWarning,
    EnvironmentWarning,
    ExecutionError,
    ExecutionTimeoutError,
    InputSizeLimitError,
    MonophylyError,
    NegativeBranchLengthError,
    PhyloDaterError,
    PhyloFormatError,
    ResourceWarning,
    ResultParsingError,
    SemanticDegradationWarning,
    TaxonomyConflictError,
    TreeValidationError,
    UnknownMethodError,
    UnknownParameterWarning,
    ValidationError,
)
from .method_interface import DatingMethod, DatingMethodRegistry
from .pipeline import DatingPipeline, ParallelPipeline, PipelineConfig

__all__ = [
    # 异常
    "PhyloDaterError",
    "ValidationError",
    "TreeValidationError",
    "AlignmentValidationError",
    "ConsistencyError",
    "CalibrationError",
    "CalibrationConflictError",
    "CalibrationResolutionError",
    "ExecutionError",
    "ExecutionTimeoutError",
    "ConvergenceWarning",
    "ResultParsingError",
    "SemanticDegradationWarning",
    "EnvironmentWarning",
    "ResourceWarning",
    "UnknownParameterWarning",
    "UnknownMethodError",
    "PhyloFormatError",
    "TaxonomyConflictError",
    "MonophylyError",
    "NegativeBranchLengthError",
    "InputSizeLimitError",
    "ConfigurationError",
    # 接口
    "DatingMethod",
    "DatingMethodRegistry",
    # 流水线
    "DatingPipeline",
    "ParallelPipeline",
    "PipelineConfig",
    # 报告
    "ComparisonReporter",
]
