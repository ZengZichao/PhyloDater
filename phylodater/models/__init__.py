"""
PhyloDater 领域模型

包含年龄约束、校准点、定年结果、系统发育树等核心领域对象
"""

from .calibration import (
    CalibrationPoint,
    FossilConfidence,
    FossilMetadata,
    NodeDefinition,
)
from .constraints import (
    AgeConstraint,
    FixedAgeConstraint,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    SkewNormalConstraint,
    SkewTConstraint,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)
from .results import ADAPTER_CAPABILITIES, CIType, DatingResult, NodeAgeEstimate
from .tree import PhylogeneticTree, strip_leading_newick_comments

__all__ = [
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
    "strip_leading_newick_comments",
]
