"""
校准点类

定义校准点及相关元数据结构
"""

import math
import warnings
from dataclasses import asdict, dataclass, field, is_dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

from .constraints import MAX_GEOLOGICAL_AGE_MA, AgeConstraint


class FossilConfidence(Enum):
    """化石鉴定置信度"""

    CONFIRMED = "confirmed"  # 已确认
    PROBABLE = "probable"  # 可能
    TENTATIVE = "tentative"  # 暂定


def _validate_stratigraphic_age(value: object, name: str) -> Optional[float]:
    """校验地层年龄字段（审阅项 C-6/B-19）。

    允许 ``None``（未填写）；否则必须是**有限**、非负、不超过地球年龄的数值。
    NaN/inf 一律拒绝：所有大小比较都对 NaN 恒为 False，静默收下等于让下游的
    "min ≤ max"检查形同虚设。
    """
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a number or None, got {type(value).__name__}")
    if not math.isfinite(value):
        raise ValueError(
            f"{name} ({value}) must be a finite number; NaN and infinity are "
            f"rejected because every age comparison silently passes"
        )
    if value < 0:
        raise ValueError(f"{name} ({value}) must be non-negative")
    if value > MAX_GEOLOGICAL_AGE_MA:
        raise ValueError(
            f"{name} ({value} Ma) exceeds geological reasonable bound "
            f"({MAX_GEOLOGICAL_AGE_MA} Ma - Earth's age)"
        )
    return float(value)


@dataclass
class FossilMetadata:
    """
    结构化化石证据信息

    替代原始的非结构化 dict，确保可重复性报告满足期刊披露要求

    地层年龄方向约定（务必读完整）：``min_stratigraphic_age`` 是**地层可能的最年
    轻年龄**（youngest possible age of the unit），``max_stratigraphic_age`` 是**最
    老可能年龄**，因此恒有 ``min ≤ max``。化石通常给出分歧年龄的**下限**。
    "应取 min 还是 max 作为保守下界"这一方向约定尚未与 YAML 文档/解析层统一核对
    （审阅报告 §七 E-1 待作者核实），所以本类只做区间内部的交叉校验，不推断方向。
    """

    fossil_name: Optional[str] = None  # 化石名称
    formation: Optional[str] = None  # 化石层位/地层单元
    geological_period: Optional[str] = None  # 地质时代
    min_stratigraphic_age: Optional[float] = None  # 地层最小年龄（最年轻可能年龄，Ma）
    max_stratigraphic_age: Optional[float] = None  # 地层最大年龄（最老可能年龄，Ma）
    reference: Optional[str] = None  # 文献引用（建议 DOI）
    confidence: Optional[FossilConfidence] = None  # 化石鉴定置信度
    notes: Optional[str] = None  # 自由文本备注

    def __post_init__(self) -> None:
        # 交叉校验（审阅项 C-6）：两个字段此前无任何大小关系/合法性检查，
        # min > max 会被静默接受并写进可复现性报告。
        self.min_stratigraphic_age = _validate_stratigraphic_age(
            self.min_stratigraphic_age, "min_stratigraphic_age"
        )
        self.max_stratigraphic_age = _validate_stratigraphic_age(
            self.max_stratigraphic_age, "max_stratigraphic_age"
        )
        if (
            self.min_stratigraphic_age is not None
            and self.max_stratigraphic_age is not None
            and self.min_stratigraphic_age > self.max_stratigraphic_age
        ):
            raise ValueError(
                f"min_stratigraphic_age ({self.min_stratigraphic_age} Ma) must be <= "
                f"max_stratigraphic_age ({self.max_stratigraphic_age} Ma): the fields "
                f"are the youngest/oldest possible ages of the fossil-bearing unit, so "
                f"a reversed pair silently inverts the reported stratigraphic range. "
                f"Values are never swapped automatically."
            )

    def to_dict(self) -> Dict[str, Any]:
        """转换为字典"""
        return {
            "fossil_name": self.fossil_name,
            "formation": self.formation,
            "geological_period": self.geological_period,
            "min_stratigraphic_age": self.min_stratigraphic_age,
            "max_stratigraphic_age": self.max_stratigraphic_age,
            "reference": self.reference,
            "confidence": self.confidence.value if self.confidence else None,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "FossilMetadata":
        """从字典创建（交叉校验见 :meth:`__post_init__`，审阅项 C-6）"""
        confidence = None
        if data.get("confidence"):
            confidence = FossilConfidence(data["confidence"])

        return cls(
            fossil_name=data.get("fossil_name"),
            formation=data.get("formation"),
            geological_period=data.get("geological_period"),
            min_stratigraphic_age=data.get("min_stratigraphic_age"),
            max_stratigraphic_age=data.get("max_stratigraphic_age"),
            reference=data.get("reference"),
            confidence=confidence,
            notes=data.get("notes"),
        )


@dataclass
class CalibrationPoint:
    """
    校准点

    系统的中心领域对象，整合各脚本中分散的校准点结构
    """

    name: str  # 生物学标识（如 LUCA、Cyanobacteriota）
    age_constraint: Optional[AgeConstraint] = (
        None  # 年龄约束对象（可选，解析后由调用方设置）
    )
    resolved_taxa: List[str] = field(default_factory=list)  # 解析的具体叶节点列表
    mrca_leaf_pair: Optional[tuple] = None  # 用于 PATHd8/treePL/r8s 的两叶节点对
    is_root_node: bool = False  # 是否定位于树根
    fossil_metadata: Optional[FossilMetadata] = None  # 结构化化石证据

    def __post_init__(self) -> None:
        # mrca_leaf_pair 不再从 resolved_taxa 自动设置，
        # 因为 (resolved_taxa[0], resolved_taxa[-1]) 是任意选取，不代表系统发育距离最远的一对。
        # 正确的 mrca_leaf_pair 应由 CalibrationResolver 通过 tree.get_mrca_terminals() 计算，
        # 或由用户在 YAML 中通过 mrca_pair 显式指定。

        # 防御性检查：如果校准点已解析但无年龄约束，发出警告
        if self.age_constraint is None and self.resolved_taxa:
            warnings.warn(
                f"CalibrationPoint '{self.name}' has resolved taxa but no age_constraint. "
                f"Callers must set age_constraint before using this calibration point in analysis.",
                UserWarning,
                stacklevel=2,
            )

        self._check_fossil_metadata_consistency()

    def _constraint_age_interval(self) -> Optional[tuple]:
        """约束隐含的年龄区间 [lo, hi]（Ma），无法确定界时对应分量为 None。"""
        from .constraints import (
            FixedAgeConstraint,
            MaximumAgeConstraint,
            SoftBoundsConstraint,
            SoftLowerBoundConstraint,
            UniformAgeConstraint,
        )

        constraint = self.age_constraint
        if constraint is None:
            return None
        if isinstance(constraint, FixedAgeConstraint):
            return (constraint.fixed_age, constraint.fixed_age)
        if isinstance(constraint, (UniformAgeConstraint, SoftBoundsConstraint)):
            return (constraint.min_age, constraint.max_age)
        if isinstance(constraint, MaximumAgeConstraint):
            return (None, constraint.max_age)
        if isinstance(constraint, SoftLowerBoundConstraint):
            return (constraint.min_age, None)
        return None

    def _check_fossil_metadata_consistency(self) -> None:
        """化石地层区间与校准区间的交叉校验（审阅项 C-6）。

        只检查在**任一方向约定下都成立**的事实：分歧年龄不可能比它所嵌入的地层还
        年轻。若校准允许的最老年龄 ``hi`` 低于地层可能的最年轻年龄 ``min_s``，说明
        两者毫无重叠——通常意味着 min/max 写反、单位写错或挂错了化石。这里发警告而
        非报错：把"地层的 min 还是 max 才是保守下界"这一方向约定钉死属于审阅报告
        §七 E-1 的作者核实事项，模型层不替作者猜。
        """
        metadata = self.fossil_metadata
        if metadata is None or metadata.min_stratigraphic_age is None:
            return
        interval = self._constraint_age_interval()
        if interval is None:
            return
        _lower, upper = interval
        if upper is None:
            return
        if upper < metadata.min_stratigraphic_age:
            warnings.warn(
                f"CalibrationPoint '{self.name}': the constraint's maximum age "
                f"({upper} Ma) is younger than the youngest possible age of the "
                f"fossil-bearing unit "
                f"(min_stratigraphic_age={metadata.min_stratigraphic_age} Ma, "
                f"max_stratigraphic_age={metadata.max_stratigraphic_age} Ma), so the "
                f"two intervals do not overlap at all. Check for reversed "
                f"min/max stratigraphic ages, a unit error (Ma vs ka/Ga) or a "
                f"mismatched fossil (see review item C-6/E-1).",
                UserWarning,
                stacklevel=3,
            )

    @property
    def has_upper_bound(self) -> bool:
        """是否存在上界约束

        审阅项 C-7：固定点 ``FixedAgeConstraint`` 此前只出现在 ``has_lower_bound``
        里，但一个点年龄在数学上同时是上界与下界——漏掉它会让任何用上界做单调性/
        冲突检查的下游完全看不见固定点校准。
        """
        from .constraints import (
            FixedAgeConstraint,
            MaximumAgeConstraint,
            SoftBoundsConstraint,
            UniformAgeConstraint,
        )

        return isinstance(
            self.age_constraint,
            (
                UniformAgeConstraint,
                MaximumAgeConstraint,
                SoftBoundsConstraint,
                FixedAgeConstraint,
            ),
        )

    @property
    def has_lower_bound(self) -> bool:
        """是否存在下界约束（固定点同时是上下界，见 ``has_upper_bound``）"""
        from .constraints import (
            FixedAgeConstraint,
            SoftBoundsConstraint,
            SoftLowerBoundConstraint,
            UniformAgeConstraint,
        )

        return isinstance(
            self.age_constraint,
            (
                UniformAgeConstraint,
                SoftLowerBoundConstraint,
                SoftBoundsConstraint,
                FixedAgeConstraint,
            ),
        )

    def to_dict(self) -> Dict[str, Any]:
        """转换为字典"""
        age_constraint_dict = None
        if self.age_constraint is not None:
            constraint = self.age_constraint
            # ``AgeConstraint`` 是 ABC，它的 8 个具体子类全部是 dataclass；
            # 但 ``dataclasses.asdict`` 只接受实例（不接受类对象），所以这里
            # 必须把这两个前提当场核住，而不是期待基类标注自动满足它们。
            if not is_dataclass(constraint) or isinstance(constraint, type):
                raise TypeError(
                    f"age_constraint {type(constraint).__name__} is not a dataclass "
                    "instance; to_dict()/from_dict() symmetry requires one"
                )
            age_constraint_dict = {
                "type": constraint.__class__.__name__,
                # 使用 dataclasses.asdict 仅序列化声明的字段（白名单语义），
                # 避免 __dict__ 泄漏私有/非字段属性，保证与 from_dict 对称。
                "params": asdict(constraint),
            }
        return {
            "name": self.name,
            "age_constraint": age_constraint_dict,
            "resolved_taxa": self.resolved_taxa,
            "mrca_leaf_pair": self.mrca_leaf_pair,
            "is_root_node": self.is_root_node,
            "fossil_metadata": (
                self.fossil_metadata.to_dict() if self.fossil_metadata else None
            ),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CalibrationPoint":
        """从字典创建"""
        # 动态创建约束对象
        age_constraint = None
        constraint_data = data.get("age_constraint")
        if constraint_data is not None:
            constraint_type = constraint_data["type"]
            constraint_params = constraint_data["params"]

            # 显式映射约束类型（避免使用 globals()）
            from .constraints import (
                FixedAgeConstraint,
                GammaPriorConstraint,
                MaximumAgeConstraint,
                SkewNormalConstraint,
                SkewTConstraint,
                SoftBoundsConstraint,
                SoftLowerBoundConstraint,
                UniformAgeConstraint,
            )

            _CONSTRAINT_CLASSES = {
                "FixedAgeConstraint": FixedAgeConstraint,
                "UniformAgeConstraint": UniformAgeConstraint,
                "SoftLowerBoundConstraint": SoftLowerBoundConstraint,
                "MaximumAgeConstraint": MaximumAgeConstraint,
                "SoftBoundsConstraint": SoftBoundsConstraint,
                "GammaPriorConstraint": GammaPriorConstraint,
                "SkewNormalConstraint": SkewNormalConstraint,
                "SkewTConstraint": SkewTConstraint,
            }

            constraint_class = _CONSTRAINT_CLASSES.get(constraint_type)
            if constraint_class is None:
                raise ValueError(f"Unknown constraint type: {constraint_type}")
            age_constraint = constraint_class(**constraint_params)

        # 解析化石元数据
        fossil_metadata = None
        if data.get("fossil_metadata"):
            fossil_metadata = FossilMetadata.from_dict(data["fossil_metadata"])

        return cls(
            name=data["name"],
            age_constraint=age_constraint,
            resolved_taxa=data.get("resolved_taxa", []),
            mrca_leaf_pair=(
                tuple(data["mrca_leaf_pair"]) if data.get("mrca_leaf_pair") else None
            ),
            is_root_node=data.get("is_root_node", False),
            fossil_metadata=fossil_metadata,
        )


@dataclass
class NodeDefinition:
    """
    节点定义

    用于 YAML 配置中定义校准节点的位置
    """

    definition_type: str = "auto"  # auto, mrca, mrca_pair
    taxa: Optional[List[str]] = None  # 用于 mrca 类型
    mrca_pair: Optional[List[str]] = None  # 用于 mrca_pair 类型

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "NodeDefinition":
        """从字典创建"""
        return cls(
            definition_type=data.get("type", "auto"),
            taxa=data.get("taxa"),
            mrca_pair=data.get("mrca_pair"),
        )
