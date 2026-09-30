"""
定年结果类

定义标准化的定年结果数据结构
"""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


class CIType(Enum):
    """置信区间类型"""

    HPD95 = "HPD95"  # 最高后验密度 95% 区间（贝叶斯）
    CI95 = "CI95"  # 95% 置信区间
    RANGE = "RANGE"  # 范围（如均值±标准差）
    NONE = None  # 无区间


@dataclass
class NodeAgeEstimate:
    """
    节点年龄估计

    单个校准节点的年龄估计结果
    """

    mean_age: float  # 均值年龄（Ma）
    median_age: Optional[float] = None  # 中位年龄（Ma，贝叶斯方法提供）
    ci_lower: Optional[float] = None  # 区间下界（Ma）
    ci_upper: Optional[float] = None  # 区间上界（Ma）
    ci_type: CIType = CIType.NONE  # 区间类型

    def __post_init__(self) -> None:
        # 确保 ci_type 是 CIType 枚举
        if isinstance(self.ci_type, str):
            self.ci_type = CIType(self.ci_type)

    def to_dict(self) -> Dict[str, Any]:
        """转换为字典"""
        return {
            "mean_age": self.mean_age,
            "median_age": self.median_age,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "ci_type": self.ci_type.value if self.ci_type else None,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "NodeAgeEstimate":
        """从字典创建"""
        return cls(
            mean_age=data["mean_age"],
            median_age=data.get("median_age"),
            ci_lower=data.get("ci_lower"),
            ci_upper=data.get("ci_upper"),
            ci_type=CIType(data.get("ci_type")) if data.get("ci_type") else CIType.NONE,
        )


@dataclass
class DatingResult:
    """
    定年结果

    单个定年方法的完整结果
    """

    method_name: str  # 适配器标识符
    run_id: str  # 运行唯一标识
    dated_tree_newick: str  # 带年龄注释的 Newick
    node_ages: Dict[str, NodeAgeEstimate] = field(default_factory=dict)  # 节点年龄估计
    raw_output_path: Optional[Path] = None  # 原始输出目录
    execution_seconds: float = 0.0  # 执行耗时
    warnings: List[Any] = field(default_factory=list)  # 语义降级警告（延迟导入类型）
    is_converged: Optional[bool] = None  # 是否收敛（贝叶斯方法）
    metadata: Dict[str, Any] = field(default_factory=dict)  # 额外元数据

    def get_age(self, node_name: str) -> Optional[float]:
        """获取指定节点的均值年龄"""
        estimate = self.node_ages.get(node_name)
        return estimate.mean_age if estimate else None

    def get_ci(self, node_name: str) -> Optional[Tuple[float, float]]:
        """获取指定节点的置信区间"""
        estimate = self.node_ages.get(node_name)
        if estimate and estimate.ci_lower is not None and estimate.ci_upper is not None:
            return (estimate.ci_lower, estimate.ci_upper)
        return None

    def add_warning(self, warning: Any) -> None:
        """添加语义降级警告"""
        self.warnings.append(warning)

    def to_dict(self) -> Dict[str, Any]:
        """转换为字典（用于序列化）"""
        return {
            "method_name": self.method_name,
            "run_id": self.run_id,
            "dated_tree_newick": self.dated_tree_newick,
            "node_ages": {
                name: estimate.to_dict() for name, estimate in self.node_ages.items()
            },
            "raw_output_path": (
                str(self.raw_output_path) if self.raw_output_path else None
            ),
            "execution_seconds": self.execution_seconds,
            "is_converged": self.is_converged,
            "warnings": [str(w) for w in self.warnings] if self.warnings else [],
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DatingResult":
        """从字典创建"""
        return cls(
            method_name=data.get("method_name", "unknown"),
            run_id=data.get("run_id", ""),
            dated_tree_newick=data.get("dated_tree_newick", ""),
            node_ages={
                name: NodeAgeEstimate.from_dict(estimate_data)
                for name, estimate_data in data.get("node_ages", {}).items()
            },
            raw_output_path=(
                Path(data["raw_output_path"]) if data.get("raw_output_path") else None
            ),
            execution_seconds=data.get("execution_seconds", 0.0),
            is_converged=data.get("is_converged"),
            warnings=data.get("warnings", []),
            metadata=data.get("metadata", {}),
        )

    def format_age_string(self, node_name: str) -> str:
        """
        格式化年龄字符串用于比较表格

        格式：
        - 点估计：{mean_age}
        - 带区间：{mean_age} [{ci_lower},{ci_upper}]_{ci_type}
        """
        estimate = self.node_ages.get(node_name)
        if not estimate:
            return "NA"

        if estimate.ci_lower is not None and estimate.ci_upper is not None:
            return (
                f"{estimate.mean_age:.1f} "
                f"[{estimate.ci_lower:.1f},{estimate.ci_upper:.1f}]_"
                f"{estimate.ci_type.value}"
            )
        else:
            return f"{estimate.mean_age:.1f}"


# 各适配器能提供的字段一览（文档参考）
ADAPTER_CAPABILITIES = {
    "mcmctree": {
        "mean_age": True,
        "median_age": True,
        "ci_lower": True,
        "ci_upper": True,
        "ci_type": CIType.HPD95,
        "is_converged": True,
    },
    "treepl": {
        "mean_age": True,
        "median_age": False,
        "ci_lower": False,
        "ci_upper": False,
        "ci_type": CIType.NONE,
        "is_converged": False,
    },
    "r8s": {
        "mean_age": True,
        "median_age": False,
        "ci_lower": False,
        "ci_upper": False,
        "ci_type": CIType.NONE,
        "is_converged": False,
    },
    "pathd8": {
        "mean_age": True,
        "median_age": False,
        "ci_lower": False,
        "ci_upper": False,
        "ci_type": CIType.NONE,
        "is_converged": False,
    },
    "lsd2": {
        "mean_age": True,
        "median_age": False,
        "ci_lower": True,
        "ci_upper": True,
        "ci_type": CIType.CI95,
        "is_converged": False,
    },
    "wlogdate": {
        "mean_age": True,
        "median_age": False,
        "ci_lower": False,  # 仅当 num_replicates > 1 时提供
        "ci_upper": False,
        "ci_type": CIType.NONE,
        "is_converged": False,
    },
    "mdcat": {
        "mean_age": True,
        "median_age": False,
        "ci_lower": True,  # 通过 bootstrap CI 计算
        "ci_upper": True,
        "ci_type": CIType.CI95,
        "is_converged": False,
    },
}
