"""
AutoCalibrator - 自动校准点添加服务

通过指定类群名称和校准时间，自动计算MRCA并添加校准信息。
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

from ..infrastructure.logging import get_logger
from ..infrastructure.safe_io import safe_writer
from ..models import CalibrationPoint, FixedAgeConstraint, PhylogeneticTree
from ..models.constraints import UniformAgeConstraint
from .calibration_resolver import CalibrationResolver
from .taxonomy_parser import TaxonomyParser


@dataclass
class AutoCalibration:
    """自动校准点配置"""

    taxon: str  # 类群名称
    age: float  # 校准时间 (Ma)
    age_max: Optional[float] = None  # 最大年龄 (用于范围校准)
    is_range: bool = False  # 是否为范围校准


@dataclass
class CalibratedTree:
    """带校准信息的树"""

    tree: PhylogeneticTree
    calibrations: List[CalibrationPoint]
    output_path: Optional[Path] = None


class AutoCalibrator:
    """
    自动校准点添加服务

    通过指定类群名称和校准时间，自动计算MRCA并添加校准信息。
    """

    def __init__(
        self, tree: PhylogeneticTree, taxonomy_parser: Optional[TaxonomyParser] = None
    ) -> None:
        """
        初始化自动校准器

        Args:
            tree: 系统发育树
            taxonomy_parser: 分类学解析器（可选）
        """
        self.tree = tree
        self.taxonomy_parser = taxonomy_parser or TaxonomyParser()
        self.resolver = CalibrationResolver(tree, self.taxonomy_parser)
        self.logger = get_logger()

    @staticmethod
    def describe_mrca_location(
        mrca_leaf_pair: Optional[Tuple[str, str]], is_root_node: bool = False
    ) -> str:
        """把校准点的节点位置渲染成人类可读串（B-20 的唯一出口）。

        ``mrca_leaf_pair`` 在 :class:`CalibrationPoint` 里是 ``Optional[tuple]``，
        而且 ``None`` **不是异常状态**：根校准（``is_root_node=True``）由
        ``RootAge`` 处理、本就没有叶对。此前三处汇报点各自无条件写
        ``mrca_leaf_pair[0]``，于是 ``--auto-calibrate "LUCA:3800"`` 这类合法输入
        在"解析成功"的那一行上崩成 ``TypeError: 'NoneType' object is not
        subscriptable``——崩溃点与真实成因隔了三层。
        """
        if mrca_leaf_pair:
            pair = tuple(mrca_leaf_pair)
            if len(pair) >= 2:
                return f"({pair[0]}, {pair[1]})"
            return f"({pair[0]})"
        return "(root)" if is_root_node else "(unresolved)"

    def parse_calibration_string(self, cal_str: str) -> AutoCalibration:
        """
        解析校准字符串

        格式:
        - "TAXON:AGE" - 固定年龄校准
        - "TAXON:AGE_MIN-AGE_MAX" - 范围校准

        Args:
            cal_str: 校准字符串

        Returns:
            AutoCalibration 实例

        Raises:
            ValueError: 格式错误
        """
        # 检查格式
        if ":" not in cal_str:
            raise ValueError(
                f"校准字符串格式错误: '{cal_str}'\n"
                f"正确格式: 类群名称:校准时间(Ma)\n"
                f"例如: Cyanobacteriota:2500 或 Bacteria:3500-4000"
            )

        parts = cal_str.split(":", 1)
        taxon = parts[0].strip()
        age_str = parts[1].strip()

        if not taxon:
            raise ValueError(f"类群名称不能为空: '{cal_str}'")

        # 解析年龄
        if "-" in age_str:
            # 范围校准
            age_parts = age_str.split("-", 1)
            try:
                age_min = float(age_parts[0].strip())
                age_max = float(age_parts[1].strip())
            except ValueError:
                raise ValueError(
                    f"年龄格式错误: '{age_str}'\n"
                    f"应为数字或范围 (如 2500 或 3500-4000)"
                )
            # 审阅项 C-28: 用户写 "TAXON:2500-2500" 想要的就是固定点校准。
            # 旧实现把它交给 UniformAgeConstraint，被 min_age >= max_age 拒绝后
            # 报出 "min_age (2500.0) must be < max_age (2500.0)" 这种讲内部变量名、
            # 不含任何可行替代写法的错误。这里直接转成 FixedAgeConstraint。
            if age_min == age_max:
                self.logger.info(
                    f"'{cal_str}' 的上下界相等 ({age_min} Ma)，按固定年龄校准处理 "
                    f"(等价于 '{taxon}:{age_min:g}')"
                )
                return AutoCalibration(taxon=taxon, age=age_min)
            return AutoCalibration(
                taxon=taxon, age=age_min, age_max=age_max, is_range=True
            )
        else:
            # 固定年龄校准
            try:
                age = float(age_str)
            except ValueError:
                raise ValueError(f"年龄格式错误: '{age_str}'\n" f"应为数字 (如 2500)")
            return AutoCalibration(taxon=taxon, age=age)

    def resolve_calibration(self, auto_cal: AutoCalibration) -> CalibrationPoint:
        """
        解析自动校准点

        Args:
            auto_cal: 自动校准点配置

        Returns:
            CalibrationPoint 实例

        Raises:
            ValueError: 解析失败
        """
        self.logger.info(f"解析校准点: {auto_cal.taxon}")

        # 使用 CalibrationResolver 解析类群
        try:
            result = self.resolver.resolve(auto_cal.taxon)
        except Exception as e:
            raise ValueError(f"无法解析类群 '{auto_cal.taxon}': {e}")

        # 创建年龄约束
        age_constraint: Union[UniformAgeConstraint, FixedAgeConstraint]
        if auto_cal.is_range:
            if auto_cal.age_max is None:
                # 标记了区间却没有上界：构造 uniform 约束会把 None 当年龄使，
                # 比当场报错恶劣得多（parse_auto_calibration 的正常产出总是同时
                # 给出两个值，走到这里只能是调用方手拼的）。
                raise ValueError(
                    f"范围校准 '{auto_cal.taxon}' 标了 is_range 但缺 age_max，"
                    "uniform 约束必须同时有上下界"
                )
            age_constraint = UniformAgeConstraint(
                min_age=auto_cal.age, max_age=auto_cal.age_max
            )
        else:
            age_constraint = FixedAgeConstraint(fixed_age=auto_cal.age)

        # 创建校准点
        calibration = CalibrationPoint(
            name=auto_cal.taxon,
            age_constraint=age_constraint,
            resolved_taxa=result.resolved_taxa,
            mrca_leaf_pair=result.mrca_leaf_pair,
            is_root_node=result.is_root_node,
        )

        if result.mrca_leaf_pair is None and not result.is_root_node:
            # 不是根校准却没有叶节点对：下游适配器（PATHd8/treePL/r8s）只能把这个
            # 校准点静默丢弃（见审阅项 B-6/C-20），这里必须提前说清楚。
            self.logger.warning(
                f"校准点 '{auto_cal.taxon}' 未解析到 MRCA 叶节点对且不是根校准，"
                f"依赖叶节点对的后端（PATHd8/treePL/r8s）将跳过该校准点"
            )

        self.logger.success(
            f"校准点 '{auto_cal.taxon}' 解析成功: "
            f"MRCA = {self.describe_mrca_location(result.mrca_leaf_pair, result.is_root_node)}, "
            f"年龄 = {auto_cal.age} Ma"
        )

        return calibration

    def auto_calibrate(self, cal_strings: List[str]) -> CalibratedTree:
        """
        自动添加校准点

        Args:
            cal_strings: 校准字符串列表

        Returns:
            CalibratedTree 实例

        Raises:
            ValueError: 解析失败
        """
        self.logger.section("自动校准点添加")

        calibrations = []
        for cal_str in cal_strings:
            auto_cal = self.parse_calibration_string(cal_str)
            calibration = self.resolve_calibration(auto_cal)
            calibrations.append(calibration)

        # 输出带校准信息的树
        self.logger.info(f"成功添加 {len(calibrations)} 个校准点")

        return CalibratedTree(tree=self.tree, calibrations=calibrations)

    def write_calibrated_tree(
        self, calibrated_tree: CalibratedTree, output_path: Path
    ) -> None:
        """
        输出带校准信息的树文件（人类可读的注释 + 原始 Newick）。

        注意：本方法写入的是"带人类可读校准注释的 Newick"——约束信息以注释行
        形式输出，并**不**把约束写回树节点（树节点标注由 mcmctree/treePL 等
        适配器的 prepare_inputs 负责）。因此这不是一个可供下游直接使用的
        "带约束校准树"，文件名中 "calibrated tree" 仅表示"附带校准信息的树文件"。
        如需下游可用的约束树，请走正常解析流程（CalibrationResolver + 适配器）。

        Args:
            calibrated_tree: 带校准信息的树
            output_path: 输出文件路径
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # 获取树的 Newick 字符串
        newick = calibrated_tree.tree.newick

        # 创建校准信息注释
        cal_info = []
        for cal in calibrated_tree.calibrations:
            mrca_pair = cal.mrca_leaf_pair
            age_constraint = cal.age_constraint

            # 复用 constraints 的集中式 describe()，覆盖全部 8 类约束
            # （此前仅 Fixed/Uniform 被识别，其余显示 "unknown"）。
            if age_constraint is not None:
                age_str = age_constraint.describe()
            else:
                age_str = "unknown"

            cal_info.append(
                f"# {cal.name}: MRCA "
                f"{self.describe_mrca_location(mrca_pair, cal.is_root_node)} = {age_str}"
            )

        # 写入文件
        with safe_writer(output_path, encoding="utf-8") as f:
            # 写入校准信息注释
            f.write("# PhyloDater Auto-Calibration\n")
            f.write(f"# Total calibrations: {len(calibrated_tree.calibrations)}\n")
            for info in cal_info:
                f.write(f"{info}\n")
            f.write("\n")
            # 写入树
            f.write(newick + "\n")

        self.logger.success(f"校准树已保存到: {output_path}")

    def get_calibration_summary(self, calibrations: List[CalibrationPoint]) -> str:
        """
        获取校准点摘要

        Args:
            calibrations: 校准点列表

        Returns:
            摘要字符串
        """
        lines = ["校准点摘要:"]

        for i, cal in enumerate(calibrations, 1):
            mrca_pair = cal.mrca_leaf_pair
            age_constraint = cal.age_constraint

            # 复用 constraints 的集中式 describe()，覆盖全部 8 类约束
            if age_constraint is None:
                age_str = "unknown"
            else:
                age_str = age_constraint.describe()

            lines.append(f"  {i}. {cal.name}")
            lines.append(
                f"     MRCA: "
                f"{self.describe_mrca_location(mrca_pair, cal.is_root_node)}"
            )
            lines.append(f"     年龄: {age_str}")
            lines.append(f"     叶节点数: {len(cal.resolved_taxa or [])}")

        return "\n".join(lines)


def parse_auto_calibrate_args(args: List[str]) -> List[str]:
    """
    解析命令行中的自动校准参数

    支持格式:
    - "TAXON:AGE" - 单个校准点
    - "TAXON:AGE_MIN-AGE_MAX" - 范围校准
    - 多个参数: "TAXON1:AGE1" "TAXON2:AGE2"
    - 同一参数内的多个校准点: "TAXON1:AGE1 TAXON2:AGE2"

    Args:
        args: 命令行参数列表

    Returns:
        校准字符串列表

    Raises:
        ValueError: 参数里出现不含 ``TAXON:AGE`` 形状的片段

    Note:
        审阅项 C-27：旧实现按空格二次切分后，**不含 ``:`` 的片段被静默丢弃**
        （实测 ``"Bacteria:3500 typo"`` 只剩 ``Bacteria:3500``，用户以为传了两个
        约束、实际只有一个生效）。同一条规则在 else 分支是会 raise 的，两条路径
        一严一松。现统一为"任何被丢弃的片段一律 raise"。
        注意按空格切分意味着**含空格的类群名**无法用这条语法表达；这类名字请写进
        校准 YAML（``--calibrations``），或用引号包住整个参数后以分号分隔。
    """
    if not args:
        return []

    cal_strings = []
    for arg in args:
        if not isinstance(arg, str):
            raise ValueError(
                f"--auto-calibrate 参数必须是字符串，收到 {type(arg).__name__}: {arg!r}"
            )
        # 支持空格分隔的多个校准点
        if " " in arg:
            parts = arg.split()
            for part in parts:
                if ":" in part:
                    cal_strings.append(part)
                else:
                    raise ValueError(
                        f"--auto-calibrate 参数 '{arg}' 中的片段 '{part}' 不是合法的"
                        f"校准描述（应为 类群名称:年龄(Ma)，如 Cyanobacteriota:2500）。\n"
                        f"该片段既没有被解析、也不会被应用，因此这里直接报错而不是"
                        f"静默丢弃。若类群名本身含空格，请改用校准文件 --calibrations。"
                    )
        else:
            cal_strings.append(arg)

    return cal_strings
