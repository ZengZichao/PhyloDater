"""
LSD2Method - LSD2 (IQ-TREE2) 定年适配器

实现基于 IQ-TREE2 的 LSD2 定年功能
- 使用 -te 参数而非 -t
- 解析 .timetree.nex 获取节点日期
- 约束语义一律经 models/constraints.to_software_format("lsd2") 集中渲染：
  点校准写 ``-T``、硬区间写 ``-老界:-新界``、单边校准写 ``-T:NA``（上界，
  年龄 ≤ T）/ ``NA:-T``（下界，年龄 ≥ T）。该单边语义经 IQ-TREE 2.3.6 /
  3.1.3 实测回读校准节点年龄确认（``-T:NA`` 解码为 age ≤ T、``NA:-T`` 解码为
  age ≥ T），凡 LSD2 无法精确表达的约束（软界/概率先验）都发出
  SemanticDegradationWarning。
- 叶节点采样日期按输入树保留，仅在确实没有日期时按"现在"（0）处理并告警
  （异时采样 / 尾型定年是 LSD2 的主要用途之一）。
"""

import re
import shlex
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, TextIO, Tuple, Union

from ..core import DatingMethod, DatingMethodRegistry
from ..core.exceptions import (
    CalibrationError,
    ConfigurationError,
    ExecutionError,
    ResultParsingError,
    SemanticDegradationWarning,
)
from ..infrastructure import ProcessRunner, get_logger
from ..infrastructure.configuration import (
    CommonConfig,
    LSD2Config,
    SoftwarePaths,
    ToolConfig,
)
from ..infrastructure.safe_io import safe_writer
from ..models import (
    AgeConstraint,
    CalibrationPoint,
    CIType,
    DatingResult,
    FixedAgeConstraint,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    NodeAgeEstimate,
    PhylogeneticTree,
    SkewNormalConstraint,
    SkewTConstraint,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
    strip_leading_newick_comments,
)

# LSD2 date 文件的取值轴：本适配器把校准年龄写成"相对现在的 Ma"（负值 = 过去），
# 因此叶节点日期也必须在同一轴上。日历式日期（2020-03-01 / 2020 / 2005.43）与之
# 不可混用，遇到即显式报错，而不是悄悄改单位。
_MA_INT_PATTERN = re.compile(r"^[+-]?\d+$")
_MA_DECIMAL_PATTERN = re.compile(r"^[+-]?\d+\.\d+$")
_CALENDAR_DATE_PATTERN = re.compile(r"^\d{4}(?:-\d{1,2}(?:-\d{1,2})?)?$")
# Newick 叶节点注解里的采样日期：``A[&&DATE=2020-03-01]`` / ``A[&date=-0.045]``
_TIP_COMMENT_PATTERN = re.compile(r"([^,()\[\]:;\s]+)\s*(\[[^\]]*\])")
_DATE_COMMENT_PATTERN = re.compile(
    r"\[(?:&&|&)?\s*(?:date|tip_date|sampling_date|collection_date|collect_date)"
    r"\s*=\s*([^\];]+)",
    re.IGNORECASE,
)


class LSD2Method(DatingMethod[LSD2Config]):
    """
    LSD2 (IQ-TREE2) 定年适配器

    特性：
    - 使用 -te 参数固定拓扑（而非 -t）
    - 约束一律经集中的 to_software_format("lsd2") 渲染（区间 / 单边 NA 语法），
      无法精确表达语义时发 SemanticDegradationWarning，绝不塌成等式约束
    - 叶节点采样日期按输入树保留，缺失才按 0（现在）处理并告警
    - 全 0 树防呆设计
    """

    def __init__(
        self,
        config: Union[ToolConfig, LSD2Config],
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
    ) -> None:
        # 入口保留双形态（历史上两种调用都存在），但在交给基类之前先归一：
        # 基类把 ``self.config`` 记为 ``LSD2Config``，适配器后续只读子配置。
        # ``isinstance`` 与原先的 ``hasattr(config, "lsd2")`` 对这两种入参等价：
        # ``ToolConfig`` 必有 ``lsd2``，``LSD2Config`` 必无。
        common: Optional[CommonConfig] = common_config
        if isinstance(config, ToolConfig):
            method_config = config.lsd2
            if common is None:
                common = config.common
        else:
            method_config = config
        super().__init__(
            method_config, output_dir, software_paths, common_config=common
        )

        self.logger = get_logger()
        self._input_files: Dict[str, Path] = {}
        self._calibrations: List[CalibrationPoint] = []
        self._root_age_cal: Optional[CalibrationPoint] = None
        self._warnings: List[Any] = (
            []
        )  # 收集需要写入 DatingResult/runtime_metadata 的警告
        self._execution_seconds: float = 0.0
        self._seq_type: Optional[str] = None

    @property
    def method_name(self) -> str:
        return "lsd2"

    def _get_executable(self) -> str:
        """获取 IQ-TREE2 可执行文件路径，优先使用 software_paths 自定义路径。"""
        path = self.get_software_path("iqtree_bin")
        return path if path else self.config.iqtree_bin

    @staticmethod
    def _detect_alignment_seq_type(alignment_path: Optional[Path]) -> Optional[str]:
        """嗅探比对序列类型，返回 ``"nucleotide"`` / ``"protein"`` / ``None``。

        只看残基字母表，不做生物学推断：核酸字母表是 ``ACGTU`` 加 IUPAC 歧义码与
        缺口符；蛋白判定只用那些**不是** IUPAC 核酸歧义码的高频氨基酸字母
        （``E Q I L F P X Z``；A/C/G/D/H/K/M/N/R/S/T/V/W/Y 既可能是氨基酸也可能
        是核酸码，不能用作判据）。读不出残基或两者都不像时返回 ``None``。
        """
        if not alignment_path or not Path(alignment_path).exists():
            return None
        nucleotide_letters = set("ACGTURWSYKVHDBNnN-*?uU")
        protein_only_letters = set("EQILFPXZ")
        seen = set()
        try:
            with open(alignment_path, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith((">", "#")):
                        continue
                    seen |= set(line.upper())
                    if seen & protein_only_letters or len(seen) > 30:
                        break
        except OSError:
            return None
        if not seen:
            return None
        seen = seen - {"0", "1"} - set(" \t.-*?u")
        if seen & protein_only_letters:
            return "protein"
        if seen and seen <= nucleotide_letters:
            return "nucleotide"
        return None

    def _resolve_model(self, alignment_path: Optional[Path]) -> str:
        """确定传给 IQ-TREE2 ``-m`` 的替换模型。

        用户显式配置的模型一律照用（只在明显与数据类型不符时报错，因为 IQ-TREE2
        此时只会丢一句 ``ERROR: File not found LG``，用户无从判断是模型选错了）；
        未配置时按比对类型自动选择：核酸 ``GTR+G``、蛋白 ``LG+G``。
        """
        seq_type = self._detect_alignment_seq_type(alignment_path)
        self._seq_type = seq_type

        explicit = (self.config.model or "").strip()
        if explicit:
            base = re.split(r"[+/]", explicit)[0].strip().upper()
            amino_acid_models = {
                "LG",
                "WAG",
                "JTT",
                "JTTDCMUT",
                "DAYHOFF",
                "MTREV",
                "MTMAM",
                "MTART",
                "MTZOA",
                "CPREV",
                "BLOSUM62",
                "BLOSUM80",
                "HONG",
                "PBREM",
                "MTZO",
                "LG45",
                "FLU",
                "INFLUENZA",
                "RTREVF93",
                "JTT1992",
            }
            if seq_type == "nucleotide" and base in amino_acid_models:
                raise ConfigurationError(
                    f"LSD2: model '{explicit}' is an amino-acid model but the alignment "
                    "looks nucleotide-based. IQ-TREE2 would abort with "
                    "'ERROR: File not found <MODEL>'. Use a nucleotide model "
                    "(e.g. 'GTR+G', 'HKY+G', 'TIM3+I+G') or leave \"model\" unset to "
                    "let PhyloDater pick one from the alignment."
                )
            if seq_type == "protein" and base in {
                "GTR",
                "HKY",
                "K2P",
                "K3PU",
                "JC",
                "F81",
                "F84",
                "TN",
                "TPM",
                "TVM",
                "TIM",
                "TIM3",
                "SYM",
            }:
                raise ConfigurationError(
                    f"LSD2: model '{explicit}' is a nucleotide model but the alignment "
                    "looks protein-based. Use an amino-acid model (e.g. 'LG+G', 'WAG+G') "
                    'or leave "model" unset.'
                )
            return explicit

        if seq_type == "protein":
            chosen = "LG+G"
        else:
            # None（无法判断）按核酸处理：PhyloDater 的主用例是 DNA 比对。
            chosen = "GTR+G"
        self.logger.info(
            f"LSD2: no substitution model configured; using '{chosen}' "
            f"based on detected alignment type ({seq_type or 'undetermined'})"
        )
        return chosen

    def validate_environment(self) -> bool:
        """检测 IQ-TREE2 是否可用"""
        runner = ProcessRunner()
        executable = self._get_executable()

        available = runner.check_executable(executable)

        if available:
            version = runner.get_version(executable, "--version")
            self.software_version = version or None
            self.logger.info(f"IQ-TREE2 detected: {version or 'unknown version'}")
        else:
            self.logger.warning(f"IQ-TREE2 not found: {executable}")

        return available

    def _get_root_representative_leaves(
        self, tree: PhylogeneticTree
    ) -> Tuple[str, str]:
        """从根节点的两个子树各取一片叶子。

        先用 dendropy（本项目的核心解析后端，在 ete3 因 Python>=3.13 移除 ``cgi``
        而不可导入的解释器上依然可用）定位根分裂；失败再试 ete3；两者都不行才退回
        "前两片叶"并显式告警——因为名字顺序上的前两片叶很可能同侧，会把根校准
        挂到错误的节点上。
        """
        newick = getattr(tree, "newick", None) or ""
        if newick:
            try:
                from dendropy import Tree as _DendropyTree

                parsed = _DendropyTree.get(
                    data=newick, schema="newick", preserve_underscores=True
                )
                children = list(parsed.seed_node.child_nodes())
                if len(children) < 2:
                    raise ValueError(f"root has only {len(children)} child node(s)")
                if len(children) > 2:
                    self.logger.warning(
                        f"LSD2: the input tree has {len(children)} children at the root "
                        "(an unrooted trifurcation), so which pair straddles the true "
                        "root is ambiguous; the root calibration is written on the "
                        "first two subtrees. Root the tree explicitly if this matters."
                    )
                leaves = []
                for child in children[:2]:
                    child_leaves = list(child.leaf_nodes())
                    if child_leaves and child_leaves[0].taxon is not None:
                        leaves.append(str(child_leaves[0].taxon.label))
                if len(leaves) == 2:
                    return leaves[0], leaves[1]
                raise ValueError("could not pick one leaf per root subtree")
            except Exception as e:
                self.logger.debug(f"dendropy root-split lookup failed: {e}")

        try:
            from ete3 import Tree

            et = Tree(strip_leading_newick_comments(newick), format=1)
            root = et.get_tree_root()
            children = [ch for ch in root.children if not ch.is_leaf()]
            if len(children) < 2:
                children = root.children
            leaves = []
            for ch in children[:2]:
                leaf = next(ch.iter_leaves(), None)
                if leaf is not None:
                    leaves.append(leaf.name)
            if len(leaves) == 2:
                return leaves[0], leaves[1]
        except Exception as e:
            self.logger.warning(f"Could not pick root representative leaves: {e}")

        names = tree.tip_names[:2]
        if len(names) < 2:
            raise CalibrationError(
                "LSD2: cannot determine the two root-representative leaves because "
                f"the tree has fewer than 2 tips (got {len(names)})."
            )
        self.logger.warning(
            "LSD2: could not identify the root split via ete3; falling back to the "
            f"first two tips ({names[0]!r}, {names[1]!r}) as the root-representative "
            "leaf pair. If these two tips are not on opposite sides of the root, the "
            "root calibration will constrain the WRONG node. Verify the generated "
            "constraints.date file, or install a working ete3/ete4."
        )
        return names[0], names[1]

    @staticmethod
    def _classify_date_token(token: str) -> Optional[str]:
        """判定日期写法的轴向。

        Returns:
            ``"ma"``      —— 可直接写入 date 文件的相对时间数值（``0``、``-0.045``、
                ``-100`` 这类"距今多少 Ma"的写法）
            ``"calendar"`` —— 日历写法（``2020`` / ``2020-03`` / ``2020-03-01`` /
                ``2005.43``），与本适配器输出的 Ma 轴节点约束不可混用
            ``None``      —— 不像日期（例如登录号尾段 ``GCA_0001`` 之后的整数），
                按"该叶没有提供日期"处理，避免把名字里的数字当成日期
        """
        token = (token or "").strip().strip('"').strip("'")
        if not token:
            return None
        if _CALENDAR_DATE_PATTERN.match(token):
            return "calendar"
        if _MA_DECIMAL_PATTERN.match(token):
            integer_part = token.lstrip("+-").split(".")[0]
            # 2005.43 这类十进历年属于日历轴；0.045 / -0.045 属于 Ma 轴
            return "calendar" if len(integer_part) == 4 else "ma"
        if _MA_INT_PATTERN.match(token):
            bare = token.lstrip("+-")
            if token.startswith("-") or bare == "0":
                return "ma"
            return "calendar" if len(bare) == 4 else None
        return None

    def _extract_tip_sampling_dates(
        self, tree: PhylogeneticTree
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """从输入树里读取叶节点的真实采样日期。

        支持两种上游写法（都是 IQ-TREE2/LSD2 生态里的通行写法）：

        1. 分类名尾随 ``|`` 的日期后缀（``--date TAXNAME`` 约定的分隔符）：
           ``EPI_ISL_12345|-0.045`` 或 ``hCoV/Wuhan|2019-12-31``。
        2. Newick 叶节点注解：``A[&&DATE=-0.045]`` / ``A[&date=-0.045]``。

        Returns:
            (numeric, calendar)：两个 ``{叶节点名: 日期串}`` 字典。``numeric`` 可直接
            写入 date 文件；``calendar`` 为轴向冲突项，由调用方显式报错。
        """
        tips = list(getattr(tree, "tip_names", None) or [])
        newick = getattr(tree, "newick", None) or ""

        comment_dates: Dict[str, str] = {}
        if newick:
            for match in _TIP_COMMENT_PATTERN.finditer(newick):
                label, comment = match.group(1), match.group(2)
                date_match = _DATE_COMMENT_PATTERN.search(comment)
                if not date_match:
                    continue
                value = date_match.group(1).strip()
                if value:
                    comment_dates.setdefault(label, value)

        numeric: Dict[str, str] = {}
        calendar: Dict[str, str] = {}
        for tip in tips:
            candidates: List[str] = []
            if "|" in tip:
                candidates.append(tip.rsplit("|", 1)[1])
            annotated = comment_dates.get(tip)
            if annotated is not None:
                candidates.append(annotated)
            for token in candidates:
                kind = self._classify_date_token(token)
                if kind == "ma":
                    numeric[tip] = token.strip().strip('"').strip("'")
                    break
                if kind == "calendar":
                    calendar[tip] = token.strip().strip('"').strip("'")
                    break
        return numeric, calendar

    def _warn_if_cannot_represent_exactly(
        self, constraint: AgeConstraint, cal_name: str, value_str: str
    ) -> None:
        """LSD2 date 文件只能表达"点 / 硬区间 / 硬单边界"。

        凡输入是软界或概率先验，都发一条 SemanticDegradationWarning 说明被写成了什么，
        使"不等式被改写"这件事在日志与结果对象里都可见（对齐 pathd8/wlogdate 的既有做法）。
        """
        if isinstance(constraint, (FixedAgeConstraint, UniformAgeConstraint)):
            # 点校准 -> -T，均匀区间 -> -老界:-新界，均为精确表达
            return

        if isinstance(constraint, SoftLowerBoundConstraint):
            msg = (
                f"LSD2 cannot represent a soft lower bound: calibration '{cal_name}' "
                f"({constraint.describe()}) was written as the HARD one-sided bound "
                f"'{value_str}'. Ages older than {constraint.min_age} Ma are now "
                "impossible instead of merely unlikely."
            )
        elif isinstance(constraint, MaximumAgeConstraint):
            msg = (
                f"LSD2 cannot represent a soft upper bound: calibration '{cal_name}' "
                f"({constraint.describe()}) was written as the HARD one-sided bound "
                f"'{value_str}'. Ages younger than {constraint.max_age} Ma are now "
                "impossible instead of merely unlikely."
            )
        elif isinstance(constraint, SoftBoundsConstraint):
            msg = (
                f"LSD2 cannot represent soft bounds: calibration '{cal_name}' "
                f"({constraint.describe()}) was written as the HARD interval "
                f"'{value_str}'; the probability mass outside the interval is lost."
            )
        elif isinstance(
            constraint, (GammaPriorConstraint, SkewNormalConstraint, SkewTConstraint)
        ):
            msg = (
                f"LSD2 cannot represent probabilistic priors: calibration "
                f"'{cal_name}' ({constraint.describe()}) was written as the hard "
                f"95%-interval '{value_str}' derived from the prior quantiles."
            )
        else:
            msg = (
                f"LSD2 has no exact representation for calibration '{cal_name}' "
                f"({constraint.describe()}); written as '{value_str}'."
            )

        self.logger.warning(msg, category=SemanticDegradationWarning)
        self._warnings.append(SemanticDegradationWarning(msg))

    def _constraint_value_for_lsd2(
        self, constraint: Optional[AgeConstraint], cal_name: str
    ) -> str:
        """把约束渲染为 LSD2 date 文件的取值串（集中走 to_software_format）。

        Raises:
            CalibrationError: 约束缺失或无法渲染时不再静默丢点，而是显式失败。
        """
        if constraint is None:
            raise CalibrationError(
                f"LSD2: calibration '{cal_name}' has no age_constraint and cannot be "
                "written to the date file."
            )
        try:
            value_str = constraint.to_software_format("lsd2")
        except Exception as e:
            raise CalibrationError(
                f"LSD2: cannot convert the age constraint of calibration "
                f"'{cal_name}' ({type(constraint).__name__}) to LSD2 format: {e}"
            ) from e
        if not value_str or not value_str.strip():
            raise CalibrationError(
                f"LSD2: age constraint of calibration '{cal_name}' rendered to an "
                "empty LSD2 date value."
            )
        value: str = value_str.strip()
        self._warn_if_cannot_represent_exactly(constraint, cal_name, value)
        return value

    def _write_tip_date_lines(
        self, f: TextIO, tree: PhylogeneticTree
    ) -> Tuple[int, int]:
        """写叶节点日期行：有日期用真实值，无日期按 0（现在）并告警。

        Returns:
            (dated, undated) 计数。
        """
        tips = list(tree.tip_names)
        numeric, calendar = self._extract_tip_sampling_dates(tree)

        if calendar:
            examples = ", ".join(
                f"{name}|...{value}" if "|" in name else f"{name}={value}"
                for name, value in list(calendar.items())[:5]
            )
            raise CalibrationError(
                f"LSD2: {len(calendar)} tip(s) carry CALENDAR sampling dates "
                f"(e.g. {examples}), but this adapter writes node calibrations as "
                "Ma relative to the present, so the two axes cannot share one date "
                "file. Convert the tip dates to the Ma axis (e.g. 'A|-0.045' for a "
                "sample taken 45 ka ago), or pass IQ-TREE's own date handling via "
                "lsd2.date_options and provide constraints in the same unit."
            )

        undated = [tip for tip in tips if tip not in numeric]
        for tip in tips:
            f.write(f"{tip}\t{numeric.get(tip, 0)}\n")

        if undated:
            preview = ", ".join(undated[:5])
            suffix = "..." if len(undated) > 5 else ""
            msg = (
                f"LSD2: {len(undated)}/{len(tips)} tip(s) carry no sampling date and "
                f"were written as 0 Ma (contemporaneous sampling): {preview}{suffix}. "
                "Heterochronous (tip-dating) analyses need the real sampling dates in "
                "the tree (name suffix after '|' or a [&&DATE=...] annotation); "
                "otherwise rate and age estimates assume a single sampling time."
            )
            self.logger.warning(msg, category=SemanticDegradationWarning)
            self._warnings.append(SemanticDegradationWarning(msg))
        elif numeric:
            self.logger.info(
                f"LSD2: preserved the sampling dates of all {len(numeric)} tip(s) "
                "from the input tree."
            )
        return len(numeric), len(undated)

    def prepare_inputs(
        self,
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        alignment_path: Optional[Path] = None,
    ) -> Dict:
        """
          生成 LSD2 输入文件（树 + IQ-TREE2 ``--date`` 采样时间/约束文件）。

          date 文件行格式（与初代 run_lsd2.py 及 IQ-TREE2 手册 §8.7-8.8 一致）：
            叶节点      ``tip<TAB>date``
            内部节点    ``tip1,tip2<TAB>value``，value 为 ``-T``（点）、
            ``-老界:-新界``（区间）、``-T:NA``（上界，年龄 ≤ T）/ ``NA:-T``
        （下界，年龄 ≥ T）（单边；语义已用 IQ-TREE 2.x/3.x 实测回读确认）。
        """
        self._calibrations = calibrations
        self._warnings = []

        # 检查内部节点校准数量
        internal_cals = [cal for cal in calibrations if not cal.is_root_node]
        if len(internal_cals) == 0:
            raise CalibrationError(
                "LSD2 requires at least one internal node calibration. "
                "A tree with only tip dates cannot infer absolute dates."
            )

        # 准备树文件
        tree_file = self.work_dir / "input_tree.nwk"
        tree.without_internal_labels().write(tree_file)

        # 生成约束文件（采样时间文件）
        constraints_file = self.work_dir / "constraints.date"

        self._root_age_cal = None
        root_cals = [c for c in calibrations if c.is_root_node]
        if len(root_cals) > 1:
            dropped = ", ".join(
                str(c.name) for c in root_cals[1:] if getattr(c, "name", None)
            )
            self.logger.warning(
                f"LSD2: {len(root_cals)} root calibrations given; only "
                f"'{root_cals[0].name}' is written and the rest "
                f"({dropped or 'unnamed'}) are ignored."
            )

        written = 0
        with safe_writer(constraints_file) as f:
            # 叶节点日期：有则用真实采样日期，无则按"现在"（0 Ma）并告警。
            self._write_tip_date_lines(f, tree)

            # 根节点约束：用跨根的两个代表叶表达根日期。
            if root_cals:
                root_cal = root_cals[0]
                self._root_age_cal = root_cal
                value_str = self._constraint_value_for_lsd2(
                    root_cal.age_constraint, root_cal.name
                )
                rt1, rt2 = self._get_root_representative_leaves(tree)
                f.write(f"{rt1},{rt2}\t{value_str}\n")
                written += 1

            for cal in calibrations:
                if cal.is_root_node:
                    continue
                if not cal.mrca_leaf_pair:
                    self.logger.warning(
                        f"Skipping calibration '{cal.name}' in LSD2: no mrca_leaf_pair"
                    )
                    continue

                value_str = self._constraint_value_for_lsd2(
                    cal.age_constraint, cal.name
                )
                tips_str = ",".join(cal.mrca_leaf_pair)
                f.write(f"{tips_str}\t{value_str}\n")
                written += 1

        if written == 0:
            raise CalibrationError(
                "LSD2: no valid calibrations written. "
                "Check that each calibration has an age constraint and a "
                "mrca_leaf_pair (two representative leaves per MRCA node)."
            )

        # 复制序列文件到工作目录（IQ-TREE --s 参数需要在 cwd 下找到）
        alignment_copy = None
        if alignment_path and alignment_path.exists():
            alignment_copy = self.work_dir / "alignment.fasta"
            import shutil

            shutil.copy2(alignment_path, alignment_copy)

        self._input_files = {
            "tree": tree_file,
            "constraints": constraints_file,
        }
        # 比对是可选输入：没有比对时不要往“文件字典”里塞 None。
        # 下游两个读取处（本适配器的 execute 与基类 cleanup 的归档）都按
        # 路径对待这些条目；旧写法会让 ``Optional[Path]`` 混进 Dict[str, Path]。
        resolved_alignment = alignment_copy or alignment_path
        if resolved_alignment is not None:
            self._input_files["alignment"] = resolved_alignment

        self.logger.info(f"Generated LSD2 constraint file: {constraints_file}")
        return self._input_files

    def execute(self) -> bool:
        """执行 LSD2 分析"""
        tree_file = self._input_files["tree"]
        constraints_file = self._input_files["constraints"]
        alignment_path = self._input_files.get("alignment")

        if not alignment_path or not Path(alignment_path).exists():
            raise ExecutionError("Missing alignment file for LSD2 analysis")

        # 构建 IQ-TREE2 命令 — 参数顺序与初代 run_lsd2.py 完全一致（已验证工作）
        # 注意：--date 工作流必须显式传 -m，否则 IQ-TREE2 会跑 modelFinder
        # 遍历全部蛋白模型（默认 1232 个），在较大的比对上极慢。
        model = self._resolve_model(alignment_path)
        cmd = [
            self._get_executable(),
            "-s",
            alignment_path.name,
            "-te",
            tree_file.name,
            "--date",
            constraints_file.name,
            "-m",
            model,
            "-pre",
            "lsd2_run",
        ]

        if self.config.date_options:
            cmd += shlex.split(self.config.date_options)

        self.logger.info(f"Executing: {' '.join(cmd[:])} ...")

        runner = ProcessRunner(
            cwd=self.work_dir, timeout=self.common_config.timeout or 7200
        )
        started = time.monotonic()
        result = runner.run(cmd)
        self._execution_seconds = time.monotonic() - started

        if result.returncode != 0:
            raise ExecutionError(
                f"LSD2 execution failed (returncode={result.returncode}): {result.stderr}"
            )

        # 输出文件在 work_dir/lsd2_run.timetree.nex/.wk/.ld
        self._input_files["timetree_nex"] = self.work_dir / "lsd2_run.timetree.nex"
        self._input_files["timetree_nwk"] = self.work_dir / "lsd2_run.timetree.nwk"
        self._input_files["timetree_lsd"] = self.work_dir / "lsd2_run.timetree.lsd"
        self._input_files["log_file"] = self.work_dir / "lsd2_run.log"

        self.logger.success("LSD2 execution completed")
        return True

    def parse_results(self) -> DatingResult:
        """解析 LSD2 输出

        LSD2/IQ-TREE2 输出的时间树是 chronogram：根节点年龄 = 0，向叶节点方向递增，
        年龄编码在分支长度的累加中。IQ-TREE2 同时在每个节点的 ``[&date="..."]``
        注解中直接写入该节点的绝对时间（负值表示距今，即 ``date = 根年龄 −
        distance_from_root``）。

        解析按以下优先级：
        1. 主路径：全树解析 ``[&date=...]`` 注解（每节点都有），直接得到节点时间（取
           ``−date`` 转为正值的 Ma）。
        2. 回退路径：若注解解析失败，按 MRCA 拓扑以 root-to-node 距离解码年龄
           （``_extract_chronogram_ages``）。

        LSD2 输出通常不含内部节点 label，故匹配校准 MRCA 的节点以校准点名作键，
        其余带 label 的节点用 label，未命名内部节点用 ``internal_node_<序号>``。
        """
        from dendropy import Tree as _DendropyTree

        dated_tree_newick = ""
        source_path = None

        # 优先从 NEXUS 解析，回退到 NWK
        nex_path = self._input_files.get("timetree_nex")
        if nex_path and Path(nex_path).exists():
            source_path = nex_path
        if not source_path:
            nwk_path = self._input_files.get("timetree_nwk")
            if nwk_path and Path(nwk_path).exists():
                source_path = nwk_path

        if not source_path:
            raise ResultParsingError("LSD2 timetree output not found")

        with open(source_path) as f:
            dated_tree_newick = f.read().strip()

        # ------------------------------------------------------------------
        # 主路径：全树解析 [&date=...] 注解。IQ-TREE2 默认 .timetree.nex 的每
        # 个节点都携带 ``[&date="..."]`` 注解，即该节点的绝对时间（负 = 距今）。
        #
        # 注意：node_ages 必须在 schema 循环**内部**重置，否则第一次尝试（nexus）
        # 中途抛异常时已填入的残片会留在字典里，而第二次尝试（newick）的同名键
        # 因 "key not in node_ages" 保护而无法覆盖上一次失败尝试的结果。
        # ------------------------------------------------------------------
        node_ages: Dict[str, NodeAgeEstimate] = {}
        parse_errors = []
        negative_date_annotations = 0
        for schema in ("nexus", "newick"):
            node_ages = {}
            try:
                chronogram = _DendropyTree.get(
                    data=dated_tree_newick, schema=schema, preserve_underscores=True
                )
                try:
                    if chronogram.is_rooted is None or chronogram.is_rooted is False:
                        chronogram.is_rooted = True
                except Exception:
                    pass

                # 用全名在树中定位各校准的 MRCA 节点（LSD2 输出保留全名）。
                # 根节点校准没有 mrca_leaf_pair，直接映射到 chronogram 的 seed_node。
                # 同时持有节点对象引用（而不只是 id），并以身份比对取回，避免
                # "只存 id、不存引用"时对象被回收后 id 被复用而静默贴错名字。
                mrca_cal_by_node: Dict[int, Tuple[Any, str]] = {}
                for cal in self._calibrations:
                    if not cal.name:
                        continue
                    if cal.is_root_node:
                        mrca_cal_by_node[id(chronogram.seed_node)] = (
                            chronogram.seed_node,
                            cal.name,
                        )
                        continue
                    if not cal.mrca_leaf_pair:
                        continue
                    taxa = []
                    for tip_name in cal.mrca_leaf_pair:
                        taxon = chronogram.taxon_namespace.get_taxon(tip_name)
                        if taxon is None:
                            break
                        taxa.append(taxon)
                    if len(taxa) == len(cal.mrca_leaf_pair):
                        try:
                            mrca = chronogram.mrca(taxa=taxa)
                            if mrca is not None:
                                mrca_cal_by_node[id(mrca)] = (mrca, cal.name)
                        except Exception:
                            pass

                total_nodes = 0
                internal_idx = 1
                for node in chronogram.preorder_node_iter():
                    total_nodes += 1
                    date_val = None
                    for ann in node.annotations:
                        ann_name = str(ann.name) if ann.name is not None else ""
                        if ann_name == "date" or ann_name.endswith(":date"):
                            try:
                                date_val = float(ann.value)
                            except Exception:
                                date_val = None
                            break
                    if date_val is None:
                        continue
                    # date 为负（距今），取负得到正的节点年龄（Ma）
                    age = -date_val
                    if age < 0:
                        negative_date_annotations += 1

                    # 过滤叶节点：PhyloDater 的节点年龄表只保留内部/校准节点
                    if node.is_leaf():
                        continue

                    entry = mrca_cal_by_node.get(id(node))
                    cal_name = (
                        entry[1] if entry is not None and entry[0] is node else ""
                    )
                    node_label = node.label or ""

                    # 优先使用校准名或节点 label；未命名内部节点使用稳定序号
                    if cal_name:
                        key = cal_name
                    elif node_label:
                        key = node_label
                    else:
                        key = f"internal_node_{internal_idx}"
                        internal_idx += 1
                    if key and key not in node_ages:
                        node_ages[key] = NodeAgeEstimate(
                            mean_age=age, ci_type=CIType.NONE
                        )

                if node_ages:
                    self.logger.info(
                        f"LSD2: annotated {len(node_ages)}/{total_nodes} node ages from [&date=...]"
                    )
                    break
                else:
                    parse_errors.append(f"{schema}: no [&date=...] ages extracted")
            except Exception as e:
                parse_errors.append(f"{schema}: {e}")
                continue

        if negative_date_annotations:
            msg = (
                f"LSD2: {negative_date_annotations} node date annotation(s) were "
                "positive, so 'age = -date' produced NEGATIVE node ages. The output "
                "dates are not on the 'negative = past' axis this adapter writes "
                "(Ma relative to the present) — check the date unit of the input "
                "calibrations/tip dates. Values are kept as computed (not abs())."
            )
            self.logger.warning(msg, category=SemanticDegradationWarning)
            self._warnings.append(SemanticDegradationWarning(msg))

        # ------------------------------------------------------------------
        # 回退路径：若注解解析未能提取到任何年龄，按 MRCA 拓扑以 root-to-node
        # 距离解码 chronogram 年龄（_extract_chronogram_ages）。
        # ------------------------------------------------------------------
        if not node_ages:
            for schema in ("nexus", "newick"):
                try:
                    chronogram = _DendropyTree.get(
                        data=dated_tree_newick,
                        schema=schema,
                        preserve_underscores=True,
                    )
                    node_ages = self._extract_chronogram_ages(
                        chronogram, self._calibrations, unit_factor=1.0
                    )
                    if node_ages:
                        break
                except Exception as e:
                    parse_errors.append(f"fallback {schema}: {e}")
                    continue

        if not node_ages:
            err_detail = "; ".join(parse_errors) if parse_errors else "unknown"
            self.logger.warning(f"Could not decode LSD2 chronogram ages ({err_detail})")

        if not node_ages:
            self.logger.warning(
                "LSD2: parsed a timetree but extracted no node ages. "
                "The chronogram may have tip-name mismatches with calibrations."
            )

        result = DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=dated_tree_newick,
            node_ages=node_ages,
            raw_output_path=self.work_dir,
            execution_seconds=self._execution_seconds,
            metadata={"execution_time": self._execution_seconds},
        )

        for warning in self._warnings:
            result.add_warning(warning)

        self.logger.success(f"Parsed LSD2 results: {len(node_ages)} node ages")
        return result


# 注册适配器
DatingMethodRegistry.register("lsd2", LSD2Method)
