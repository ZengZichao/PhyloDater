"""
DeepValidator - 输入格式深层验证与错误定位

提供树文件和序列文件的深层验证功能：
- 树文件：括号平衡、分支长度非负、节点名检查、格式支持、**超度量性/时间树判定**
- 序列文件：字母表检查、ID唯一性、长度一致性、格式识别
"""

import math
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple, Union

from ..infrastructure.logging import get_logger

if TYPE_CHECKING:  # 只供标注，不引入运行时依赖
    from ..models import PhylogeneticTree

# Newick 合法节点名规则：
#  - 引号包裹（可含空格/括号/逗号等特殊字符）："..."（内部 "" 视为转义引号）
#  - 裸名：不含结构字符 (空白 , : ( ) ; [ ])；允许以数字/字母/下划线开头
# 允许数字开头的标签（Newick 规范未禁止），并支持引号包裹的标签。
_NODE_NAME_CHARCLASS = r"[^\\s,:();\[\]]"
_NODE_NAME_RE = re.compile(r'[,(]\s*("(?:[^"]|"")*"|' + _NODE_NAME_CHARCLASS + r"+)")
_INTERNAL_NODE_RE = re.compile(
    r'\)\s*("(?:[^"]|"")*"|' + _NODE_NAME_CHARCLASS + r"+)\s*(?=:)"
)
_SELF_LOOP_RE = re.compile(
    r'("(?:[^"]|"")*"|' + _NODE_NAME_CHARCLASS + r"+):0(?:\.0+)?,\s*\1(?=[,):])"
)


def _strip_node_name(name: str) -> str:
    """去除节点名可能的引号并还原内部转义（"" -> "）。"""
    if len(name) >= 2 and name.startswith('"') and name.endswith('"'):
        return name[1:-1].replace('""', '"')
    return name


class SequenceAlphabet(Enum):
    """序列字母表类型"""

    DNA = "dna"
    RNA = "rna"
    PROTEIN = "protein"
    UNKNOWN = "unknown"


class SequenceFormat(Enum):
    """序列文件格式"""

    FASTA = "fasta"
    FASTQ = "fastq"
    UNKNOWN = "unknown"


@dataclass
class TreeValidationDetail:
    """树文件验证详情"""

    is_valid: bool
    format_detected: Optional[str] = None
    tree_count: int = 0
    tree_summaries: List[Dict] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class SequenceValidationDetail:
    """序列文件验证详情"""

    is_valid: bool
    format_detected: Optional[str] = None
    alphabet_detected: Optional[str] = None
    num_sequences: int = 0
    sequence_length: Optional[int] = None
    is_aligned: bool = True
    duplicate_ids: List[str] = field(default_factory=list)
    invalid_chars: List[Dict] = field(default_factory=list)  # [{line, char, position}]
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# 超度量性 / "是不是时间树"（审阅项 B-18）
# --------------------------------------------------------------------------- #

#: root-to-tip 距离相对极差的默认容差（1e-6 相对值，与审阅报告的建议一致）
ULTRAMETRIC_DEFAULT_TOLERANCE = 1e-6

#: 这些方法的设计前提**允许**非超度量输入（速率平滑/相对定年），
#: 因此只能给 info，不能拦用户。
ULTRAMETRIC_TOLERANT_METHODS = frozenset(
    {"treepl", "pathd8", "r8s", "pyr8s", "mdcat", "md-cat", "wlogdate"}
)

#: 这些方法（或本项目的相应通路）把"树已是时间树"当作前提：
#: LSD2 最小二乘定年、IQ-TREE2 尖端定年、MCMCTree 全局钟 clock=1，
#: 以及 ``viz`` 层 ``_compute_node_age`` 的 ``root_age - dist(root, node)`` 公式。
TIME_TREE_ASSUMING_METHODS = frozenset(
    {"lsd2", "iqtree", "iqtree2", "mcmctree", "viz", "visualization"}
)


@dataclass
class UltrametricDetail:
    """一棵树的超度量性（时间树）判定结果"""

    assessable: bool
    is_time_tree: bool
    n_tips: int = 0
    min_height: Optional[float] = None
    max_height: Optional[float] = None
    relative_spread: Optional[float] = None
    tolerance: float = ULTRAMETRIC_DEFAULT_TOLERANCE
    worst_tips: List[Tuple[str, float]] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def summarize(self) -> str:
        """一句话描述（面向日志/终端）。"""
        if not self.assessable:
            return "无法判定输入树是否为时间树（超度量性）：" + (
                "; ".join(self.notes) or "树不可解析或缺少分支长度"
            )
        if self.is_time_tree:
            return (
                f"输入树是时间树（{self.n_tips} 个叶节点的 root-to-tip 距离相对极差 "
                f"{self.relative_spread:.3g} <= 容差 {self.tolerance:.3g}）"
            )
        offenders = ", ".join(f"{name}({dev:.3g})" for name, dev in self.worst_tips[:5])
        return (
            f"输入树**不是**时间树（非超度量）：{self.n_tips} 个叶节点的 root-to-tip "
            f"距离从 {self.min_height:.6g} 到 {self.max_height:.6g}，相对极差 "
            f"{self.relative_spread:.3g} > 容差 {self.tolerance:.3g}；"
            f"偏离最大者: {offenders}"
        )


# 定义合法的序列字符集
DNA_CHARS = set("ATCGUNatcgun-?.")  # 包含简并碱基和gap
RNA_CHARS = set("ACGUNacgun-?.")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx-?.*")


class DeepValidator:
    """
    深度输入验证器

    提供树文件和序列文件的深层验证与错误定位功能。
    """

    MAX_INPUT_SIZE = 500 * 1024 * 1024

    def __init__(self, ignore_malformed: bool = False) -> None:
        self.logger = get_logger()
        self.ignore_malformed = ignore_malformed

    # ==================== 树文件验证 ====================

    def validate_tree_deep(
        self, file_path: Path, ignore_malformed: Optional[bool] = None
    ) -> TreeValidationDetail:
        """
        深度验证树文件

        检查项目：
        1. 括号平衡
        2. 引号转义
        3. 分支长度非负
        4. 空节点名检测
        5. 重复节点名检测
        6. 多棵树检测和概要

        Args:
            file_path: 树文件路径
            ignore_malformed: 是否将恶意字符/格式错误降级为警告（None 表示使用实例默认值）

        Returns:
            TreeValidationDetail: 验证详情
        """
        effective_ignore = (
            ignore_malformed if ignore_malformed is not None else self.ignore_malformed
        )
        if not file_path.exists():
            return TreeValidationDetail(
                is_valid=False, errors=[f"树文件不存在: {file_path}"]
            )

        if file_path.stat().st_size == 0:
            return TreeValidationDetail(
                is_valid=False, errors=[f"树文件为空: {file_path}"]
            )

        if file_path.stat().st_size > self.MAX_INPUT_SIZE:
            from ..core.exceptions import InputSizeLimitError

            raise InputSizeLimitError(
                f"树文件过大 ({file_path.stat().st_size / (1024*1024):.0f} MB), "
                f"超过限制 ({self.MAX_INPUT_SIZE // (1024*1024)} MB)"
            )

        # 检测格式
        fmt = self._detect_tree_format(file_path)

        # 读取文件内容
        try:
            content = self._read_file_safe(file_path)
        except Exception as e:
            return TreeValidationDetail(is_valid=False, errors=[f"读取树文件失败: {e}"])

        # 根据格式进行验证。effective_ignore 必须一路传到恶意字符检查：
        # 旧实现在这里算出它却从不使用，于是按调用点传入的 ignore_malformed
        # 被静默丢弃，行为始终退回实例默认值（ruff F841 指向的正是这个缺陷）。
        if fmt == "newick" or fmt == "nhx":
            return self._validate_newick_deep(
                content, file_path, fmt, ignore_malformed=effective_ignore
            )
        elif fmt == "nexus":
            return self._validate_nexus_deep(
                content, file_path, ignore_malformed=effective_ignore
            )
        else:
            from ..core.exceptions import PhyloFormatError

            raise PhyloFormatError(
                f"不支持的树文件格式: {file_path}",
                suggestion="支持的格式: Newick (.nwk, .tree), Nexus (.nex, .nexus)",
            )

    def _detect_tree_format(self, file_path: Path) -> str:
        """检测树文件格式"""
        ext = file_path.suffix.lower()

        # 根据扩展名判断
        if ext in [".nwk", ".newick", ".tree", ".treefile", ".tre"]:
            return "newick"
        elif ext in [".nhx"]:
            return "nhx"
        elif ext in [".nxs", ".nexus", ".nex"]:
            return "nexus"

        # 根据内容判断
        try:
            with open(file_path, "r", encoding="utf-8", errors="strict") as f:
                first_line = f.readline().strip()

            if first_line.startswith("#NEXUS") or first_line.startswith("#nexus"):
                return "nexus"
            elif first_line.startswith("(") or ";" in first_line:
                return "newick"
        except UnicodeDecodeError as e:
            self.logger.warning(
                f"文件编码问题，无法以 UTF-8 严格模式读取: {file_path} - {e}"
            )
        except OSError:
            pass

        return "unknown"

    def _read_file_safe(self, file_path: Path) -> str:
        """安全读取文件内容"""
        encodings = ["utf-8-sig", "utf-8", "latin-1", "ascii"]

        for encoding in encodings:
            try:
                with open(file_path, "r", encoding=encoding) as f:
                    return f.read()
            except UnicodeDecodeError:
                continue

        raise ValueError(f"无法以任何编码读取文件: {file_path}")

    def _validate_newick_deep(
        self,
        content: str,
        file_path: Path,
        fmt: str,
        ignore_malformed: Optional[bool] = None,
    ) -> TreeValidationDetail:
        """深度验证 Newick/NHX 格式树文件"""
        errors = []
        warnings = []
        tree_summaries = []

        # 移除注释（NHX格式的 [&...] 注释）
        clean_content = re.sub(r"\[&&[^\]]*\]", "", content)

        # 按分号分割多棵树
        tree_strings = [t.strip() for t in clean_content.split(";") if t.strip()]

        if not tree_strings:
            return TreeValidationDetail(
                is_valid=False, errors=[f"树文件中没有找到有效的树: {file_path}"]
            )

        tree_count = len(tree_strings)
        if tree_count > 1:
            warnings.append(f"树文件包含 {tree_count} 棵树")

        for i, tree_str in enumerate(tree_strings):
            tree_num = i + 1
            prefix = f"树 {tree_num}" if tree_count > 1 else "树"

            # 1. 括号平衡检查
            bracket_errors = self._check_bracket_balance(tree_str, prefix)
            errors.extend(bracket_errors)

            # 2. 分支长度非负检查 (CRITICAL 级别)
            length_critical, length_warn = self._check_branch_lengths(tree_str, prefix)
            if length_critical:
                from ..core.exceptions import NegativeBranchLengthError

                raise NegativeBranchLengthError(
                    f"{prefix}包含 {len(length_critical)} 个负分支长度",
                    suggestion="请检查树文件中分支长度是否正确，或使用 --skip-length-check 跳过此检查。",
                )
            warnings.extend(length_warn)

            # 3. 节点名检查
            node_errors, node_warnings = self._check_node_names(
                tree_str, prefix, ignore_malformed=ignore_malformed
            )
            errors.extend(node_errors)
            warnings.extend(node_warnings)

            # 4. 多根节点和自循环边检测 (CRITICAL)
            critical_errors = self._check_multi_root_and_self_loop(tree_str, prefix)
            if critical_errors:
                from ..core.exceptions import TreeValidationError

                raise TreeValidationError(
                    "；".join(critical_errors),
                    suggestion="请检查树文件是否包含多棵独立树或自循环边（节点连接到自身）。",
                )

            # 5. 统计节点信息
            summary = self._get_tree_summary(tree_str, tree_num)

            # 6. 超度量性/时间树判定（审阅项 B-18）：只报告、不拦下——
            #    treePL/PATHd8/r8s 本来就能吃非超度量树，而 LSD2/尖端定年/可视化
            #    把它当前提。真正的按方法分诊在 TreeValidator.full_validation 里。
            ultra_warnings = self._append_ultrametric_summary(summary, tree_str, prefix)
            warnings.extend(ultra_warnings)
            tree_summaries.append(summary)

            if tree_count > 1:
                warnings.append(
                    f"  树 {tree_num}: {summary['terminals']} 个末端节点, "
                    f"{summary['internal']} 个内部节点"
                )

        return TreeValidationDetail(
            is_valid=len(errors) == 0,
            format_detected=fmt,
            tree_count=tree_count,
            tree_summaries=tree_summaries,
            errors=errors,
            warnings=warnings,
        )

    def _validate_nexus_deep(
        self,
        content: str,
        file_path: Path,
        ignore_malformed: Optional[bool] = None,
    ) -> TreeValidationDetail:
        """深度验证 Nexus 格式树文件"""
        errors = []
        warnings = []
        tree_summaries = []

        # 检查是否以 #NEXUS 开头
        if not content.strip().upper().startswith("#NEXUS"):
            return TreeValidationDetail(
                is_valid=False,
                format_detected="nexus",
                errors=[f"Nexus 文件必须以 '#NEXUS' 开头: {file_path}"],
            )

        # 提取 TREE 块
        tree_block_pattern = re.compile(
            r"TREE\s+(?:\*\s+)?(\w+)\s*=\s*(\[&[^\]]*\]\s*)?\(?[^;]+;",
            re.IGNORECASE | re.DOTALL,
        )
        tree_matches = tree_block_pattern.findall(content)

        if not tree_matches:
            # 尝试更宽松的匹配
            tree_pattern = re.compile(r"TREE\s+.*?;", re.IGNORECASE | re.DOTALL)
            tree_matches_raw = tree_pattern.findall(content)

            if not tree_matches_raw:
                return TreeValidationDetail(
                    is_valid=False,
                    format_detected="nexus",
                    errors=[f"Nexus 文件中没有找到 TREE 块: {file_path}"],
                )

            tree_count = len(tree_matches_raw)
        else:
            tree_count = len(tree_matches)

        if tree_count > 1:
            warnings.append(f"Nexus 文件包含 {tree_count} 棵树")

        # 提取并验证每棵树
        tree_pattern = re.compile(
            r"TREE\s+(?:\*\s+)?(\w+)\s*=\s*(.*?);", re.IGNORECASE | re.DOTALL
        )
        for i, match in enumerate(tree_pattern.finditer(content)):
            tree_name = match.group(1)
            tree_str = match.group(2).strip()

            # 移除可能的注释
            tree_str = re.sub(r"\[&&[^\]]*\]", "", tree_str)
            tree_str = tree_str.strip()

            if not tree_str:
                warnings.append(f"树 '{tree_name}' 为空")
                continue

            prefix = f"树 '{tree_name}'" if tree_count > 1 else "树"

            # 验证树内容
            bracket_errors = self._check_bracket_balance(tree_str, prefix)
            errors.extend(bracket_errors)

            length_critical, length_warn = self._check_branch_lengths(tree_str, prefix)
            if length_critical:
                from ..core.exceptions import NegativeBranchLengthError

                raise NegativeBranchLengthError(
                    f"{prefix}包含 {len(length_critical)} 个负分支长度",
                    suggestion="请检查树文件中分支长度是否正确，或使用 --skip-length-check 跳过此检查。",
                )
            warnings.extend(length_warn)

            node_errors, node_warnings = self._check_node_names(
                tree_str, prefix, ignore_malformed=ignore_malformed
            )
            errors.extend(node_errors)
            warnings.extend(node_warnings)

            # 多根节点和自循环边检测 (CRITICAL)
            critical_errors = self._check_multi_root_and_self_loop(tree_str, prefix)
            if critical_errors:
                from ..core.exceptions import TreeValidationError

                raise TreeValidationError(
                    "；".join(critical_errors),
                    suggestion="请检查树文件是否包含多棵独立树或自循环边（节点连接到自身）。",
                )

            summary = self._get_tree_summary(tree_str, i + 1)
            summary["name"] = tree_name
            warnings.extend(self._append_ultrametric_summary(summary, tree_str, prefix))
            tree_summaries.append(summary)

        return TreeValidationDetail(
            is_valid=len(errors) == 0,
            format_detected="nexus",
            tree_count=tree_count,
            tree_summaries=tree_summaries,
            errors=errors,
            warnings=warnings,
        )

    def _check_branch_lengths(
        self, tree_str: str, prefix: str
    ) -> Tuple[List[str], List[str]]:
        """检查分支长度：负值 = CRITICAL，零长枝/异常大值 = WARNING。

        审阅项 C-25：此前用正则 ``r":(-?\\d+\\.?\\d*...)"`` 直接扫**原始 Newick
        文本**，实测既有假阳性（``("A:-3.5":0.25,B:0.25);`` 里引号内的叶名被当成
        负分支长度，合法输入被判 CRITICAL 并拦下）也有假阴性（``:-.5``、``:.5``、
        ``:+0.5`` 都不匹配，而 ABN Newick 的 ``<real>`` 明确允许 ``.<digits>``，
        也就是**真正负的** ``:-.5`` 会被漏过）。现在长度一律从**已解析的树**里读。

        审阅项 C-26：``warnings`` 通道此前在函数体内从未被追加（恒空），
        调用方以为拿到了一档"软问题"。现在零长枝与超长枝会真的走这条通道。

        Returns:
            (critical_errors, warnings): CRITICAL 级别错误和警告
        """
        critical_errors: List[str] = []
        warnings: List[str] = []

        lengths = self._branch_lengths_from_parsed_tree(tree_str)
        if lengths is None:
            # 解析不出来（通常是更基本的语法问题，括号平衡等检查会同时报错）。
            # 这里退回旧的文本扫描，但**明确告知**这一档检查是降级进行的。
            warnings.append(
                f"{prefix}: 无法从解析后的树对象读取分支长度，已退回文本层扫描"
                f"（对 ':.5' / 引号内叶名这类写法可能误判或漏判）"
            )
            lengths = self._branch_lengths_from_text(tree_str)

        zero_count = 0
        for label, length in lengths:
            if length is None:
                continue
            if length < 0:
                critical_errors.append(
                    f"[CRITICAL] {prefix}: 负分支长度 {length} 在节点 {label}"
                )
            elif length == 0:
                zero_count += 1
            elif length > 1e6:
                # 单位混淆（把 Ma 当成 substitutions/site 之类）最常见的表现
                warnings.append(
                    f"{prefix}: 节点 {label} 的分支长度异常大 ({length:g})，"
                    f"请确认单位与输入类型（时间树应使用与所定年单位一致的长度）"
                )
        if zero_count:
            warnings.append(
                f"{prefix}: 发现 {zero_count} 条零长度分支：超度量性/定年区间可能被"
                f"人为压平（部分后端会拒绝零长枝，或把该节点年龄算成与父节点相同）"
            )

        return critical_errors, warnings

    def _parse_newick_for_structure(
        self, tree_str: str
    ) -> Optional[Tuple[str, Any, Any]]:
        """把 Newick 文本解析成 ``(后端名, 根节点, 树对象)``。

        后两项标成 ``Any`` 而不是 ``object``：它们其实是 DendroPy 的
        ``Tree``/``Node`` 或 Bio.Phylo 的 ``Tree``/``Clade``（两个库都没有
        可用的类型标注），而调用方必须按后端类型分头用它们的方法。标成
        ``object`` 会让那些方法调用在静态视角下全部不存在。

        先试 DendroPy（本项目的核心解析后端），失败再试 Bio.Phylo —— 引号包裹的
        节点名（``("sp A":0.25,B:0.25)``）在 DendroPy 5.0 的 Newick 词法里会被拒，
        而那是完全合法的 Newick（审阅项 C-25 的假阳性就来自引号内的名字）。
        两个后端都失败时返回 None，调用方据此降级。
        """
        text = tree_str.strip()
        if not text:
            return None
        if not text.endswith(";"):
            text += ";"
        errors: List[str] = []

        try:
            from dendropy import Tree as DendropyTree

            parsed = DendropyTree.get(
                data=text,
                schema="newick",
                preserve_underscores=True,
                rooting="force-rooted",
            )
            return ("dendropy", parsed.seed_node, parsed)
        except Exception as e:
            errors.append(f"dendropy: {e}")

        try:
            from io import StringIO

            from Bio import Phylo as BioPhylo

            parsed = BioPhylo.read(StringIO(text), "newick")
            return ("biopython", parsed.root, parsed)
        except Exception as e:
            errors.append(f"biopython: {e}")

        self.logger.debug(f"树结构解析失败（双后端）: {' | '.join(errors)}")
        return None

    def _branch_lengths_from_parsed_tree(
        self, tree_str: str
    ) -> Optional[List[Tuple[str, Optional[float]]]]:
        """从**已解析的树**读取 (节点名, 分支长度)。解析失败返回 None。"""
        parsed = self._parse_newick_for_structure(tree_str)
        if parsed is None:
            return None
        backend, root, tree_obj = parsed

        out: List[Tuple[str, Optional[float]]] = []
        if backend == "dendropy":
            for node in tree_obj.postorder_node_iter():
                if node is root:
                    continue
                label = str(node.taxon.label) if node.taxon is not None else node.label
                out.append((label or "internal", node.edge.length))
            return out

        for clade in root.get_terminals() + root.get_nonterminals():
            if clade is root:
                continue
            out.append((clade.name or "internal", clade.branch_length))
        return out

    @staticmethod
    def _branch_lengths_from_text(
        tree_str: str,
    ) -> List[Tuple[str, Optional[float]]]:
        """文本层兜底扫描（仅在树无法解析时使用；见 :meth:`_check_branch_lengths`）。

        比旧实现多支持 ``:.5`` / ``:-.5`` 这类省略整数位的形式（ABN Newick 的
        ``<real>`` 允许），并跳过引号包裹的节点名，减少假阳性。
        """
        results: List[Tuple[str, Optional[float]]] = []
        pattern = re.compile(r":(-?(?:\d+)?(?:\.\d+)?(?:[eE][-+]?\d+)?)")
        stripped = re.sub(r'"(?:[^"]|"")*"', '"NAME"', tree_str)
        for match in pattern.finditer(stripped):
            try:
                results.append((f"位置 {match.start()}", float(match.group(1))))
            except ValueError:
                continue
        return results

    # ---------------- 超度量性 / 时间树判定（B-18） ---------------- #

    def check_ultrametricity(
        self,
        tree: Union["PhylogeneticTree", Path, str, None],
        tolerance: Optional[float] = None,
    ) -> UltrametricDetail:
        """判定输入树是否为**时间树**（超度量：所有叶节点到根的距离相等）。

        审阅项 B-18：整个包此前没有任何一处做这个检查（``ultrametric`` 等关键词
        在 53 个 .py 文件里零命中），而 LSD2/IQ-TREE2 尖端定年、MCMCTree 全局钟
        与 ``viz/tree_plot._compute_node_age`` 的 ``root_age - dist(root, node)``
        公式都把它当前提。

        Args:
            tree: :class:`~phylodater.models.tree.PhylogeneticTree`、Newick 字符串
                或树文件路径（只读，不要求任何外部树后端）
            tolerance: root-to-tip **相对极差** ``(max-min)/max`` 的容差，
                默认 :data:`ULTRAMETRIC_DEFAULT_TOLERANCE`

        Returns:
            UltrametricDetail
        """
        tol = ULTRAMETRIC_DEFAULT_TOLERANCE if tolerance is None else float(tolerance)
        newick, notes = self._as_newick_text(tree)
        if newick is None:
            return UltrametricDetail(
                assessable=False, is_time_tree=False, tolerance=tol, notes=notes
            )

        heights, notes2 = self._root_to_tip_heights(newick)
        notes = notes + notes2
        if heights is None:
            return UltrametricDetail(
                assessable=False,
                is_time_tree=False,
                tolerance=tol,
                notes=notes,
            )

        min_height = min(heights.values())
        max_height = max(heights.values())
        if max_height <= 0:
            notes.append(
                "所有 root-to-tip 距离都为 0（无分支长度或纯拓扑树），无法评估超度量性"
            )
            return UltrametricDetail(
                assessable=False,
                is_time_tree=False,
                n_tips=len(heights),
                tolerance=tol,
                notes=notes,
            )

        spread = (max_height - min_height) / max_height
        offenders = sorted(
            ((name, (max_height - h) / max_height) for name, h in heights.items()),
            key=lambda item: item[1],
            reverse=True,
        )
        return UltrametricDetail(
            assessable=True,
            is_time_tree=spread <= tol,
            n_tips=len(heights),
            min_height=min_height,
            max_height=max_height,
            relative_spread=spread,
            tolerance=tol,
            worst_tips=offenders[:5],
            notes=notes,
        )

    def _as_newick_text(
        self, tree: Union["PhylogeneticTree", Path, str, None]
    ) -> Tuple[Optional[str], List[str]]:
        """把 PhylogeneticTree / 路径 / 字符串统一成 Newick 文本。"""
        notes: List[str] = []
        if tree is None:
            return None, ["输入树为空"]
        newick = getattr(tree, "newick", None)
        if isinstance(newick, str) and newick.strip():
            return newick, notes
        if isinstance(tree, Path):
            try:
                return self._read_file_safe(tree), notes
            except Exception as e:
                return None, [f"无法读取树文件 {tree}: {e}"]
        if isinstance(tree, str):
            return tree, notes
        # Mock/其它对象：没有可用的 newick 文本
        return None, ["输入对象不提供 newick 文本，无法判定超度量性"]

    def _root_to_tip_heights(
        self, newick: str
    ) -> Tuple[Optional[Dict[str, float]], List[str]]:
        """根到每个叶节点的距离（走解析后的树，不扫文本；双后端）。"""
        notes: List[str] = []
        text = newick.strip()
        if not text:
            return None, ["Newick 文本为空"]

        parsed = self._parse_newick_for_structure(text)
        if parsed is None:
            return None, ["树无法解析，无法判定超度量性（DendroPy/BioPython 均失败）"]
        backend, root, _tree = parsed

        if backend == "dendropy":
            is_leaf = lambda node: node.is_leaf()  # noqa: E731
            children_of = lambda node: node.child_nodes()  # noqa: E731
            length_of = lambda node: node.edge.length  # noqa: E731
            label_of = lambda node: (  # noqa: E731
                str(node.taxon.label) if node.taxon is not None else (node.label or "?")
            )
        else:
            is_leaf = lambda clade: clade.is_terminal()  # noqa: E731
            children_of = lambda clade: list(clade.clades)  # noqa: E731
            length_of = lambda clade: clade.branch_length  # noqa: E731
            label_of = lambda clade: (clade.name or "?")  # noqa: E731

        heights: Dict[str, float] = {}
        missing = 0
        non_finite = 0
        # 显式栈，避免超深树触发 RecursionError
        stack: List[Tuple[object, float, int]] = [(root, 0.0, 0)]
        while stack:
            node, acc, depth = stack.pop()
            if depth > 10000:
                return None, ["树过深，超度量性判定被跳过（深度上限 10000）"]
            if is_leaf(node):
                heights[label_of(node)] = acc
                continue
            for child in children_of(node):
                length = length_of(child)
                if length is None:
                    missing += 1
                    length = 0.0
                elif not math.isfinite(length):
                    # NaN/inf 让 root-to-tip 距离失去意义（B-19 的同类形态）
                    non_finite += 1
                stack.append((child, acc + length, depth + 1))

        if len(heights) < 2:
            return None, ["叶节点少于 2 个，无法判定超度量性"]
        if non_finite:
            notes.append(
                f"有 {non_finite} 条分支长度不是有限数值（NaN/inf），"
                f"root-to-tip 距离不可信"
            )
            return None, notes
        if missing:
            notes.append(
                f"有 {missing} 条分支没有长度信息（按 0 处理）——"
                f"这说明输入很可能不是时间树，甚至是纯拓扑树"
            )
        if not heights:
            return None, ["未能计算出任何 root-to-tip 距离"]
        return heights, notes

    @staticmethod
    def ultrametric_severity(
        method: Optional[str], require_time_tree: bool = False
    ) -> str:
        """按方法决定"输入树非超度量"该以什么级别报告。

        - ``info``：treePL / PATHd8 / r8s / MD-Cat / wLogDate 等**设计上就允许**
          非超度量输入（速率平滑或相对定年），不该拦用户；
        - ``warning``：LSD2 / IQ-TREE2 尖端定年、MCMCTree 全局钟以及 viz 的
          节点年龄公式把这些当前提；
        - ``error``：调用方显式 ``require_time_tree=True``。
        """
        if require_time_tree:
            return "error"
        key = (method or "").strip().lower()
        if not key:
            return "warning"
        if key in ULTRAMETRIC_TOLERANT_METHODS:
            return "info"
        if key in TIME_TREE_ASSUMING_METHODS:
            return "warning"
        return "warning"

    def validate_ultrametricity(
        self,
        tree: Union["PhylogeneticTree", Path, str, None],
        method: Optional[str] = None,
        tolerance: Optional[float] = None,
        require_time_tree: bool = False,
    ) -> Tuple[List[str], List[str], List[str]]:
        """超度量性检查的报告分诊：返回 ``(errors, warnings, infos)``。

        默认**永不**产生 error（B-18 要求 WARN 而非 hard-fail）；只有显式
        ``require_time_tree=True`` 时才升级为 error。
        """
        detail = self.check_ultrametricity(tree, tolerance=tolerance)
        message = detail.summarize()
        severity = self.ultrametric_severity(method, require_time_tree)

        if not detail.assessable:
            # 无法判定：只报 warning（不能说"没问题"，也不能冤枉用户）
            return [], [f"{message}（方法: {method or '未指定'}）"], []
        if detail.is_time_tree:
            return [], [], [message]

        label = f"（方法: {method or '未指定'}）"
        if severity == "info":
            return (
                [],
                [],
                [
                    f"{message}{label} — 该方法允许非超度量输入（速率平滑/相对定年），"
                    f"仅作告知"
                ],
            )
        if severity == "error":
            return [f"{message}{label}"], [], []
        return (
            [],
            [
                f"{message}{label} — 该通路把'输入树已是时间树'当前提"
                f"（LSD2/IQ-TREE2 尖端定年、MCMCTree 全局钟、以及可视化层"
                f" root_age − 根距 的节点年龄公式）。非超度量输入会让末端枝被强拉到"
                f" present、内部节点与其子枝在时间轴上错位。"
                f"请先用 treePL/PATHd8/r8s/bigg-me 一类方法或 --rooting 生成时间树，"
                f"或提高 calibration 的容差设置"
            ],
            [],
        )

    def _append_ultrametric_summary(
        self, summary: Dict, tree_str: str, prefix: str
    ) -> List[str]:
        """把超度量性判定写进树摘要，并返回要上报的 warning 列表。

        判定不出来时**不会**声称"是时间树"，也不会因此判输入无效（B-18 要求
        对允许非超度量的方法只提示）。
        """
        detail = self.check_ultrametricity(tree_str)
        summary["is_ultrametric"] = (
            None if not detail.assessable else detail.is_time_tree
        )
        summary["root_to_tip_relative_spread"] = detail.relative_spread
        summary["root_to_tip_range"] = (
            None
            if detail.min_height is None or detail.max_height is None
            else (detail.min_height, detail.max_height)
        )

        if not detail.assessable:
            return [f"{prefix}: {detail.summarize()}"]
        if detail.is_time_tree:
            return []
        return [
            f"{prefix}: {detail.summarize()}",
            f"{prefix}: 输入树不是时间树——treePL/PATHd8/r8s 可以直接使用它，"
            f"但 LSD2/IQ-TREE2 尖端定年、MCMCTree 全局钟与可视化层的"
            f" 'root_age − 根距' 节点年龄公式都以它超度量为前提",
        ]

    def _check_bracket_balance(self, tree_str: str, prefix: str) -> List[str]:
        """检查括号平衡"""
        errors = []

        # 检查括号平衡（考虑引号内的括号）
        in_quote = False
        quote_char = None
        paren_depth = 0

        for i, char in enumerate(tree_str):
            if in_quote:
                if char == quote_char and (i == 0 or tree_str[i - 1] != "\\"):
                    in_quote = False
                continue

            if char in ("'", '"'):
                in_quote = True
                quote_char = char
            elif char == "(":
                paren_depth += 1
            elif char == ")":
                paren_depth -= 1
                if paren_depth < 0:
                    errors.append(f"{prefix}: 括号不匹配 - 多余的右括号 ')' 在位置 {i}")
                    break

        if in_quote:
            errors.append(f"{prefix}: 引号未正确关闭")
        elif paren_depth > 0:
            errors.append(f"{prefix}: 括号不匹配 - 缺少 {paren_depth} 个右括号 ')'")
        elif paren_depth < 0:
            errors.append(f"{prefix}: 括号不匹配 - 多余 {-paren_depth} 个右括号 ')'")

        return errors

    def _check_node_names(
        self, tree_str: str, prefix: str, ignore_malformed: Optional[bool] = None
    ) -> Tuple[List[str], List[str]]:
        """检查节点名称"""
        errors = []
        warnings = []

        # 提取所有节点名（在括号外、逗号前或冒号前的名称）
        # 简化模式：匹配 ,(name): 或 ,(name), 或 (name): 等
        names = []
        for match in _NODE_NAME_RE.finditer(tree_str):
            name = _strip_node_name(match.group(1)).strip()
            if name:
                names.append(name)

        # 检查空节点名
        # 匹配连续的逗号或括号，表示空节点名
        empty_name_pattern = re.compile(r",[,\s]*,|\(\s*[:,)]|\,\s*\)|\(\s*\)")
        if empty_name_pattern.search(tree_str):
            errors.append(f"{prefix}: 发现空节点名")

        # 区分末端节点和内部节点的重复名
        # 末端节点名：在 Newick 中不出现在 ')' 之后的名称
        # 内部节点名：出现在 ')' 之后、':' 之前的名称
        internal_names = []
        for match in _INTERNAL_NODE_RE.finditer(tree_str):
            iname = _strip_node_name(match.group(1)).strip()
            if iname:
                internal_names.append(iname)

        # 末端节点名 = 所有名称 - 内部节点名
        tip_names = [n for n in names if n not in internal_names]

        # 检查末端节点重复（ERROR）
        seen_tips = set()
        duplicate_tips = set()
        for name in tip_names:
            if name in seen_tips:
                duplicate_tips.add(name)
            seen_tips.add(name)

        if duplicate_tips:
            dup_list = sorted(list(duplicate_tips))
            errors.append(
                f"{prefix}: 发现 {len(duplicate_tips)} 个重复的末端节点名 (ERROR): "
                f"{', '.join(dup_list[:5])}" + ("..." if len(dup_list) > 5 else "")
            )

        # 检查内部节点重复（WARNING，含嵌套层级）
        seen_internals = set()
        duplicate_internals = set()
        for name in internal_names:
            if name in seen_internals:
                duplicate_internals.add(name)
            seen_internals.add(name)

        if duplicate_internals:
            # 计算每个重复内部节点名出现的括号嵌套层级
            dup_depths: Dict[str, Set[int]] = {n: set() for n in duplicate_internals}
            for match in _INTERNAL_NODE_RE.finditer(tree_str):
                iname = match.group(1).strip()
                if iname in duplicate_internals:
                    pos = match.start()
                    # 计算该位置的括号嵌套深度（未闭合的 '(' 数量）
                    depth = 0
                    for i in range(pos):
                        if tree_str[i] == "(":
                            depth += 1
                        elif tree_str[i] == ")":
                            depth -= 1
                    dup_depths[iname].add(depth)
            dup_list = sorted(list(duplicate_internals))
            depth_info = ", ".join(
                f"{n}(嵌套层级:{'/'.join(str(d) for d in sorted(dup_depths[n]) if d >= 0) or '0'})"
                for n in dup_list[:5]
            )
            warnings.append(
                f"{prefix}: 发现 {len(duplicate_internals)} 个重复的内部节点名 (WARNING): "
                f"{depth_info}" + ("..." if len(dup_list) > 5 else "")
            )

        # 检查超长节点名
        for match in _NODE_NAME_RE.finditer(tree_str):
            name = _strip_node_name(match.group(1))
            if len(name) > 4096:
                warnings.append(
                    f"{prefix}: 节点名长度超过 4096 字节 ({len(name)} 字节)，"
                    "可能影响部分下游软件兼容性"
                )
                break

        # 检查恶意符号注入
        malicious_errors, malicious_warnings = self._check_malicious_chars(
            tree_str, prefix, ignore_malformed=ignore_malformed
        )
        errors.extend(malicious_errors)
        warnings.extend(malicious_warnings)

        return errors, warnings

    def _check_multi_root_and_self_loop(self, tree_str: str, prefix: str) -> List[str]:
        """
        检测多根节点和自循环边（CRITICAL 级别）

        检查项目：
        1. 多根检测：Newick 字符串中存在多个以 ``;`` 结尾的顶层括号组（多棵独立树）
        2. 自循环边检测：存在零长度自环模式 ``name:0.0,name``

        Args:
            tree_str: 单棵树的 Newick 字符串
            prefix: 错误前缀（如 "树 1"）

        Returns:
            List[str]: CRITICAL 错误消息列表（带 ``[CRITICAL]`` 前缀）
        """
        critical_errors = []

        # 多根检测：统计顶层树数量（以分号结尾的独立树）
        stripped = tree_str.strip()
        top_level_trees = [t for t in re.split(r";\s*", stripped) if t.strip()]
        if len(top_level_trees) > 1:
            critical_errors.append(
                f"[CRITICAL] {prefix}: 检测到多根节点（{len(top_level_trees)} 棵独立树），"
                f"Newick 格式应仅含一棵树"
            )

        # 自循环边检测：零长度自环模式 (name:0.0,name)
        self_loop_pattern = _SELF_LOOP_RE
        self_loop_matches = list(self_loop_pattern.finditer(tree_str))
        if self_loop_matches:
            positions = [str(m.start()) for m in self_loop_matches[:5]]
            critical_errors.append(
                f"[CRITICAL] {prefix}: 检测到自循环边（零长度自环）({len(self_loop_matches)} 处): "
                f"位置 {', '.join(positions)}"
                + ("..." if len(self_loop_matches) > 5 else "")
            )

        return critical_errors

    def _check_malicious_chars(
        self, text: str, prefix: str, ignore_malformed: Optional[bool] = None
    ) -> Tuple[List[str], List[str]]:
        """
        检测恶意符号注入

        检查项目：
        1. 控制字符（\x00-\x1f, \x7f）
        2. Unicode 双向文本控制字符（\u202a-\u202e, \u2066-\u2069）
        3. 零宽字符（\u200b-\u200f, \u2028-\u2029, \u2060-\u2064, \ufeff）
        4. 其他可疑 Unicode 字符

        ``ignore_malformed`` 为 True 时，检测到的恶意字符降级为警告（写入 warnings
        而非 errors）；传 ``None`` 表示沿用实例级默认值 ``self.ignore_malformed``。

        Returns:
            (errors, warnings): 错误列表和警告列表
        """
        errors = []
        warnings = []
        effective = (
            self.ignore_malformed if ignore_malformed is None else ignore_malformed
        )

        # 控制字符模式（除了常见的 \t, \n, \r）
        control_char_pattern = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]")

        # Unicode 双向文本控制字符
        bidi_pattern = re.compile(r"[\u202a-\u202e\u2066-\u2069]")

        # 零宽字符
        zero_width_pattern = re.compile(
            r"[\u200b-\u200f\u2028-\u2029\u2060-\u2064\ufeff]"
        )

        # 检查控制字符
        control_matches = list(control_char_pattern.finditer(text))
        if control_matches:
            positions = [
                f"0x{ord(m.group()):02x}@{m.start()}" for m in control_matches[:5]
            ]
            msg = (
                f"{prefix}: 发现控制字符注入 ({len(control_matches)} 处): "
                f"{', '.join(positions)}" + ("..." if len(control_matches) > 5 else "")
            )
            if effective:
                warnings.append(msg)
            else:
                errors.append(msg)

        # 检查双向文本控制字符
        bidi_matches = list(bidi_pattern.finditer(text))
        if bidi_matches:
            positions = [
                f"\\u{ord(m.group()):04x}@{m.start()}" for m in bidi_matches[:5]
            ]
            msg = (
                f"{prefix}: 发现 Unicode 双向文本控制字符 ({len(bidi_matches)} 处): "
                f"{', '.join(positions)}" + ("..." if len(bidi_matches) > 5 else "")
            )
            if effective:
                warnings.append(msg)
            else:
                errors.append(msg)

        # 检查零宽字符
        zero_width_matches = list(zero_width_pattern.finditer(text))
        if zero_width_matches:
            positions = [
                f"\\u{ord(m.group()):04x}@{m.start()}" for m in zero_width_matches[:5]
            ]
            msg = (
                f"{prefix}: 发现零宽字符注入 ({len(zero_width_matches)} 处): "
                f"{', '.join(positions)}"
                + ("..." if len(zero_width_matches) > 5 else "")
            )
            if effective:
                warnings.append(msg)
            else:
                errors.append(msg)

        return errors, warnings

    def _get_tree_summary(self, tree_str: str, tree_num: int) -> Dict:
        """获取树的统计摘要"""
        # 统计叶节点数（简化：统计逗号分隔的名称数量 + 1）
        # 更准确的方法是解析树结构

        # 统计括号对数（大致估计内部节点数）
        internal_count = tree_str.count("(")

        # 统计叶节点（名称数量）
        name_pattern = re.compile(r"[,(]([A-Za-z_][A-Za-z0-9_.]*?)(?=[,:)\s])")
        names = [m.group(1) for m in name_pattern.finditer(tree_str) if m.group(1)]
        terminal_count = len(set(names)) if names else 0

        # 如果没有解析到名称，使用更简单的方法
        if terminal_count == 0:
            # 统计逗号数量 + 1 作为叶节点数的粗略估计
            terminal_count = tree_str.count(",") + 1

        return {
            "tree_num": tree_num,
            "terminals": terminal_count,
            "internal": internal_count,
            "length": len(tree_str),
        }

    # ==================== 序列文件验证 ====================

    def validate_alignment_deep(
        self,
        file_path: Path,
        expected_alphabet: Optional[SequenceAlphabet] = None,
        require_aligned: bool = True,
        ignore_malformed: Optional[bool] = None,
    ) -> SequenceValidationDetail:
        """
        深度验证序列文件

        检查项目：
        1. 格式自动识别（FASTA/FASTQ）
        2. 序列字母表检查（DNA/RNA/protein）
        3. 无效字符定位（行号）
        4. 序列ID唯一性检查
        5. 序列长度一致性检查

        Args:
            file_path: 序列文件路径
            expected_alphabet: 期望的字母表类型（None 表示自动检测）
            require_aligned: 是否要求序列已对齐
            ignore_malformed: 是否将恶意字符/格式错误降级为警告（None 表示使用实例默认值）

        Returns:
            SequenceValidationDetail: 验证详情
        """
        effective_ignore = (
            ignore_malformed if ignore_malformed is not None else self.ignore_malformed
        )
        if not file_path.exists():
            return SequenceValidationDetail(
                is_valid=False, errors=[f"序列文件不存在: {file_path}"]
            )

        if file_path.stat().st_size == 0:
            return SequenceValidationDetail(
                is_valid=False, errors=[f"序列文件为空: {file_path}"]
            )

        if file_path.stat().st_size > self.MAX_INPUT_SIZE:
            from ..core.exceptions import InputSizeLimitError

            raise InputSizeLimitError(
                f"序列文件过大 ({file_path.stat().st_size / (1024*1024):.0f} MB), "
                f"超过限制 ({self.MAX_INPUT_SIZE // (1024*1024)} MB)"
            )

        # 检测格式
        fmt = self._detect_sequence_format(file_path)

        if fmt == "fasta":
            return self._validate_fasta_deep(
                file_path, expected_alphabet, require_aligned
            )
        elif fmt == "fastq":
            return self._validate_fastq_deep(file_path, expected_alphabet)
        elif fmt in ("phylip", "phylip_interleaved", "phylip_sequential"):
            if effective_ignore:
                return SequenceValidationDetail(
                    is_valid=True,
                    format_detected=fmt,
                    warnings=[
                        f"PHYLIP 格式 ({fmt}) 不被支持，但 --ignore-malformed 已启用，"
                        f"跳过深度验证: {file_path}"
                    ],
                )
            from ..core.exceptions import PhyloFormatError

            raise PhyloFormatError(
                f"PHYLIP 格式 ({fmt}) 不被支持: {file_path}",
                suggestion="请将序列文件转换为 FASTA 格式。可使用 seqmagick 或 bioawk 工具进行格式转换。",
            )
        elif fmt == "stockholm":
            if effective_ignore:
                return SequenceValidationDetail(
                    is_valid=True,
                    format_detected=fmt,
                    warnings=[
                        f"Stockholm 格式不被支持，但 --ignore-malformed 已启用，"
                        f"跳过深度验证: {file_path}"
                    ],
                )
            from ..core.exceptions import PhyloFormatError

            raise PhyloFormatError(
                f"Stockholm 格式不被支持: {file_path}",
                suggestion="请将序列文件转换为 FASTA 格式。可使用 esl-reformat 或 seqmagick 工具进行格式转换。",
            )
        else:
            from ..core.exceptions import PhyloFormatError

            raise PhyloFormatError(
                f"无法识别的序列文件格式: {file_path}",
                suggestion="支持的格式: FASTA (.fasta, .fa), FASTQ (.fastq, .fq)",
            )

    def _detect_sequence_format(self, file_path: Path) -> str:
        """检测序列文件格式"""
        ext = file_path.suffix.lower()

        # 根据扩展名判断
        if ext in [".fasta", ".fa", ".fna", ".faa", ".fas", ".afa"]:
            return "fasta"
        elif ext in [".fastq", ".fq"]:
            return "fastq"

        # 根据首字符判断
        try:
            with open(file_path, "r", encoding="utf-8", errors="strict") as f:
                first_char = f.read(1)

            if first_char == ">":
                return "fasta"
            elif first_char == "@":
                return "fastq"
        except UnicodeDecodeError as e:
            self.logger.warning(
                f"文件编码问题，无法以 UTF-8 严格模式读取: {file_path} - {e}"
            )
        except OSError:
            pass

        return "unknown"

    def _validate_fasta_deep(
        self,
        file_path: Path,
        expected_alphabet: Optional[SequenceAlphabet],
        require_aligned: bool,
    ) -> SequenceValidationDetail:
        """深度验证 FASTA 格式"""
        errors = []
        warnings = []
        duplicate_ids = []
        invalid_chars = []
        sequence_lengths = []
        detected_alphabet_counts = {"dna": 0, "rna": 0, "protein": 0, "unknown": 0}

        seen_ids: Set[str] = set()
        current_id = None
        current_seq_lines: List[str] = []
        line_number = 0

        try:
            with open(file_path, "r", encoding="utf-8", errors="strict") as f:
                for line in f:
                    line_number += 1
                    line = line.strip()

                    if not line:
                        continue

                    if line.startswith(">"):
                        # 保存前一个序列
                        if current_id is not None:
                            seq = "".join(current_seq_lines)
                            sequence_lengths.append(len(seq))

                            # 检查字母表
                            seq_alphabet = self._detect_sequence_alphabet(seq)
                            detected_alphabet_counts[seq_alphabet] = (
                                detected_alphabet_counts.get(seq_alphabet, 0) + 1
                            )

                            # 检查无效字符
                            if expected_alphabet:
                                invalid = self._find_invalid_chars(
                                    seq, expected_alphabet
                                )
                                for inv_char, inv_pos in invalid:
                                    invalid_chars.append(
                                        {
                                            "line": line_number
                                            - len(current_seq_lines),
                                            "id": current_id,
                                            "char": inv_char,
                                            "position": inv_pos,
                                        }
                                    )

                        # 解析新序列ID
                        header = line[1:].strip()
                        current_id = header.split()[0] if header else ""

                        # 检查ID唯一性
                        if current_id in seen_ids:
                            duplicate_ids.append(current_id)
                        seen_ids.add(current_id)

                        current_seq_lines = []
                    else:
                        current_seq_lines.append(line)

                # 处理最后一个序列
                if current_id is not None:
                    seq = "".join(current_seq_lines)
                    sequence_lengths.append(len(seq))

                    seq_alphabet = self._detect_sequence_alphabet(seq)
                    detected_alphabet_counts[seq_alphabet] = (
                        detected_alphabet_counts.get(seq_alphabet, 0) + 1
                    )

                    if expected_alphabet:
                        invalid = self._find_invalid_chars(seq, expected_alphabet)
                        for inv_char, inv_pos in invalid:
                            invalid_chars.append(
                                {
                                    "line": line_number - len(current_seq_lines) + 1,
                                    "id": current_id,
                                    "char": inv_char,
                                    "position": inv_pos,
                                }
                            )

        except UnicodeDecodeError as e:
            self.logger.warning(
                f"FASTA 文件编码问题，无法以 UTF-8 严格模式读取: {file_path} - {e}"
            )
            return SequenceValidationDetail(
                is_valid=False, errors=[f"解析 FASTA 文件失败（编码错误）: {e}"]
            )
        except Exception as e:
            return SequenceValidationDetail(
                is_valid=False, errors=[f"解析 FASTA 文件失败: {e}"]
            )

        num_sequences = len(sequence_lengths)
        if num_sequences == 0:
            return SequenceValidationDetail(
                is_valid=False,
                format_detected="fasta",
                errors=[f"FASTA 文件中没有找到序列: {file_path}"],
            )

        # 检测主要字母表：``dict.get`` 的返回类型带 None，用作 key 函数时
        # mypy 无法确认它可比；这里的键一定在字典里，直接下标更贴合事实。
        detected_alphabet = max(
            detected_alphabet_counts,
            key=lambda k: detected_alphabet_counts[k],
        )

        # 验证字母表一致性
        if expected_alphabet:
            expected_str = expected_alphabet.value
            if detected_alphabet != expected_str and detected_alphabet != "unknown":
                warnings.append(
                    f"检测到的字母表 ({detected_alphabet}) 与期望 ({expected_str}) 不一致"
                )

        # 检查重复ID
        if duplicate_ids:
            errors.append(
                f"发现 {len(duplicate_ids)} 个重复的序列ID: "
                f"{', '.join(duplicate_ids[:5])}"
                + ("..." if len(duplicate_ids) > 5 else "")
            )

        # 检查无效字符
        if invalid_chars:
            error_lines = []
            for inv in invalid_chars[:10]:  # 最多显示10个
                error_lines.append(
                    f"  行 {inv['line']}: 序列 '{inv['id']}' 中发现无效字符 "
                    f"'{inv['char']}' 在位置 {inv['position']}"
                )
            errors.append(
                f"发现 {len(invalid_chars)} 个无效字符:\n" + "\n".join(error_lines)
            )

        # 检查序列长度一致性
        is_aligned = len(set(sequence_lengths)) <= 1
        if require_aligned and not is_aligned:
            unique_lengths = sorted(set(sequence_lengths))
            warnings.append(
                f"序列长度不一致 (发现 {len(unique_lengths)} 种不同长度): "
                f"{unique_lengths[:5]}" + ("..." if len(unique_lengths) > 5 else "")
            )

        return SequenceValidationDetail(
            is_valid=len(errors) == 0,
            format_detected="fasta",
            alphabet_detected=detected_alphabet,
            num_sequences=num_sequences,
            sequence_length=sequence_lengths[0] if is_aligned else None,
            is_aligned=is_aligned,
            duplicate_ids=duplicate_ids,
            invalid_chars=invalid_chars,
            errors=errors,
            warnings=warnings,
        )

    def _validate_fastq_deep(
        self, file_path: Path, expected_alphabet: Optional[SequenceAlphabet]
    ) -> SequenceValidationDetail:
        """深度验证 FASTQ 格式"""
        errors = []
        warnings: List[str] = []
        duplicate_ids = []
        invalid_chars = []
        sequence_lengths = []
        detected_alphabet_counts = {"dna": 0, "rna": 0, "protein": 0, "unknown": 0}

        seen_ids: Set[str] = set()
        line_number = 0
        record_lines = []

        try:
            with open(file_path, "r", encoding="utf-8", errors="strict") as f:
                for line in f:
                    line_number += 1
                    line = line.rstrip("\n")
                    record_lines.append(line)

                    # FASTQ 记录由4行组成
                    if len(record_lines) == 4:
                        header = record_lines[0]
                        seq = record_lines[1]
                        plus = record_lines[2]
                        quality = record_lines[3]

                        # 验证 FASTQ 格式
                        if not header.startswith("@"):
                            errors.append(
                                f"行 {line_number - 3}: FASTQ 记录必须以 '@' 开头"
                            )
                        elif plus != "+":
                            errors.append(
                                f"行 {line_number - 1}: FASTQ 分隔符必须是 '+'"
                            )
                        elif len(seq) != len(quality):
                            errors.append(
                                f"行 {line_number - 2}: 序列长度 ({len(seq)}) "
                                f"与质量分数长度 ({len(quality)}) 不一致"
                            )
                        else:
                            # 解析ID
                            seq_id = header[1:].split()[0]

                            # 检查ID唯一性
                            if seq_id in seen_ids:
                                duplicate_ids.append(seq_id)
                            seen_ids.add(seq_id)

                            # 记录序列长度
                            sequence_lengths.append(len(seq))

                            # 检查字母表
                            seq_alphabet = self._detect_sequence_alphabet(seq)
                            detected_alphabet_counts[seq_alphabet] = (
                                detected_alphabet_counts.get(seq_alphabet, 0) + 1
                            )

                            # 检查无效字符
                            if expected_alphabet:
                                invalid = self._find_invalid_chars(
                                    seq, expected_alphabet
                                )
                                for inv_char, inv_pos in invalid:
                                    invalid_chars.append(
                                        {
                                            "line": line_number - 2,
                                            "id": seq_id,
                                            "char": inv_char,
                                            "position": inv_pos,
                                        }
                                    )

                        record_lines = []

        except UnicodeDecodeError as e:
            self.logger.warning(
                f"FASTQ 文件编码问题，无法以 UTF-8 严格模式读取: {file_path} - {e}"
            )
            return SequenceValidationDetail(
                is_valid=False, errors=[f"解析 FASTQ 文件失败（编码错误）: {e}"]
            )
        except Exception as e:
            return SequenceValidationDetail(
                is_valid=False, errors=[f"解析 FASTQ 文件失败: {e}"]
            )

        num_sequences = len(sequence_lengths)
        if num_sequences == 0:
            return SequenceValidationDetail(
                is_valid=False,
                format_detected="fastq",
                errors=[f"FASTQ 文件中没有找到有效记录: {file_path}"],
            )

        # 检测主要字母表（同上：用下标而不是 .get 做 key）
        detected_alphabet = max(
            detected_alphabet_counts,
            key=lambda k: detected_alphabet_counts[k],
        )

        # 检查重复ID
        if duplicate_ids:
            errors.append(
                f"发现 {len(duplicate_ids)} 个重复的序列ID: "
                f"{', '.join(duplicate_ids[:5])}"
                + ("..." if len(duplicate_ids) > 5 else "")
            )

        # 检查无效字符
        if invalid_chars:
            error_lines = []
            for inv in invalid_chars[:10]:
                error_lines.append(
                    f"  行 {inv['line']}: 序列 '{inv['id']}' 中发现无效字符 "
                    f"'{inv['char']}' 在位置 {inv['position']}"
                )
            errors.append(
                f"发现 {len(invalid_chars)} 个无效字符:\n" + "\n".join(error_lines)
            )

        return SequenceValidationDetail(
            is_valid=len(errors) == 0,
            format_detected="fastq",
            alphabet_detected=detected_alphabet,
            num_sequences=num_sequences,
            sequence_length=(
                sequence_lengths[0] if len(set(sequence_lengths)) == 1 else None
            ),
            is_aligned=False,  # FASTQ 通常不是对齐的
            duplicate_ids=duplicate_ids,
            invalid_chars=invalid_chars,
            errors=errors,
            warnings=warnings,
        )

    def _detect_sequence_alphabet(self, sequence: str) -> str:
        """检测序列的字母表类型"""
        if not sequence:
            return "unknown"

        # 统计字符
        chars = set(sequence.upper())

        # 检查是否全是 DNA 字符
        if chars.issubset(set("ATCGUN-?.*")):
            if "U" in chars:
                return "rna"
            return "dna"

        # 检查是否是蛋白质
        if chars.issubset(PROTEIN_CHARS):
            return "protein"

        return "unknown"

    def _find_invalid_chars(
        self, sequence: str, expected_alphabet: SequenceAlphabet
    ) -> List[Tuple[str, int]]:
        """查找序列中的无效字符"""
        invalid: List[Tuple[str, int]] = []

        if expected_alphabet == SequenceAlphabet.DNA:
            valid_chars = set("ATCGUNatcgun-?.")
        elif expected_alphabet == SequenceAlphabet.RNA:
            valid_chars = set("ACGUNacgun-?.")
        elif expected_alphabet == SequenceAlphabet.PROTEIN:
            valid_chars = PROTEIN_CHARS
        else:
            return invalid

        for i, char in enumerate(sequence):
            if char not in valid_chars:
                invalid.append((char, i + 1))  # 位置从1开始

        return invalid

    # ==================== 分类学表格验证 ====================

    @dataclass
    class TaxonomyValidationDetail:
        """分类学表格验证详情"""

        is_valid: bool
        circular_dependencies: List[Tuple[str, str]] = field(default_factory=list)
        errors: List[str] = field(default_factory=list)
        warnings: List[str] = field(default_factory=list)
        #: 实际进入循环依赖判定的行数（0 = 本表没有任何可解析的分类学字符串，
        #: 检查覆盖为零，必须说清楚，不能让调用方以为"查过且干净"）
        entries_checked: int = 0

    def validate_taxonomy_file(self, file_path: Path) -> "TaxonomyValidationDetail":
        """
        验证分类学表格文件

        检查项目：
        1. 文件存在性和大小
        2. 循环依赖检测
        3. 格式验证

        审阅项 C-39：本方法此前**没有任何生产调用点**（CLI 只调
        :meth:`InputValidator.validate_taxonomy_file` 的浅检查），现在由该方法在
        结构检查通过后调用一次，因此这里的返回契约（errors / warnings /
        entries_checked）都会被翻译成验证门的结果。

        Args:
            file_path: 分类学表格文件路径

        Returns:
            TaxonomyValidationDetail: 验证详情

        Raises:
            TaxonomyConflictError: 检测到分类学名称在不同层级间循环引用
        """
        if not file_path.exists():
            return self.TaxonomyValidationDetail(
                is_valid=False, errors=[f"分类学文件不存在: {file_path}"]
            )

        if file_path.stat().st_size == 0:
            return self.TaxonomyValidationDetail(
                is_valid=False, errors=[f"分类学文件为空 (CRITICAL): {file_path}"]
            )

        from ..core.exceptions import TaxonomyConflictError

        try:
            # 读取文件（与浅检查一致地按 UTF-8 读；失败由下方 except 汇总成 errors）
            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            if not lines:
                return self.TaxonomyValidationDetail(
                    is_valid=False, errors=[f"分类学文件为空: {file_path}"]
                )

            # 解析分类学数据
            taxonomy_data = self._parse_taxonomy_lines(lines)

            # 检测循环依赖
            circular_deps = self._detect_circular_dependencies(taxonomy_data)

            errors: List[str] = []
            warnings = []

            if not taxonomy_data:
                # 一行都没解析出来 = 这项检查实际上什么都没做（多列 TSV 表格、
                # 或分类学列不含 'd__X;p__Y' 式字符串时就是这个形态）
                warnings.append(
                    f"未从 {file_path} 解析出任何 'd__X;p__Y' 式分类学字符串，"
                    f"循环依赖检查本次**未覆盖任何行**（多列表格需把分类学写成 "
                    f"d__…;p__… 字符串形式，否则冲突无法被发现）"
                )

            if circular_deps:
                dep_descriptions = [f"'{p}' <-> '{c}'" for p, c in circular_deps]
                raise TaxonomyConflictError(
                    f"检测到 {len(circular_deps)} 个分类学循环依赖: {', '.join(dep_descriptions[:5])}",
                    suggestion="请检查分类学表格中是否存在名称在不同层级间循环引用的情况。",
                )

            return self.TaxonomyValidationDetail(
                is_valid=len(errors) == 0,
                circular_dependencies=circular_deps,
                errors=errors,
                warnings=warnings,
                entries_checked=len(taxonomy_data),
            )

        except TaxonomyConflictError:
            raise
        except Exception as e:
            return self.TaxonomyValidationDetail(
                is_valid=False, errors=[f"解析分类学文件失败: {e}"]
            )

    def _parse_taxonomy_lines(self, lines: List[str]) -> Dict[str, Dict[str, str]]:
        """
        解析分类学表格行

        Args:
            lines: 文件行列表

        Returns:
            {物种名: {级别: 值}} 字典
        """
        data = {}

        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # 检测分隔符
            if "\t" in line:
                parts = line.split("\t")
            elif "," in line:
                parts = line.split(",")
            else:
                continue

            if len(parts) < 2:
                continue

            name = parts[0].strip()
            taxonomy_str = parts[1].strip() if len(parts) > 1 else ""

            if not name or not taxonomy_str:
                continue

            # 解析分类学字符串
            taxonomy = self._parse_taxonomy_string(taxonomy_str)
            if taxonomy:
                data[name] = taxonomy

        return data

    def _parse_taxonomy_string(self, taxonomy_str: str) -> Dict[str, str]:
        """
        解析分类学字符串

        支持格式：
        - d__Bacteria;p__Proteobacteria;...
        - _d_Bacteria_p_Proteobacteria_...
        """
        result = {}

        # 双下划线格式
        pattern = re.compile(r"([a-z]{1,2})__([^;]+)")
        for match in pattern.finditer(taxonomy_str):
            prefix = match.group(1)
            value = match.group(2).strip()
            if value:
                result[prefix] = value

        # 如果没有匹配，尝试单下划线格式
        if not result:
            pattern = re.compile(r"_([a-z])_([^_]+)")
            for match in pattern.finditer(taxonomy_str):
                prefix = match.group(1)
                value = match.group(2).strip()
                if value:
                    result[prefix] = value

        return result

    def _detect_circular_dependencies(
        self, taxonomy_data: Dict[str, Dict[str, str]]
    ) -> List[Tuple[str, str]]:
        """
        检测分类学表格中的循环依赖

        例如：
        - 行1: d__A;p__B
        - 行2: d__B;p__A
        这构成了 A -> B -> A 的循环

        Args:
            taxonomy_data: 分类学数据

        Returns:
            循环依赖列表 [(parent, child), ...]
        """
        circular_deps = []

        # 构建父子关系图
        # 格式: {父级值: [子级值列表]}
        parent_child_map: Dict[str, Set[str]] = {}

        for name, taxonomy in taxonomy_data.items():
            levels = sorted(taxonomy.keys())
            for i in range(len(levels) - 1):
                parent_level = levels[i]
                child_level = levels[i + 1]
                parent_value = taxonomy.get(parent_level, "")
                child_value = taxonomy.get(child_level, "")

                if parent_value and child_value:
                    if parent_value not in parent_child_map:
                        parent_child_map[parent_value] = set()
                    parent_child_map[parent_value].add(child_value)

        # 检测循环：如果 A -> B 且 B -> A，则存在循环
        visited_pairs = set()
        for parent, children in parent_child_map.items():
            for child in children:
                # 检查是否存在反向关系
                if child in parent_child_map and parent in parent_child_map[child]:
                    # 避免重复报告
                    pair = tuple(sorted([parent, child]))
                    if pair not in visited_pairs:
                        visited_pairs.add(pair)
                        circular_deps.append((parent, child))

        return circular_deps
