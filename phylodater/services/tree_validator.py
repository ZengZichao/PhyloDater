"""
TreeValidator - 树结构验证服务

整合多个库的检测能力，形成最健壮的验证方案
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Tuple

from ..core.exceptions import AlignmentValidationError
from ..infrastructure.logging import get_logger
from ..models import PhylogeneticTree

if TYPE_CHECKING:  # 懒加载依赖，只在标注里出现
    from .deep_validator import DeepValidator


@dataclass
class ValidationReport:
    """验证报告"""

    is_valid: bool
    errors: List[str]
    warnings: List[str]


class TreeValidator:
    """树结构验证器"""

    def __init__(self) -> None:
        self.logger = get_logger()
        # 延迟到使用时创建，避免与 deep_validator 的构造顺序互相牵制
        self._deep_validator: Optional["DeepValidator"] = None

    @property
    def deep_validator(self) -> "DeepValidator":
        """:class:`~phylodater.services.deep_validator.DeepValidator` 懒加载实例。

        超度量性/时间树判定（B-18）复用它的实现，避免两个校验层各写一份
        root-to-tip 计算。
        """
        if self._deep_validator is None:
            from .deep_validator import DeepValidator

            self._deep_validator = DeepValidator()
        return self._deep_validator

    def validate_rooted(self, tree: PhylogeneticTree) -> bool:
        """
        验证树是否已定根

        ``PhylogeneticTree.is_rooted`` 是**属性**（返回 bool），旧实现写成
        ``tree.is_rooted()``，对真实树必然抛 ``TypeError: 'bool' object is not
        callable``；这里同时兼容"以方法暴露 is_rooted 的对象/测试替身"。
        """
        rooted = getattr(tree, "is_rooted", None)
        if rooted is None:
            info = tree.get_rooting_info()
            return bool(info.get("is_rooted", False))
        if callable(rooted):
            rooted = rooted()
        return bool(rooted)

    def validate_binary(self, tree: PhylogeneticTree) -> bool:
        """验证树是否为二叉树"""
        return tree.is_binary()

    def validate_ultrametric(
        self,
        tree: PhylogeneticTree,
        method: Optional[str] = None,
        tolerance: Optional[float] = None,
        require_time_tree: bool = False,
    ) -> ValidationReport:
        """输入树是否为**时间树**（超度量）——审阅项 B-18。

        整个包此前没有任何一处检验这件事，而 LSD2/IQ-TREE2 尖端定年、
        MCMCTree 全局钟（clock=1，分支长度被当作时间）以及可视化层的
        ``root_age - dist(root, node)`` 节点年龄公式都把它当前提。

        分档（按方法）：
        - treePL / PATHd8 / r8s / MD-Cat / wLogDate：**info**（这些方法设计上
          就允许非超度量输入，只告知不拦）；
        - LSD2 / IQ-TREE2 / MCMCTree / 可视化通路：**warning**（默认也不 hard-fail）；
        - ``require_time_tree=True``：升级为 **error**。

        Returns:
            ValidationReport: errors 只在 ``require_time_tree=True`` 时才非空
        """
        errors, warnings, infos = self.deep_validator.validate_ultrametricity(
            tree,
            method=method,
            tolerance=tolerance,
            require_time_tree=require_time_tree,
        )
        for note in infos:
            self.logger.info(note)
        for warning in warnings:
            self.logger.warning(warning)
        for error in errors:
            self.logger.error(error)
        return ValidationReport(is_valid=not errors, errors=errors, warnings=warnings)

    def validate_branch_lengths(
        self, tree: PhylogeneticTree
    ) -> Tuple[bool, List[str], List[str]]:
        """验证分支长度（返回: 是否有效, 错误列表, 警告列表）"""
        return tree.validate_branch_lengths()

    def validate_sequence_consistency(
        self, alignment_path: Path, tree: PhylogeneticTree
    ) -> ValidationReport:
        """
        验证树与比对文件的一致性

        检查项目：
        1. 物种数量是否一致
        2. 序列 ID 集合与树的 Tip Label 集合是否完全双向匹配（Set Equality）
        3. 如有差异，打印具体的差异 ID 列表

        支持格式：FASTA、PHYLIP、Stockholm、Clustal
        """
        errors = []
        warnings: List[str] = []

        # 获取比对文件的序列名称（支持多格式）
        try:
            seq_names = self._extract_alignment_sequence_names(alignment_path)
        except Exception as e:
            errors.append(f"Failed to parse alignment file: {e}")
            return ValidationReport(False, errors, warnings)

        # 获取树的叶节点名称
        tree_names = set()
        for name in tree.tip_names:
            tree_names.add(self._normalize_name(name))

        # 1. 检查物种数量一致性
        if len(tree_names) != len(seq_names):
            errors.append(
                f"Species count mismatch: tree has {len(tree_names)} tips, "
                f"alignment has {len(seq_names)} sequences"
            )

        # 2. 检查双向集合匹配（Set Equality）
        in_tree_not_in_aln = tree_names - seq_names
        in_aln_not_in_tree = seq_names - tree_names

        if in_tree_not_in_aln:
            missing_list = sorted(list(in_tree_not_in_aln))
            errors.append(f"{len(missing_list)} tip(s) in tree not found in alignment:")
            for name in missing_list:
                errors.append(f"  - {name}")
            suggestions = self._suggest_name_corrections(missing_list, seq_names)
            if suggestions:
                errors.append(f"  Possible name corrections: {suggestions}")

        if in_aln_not_in_tree:
            extra_list = sorted(list(in_aln_not_in_tree))
            errors.append(
                f"{len(extra_list)} sequence(s) in alignment not found in tree:"
            )
            for name in extra_list:
                errors.append(f"  - {name}")

        is_valid = len(errors) == 0
        return ValidationReport(is_valid, errors, warnings)

    def _extract_alignment_sequence_names(self, alignment_path: Path) -> set:
        """
        从比对文件中提取序列名称（支持 FASTA、PHYLIP、Stockholm、Clustal 格式）

        Returns:
            规范化后的序列名称集合
        """
        fmt = self._detect_alignment_format(alignment_path)

        if fmt == "phylip":
            return self._extract_phylip_names(alignment_path)

        from Bio import SeqIO

        seq_names = set()
        for record in SeqIO.parse(alignment_path, fmt):
            name = record.id.strip()
            name = self._normalize_name(name)
            seq_names.add(name)

        if not seq_names:
            raise AlignmentValidationError(
                f"No sequences found in alignment file: {alignment_path}"
            )

        return seq_names

    def _extract_phylip_names(self, file_path: Path) -> set:
        """
        从 PHYLIP 格式文件中提取序列名称

        PHYLIP 格式特殊处理：BioPython 的 PHYLIP 解析器可能将序列包含在 ID 中，
        因此手动解析以确保只提取名称部分。
        """
        seq_names = set()

        with open(file_path, "r") as f:
            lines = f.readlines()

        if not lines:
            raise AlignmentValidationError(f"Empty PHYLIP file: {file_path}")

        # 解析首行获取序列数量
        header = lines[0].strip()
        match = re.match(r"^\s*(\d+)\s+(\d+)", header)
        if not match:
            raise AlignmentValidationError(f"Invalid PHYLIP header: {header}")

        num_seqs = int(match.group(1))

        # 解析序列名称（跳过首行）
        for line in lines[1 : num_seqs + 1]:
            line = line.strip()
            if not line:
                continue
            # PHYLIP 格式：名称和序列用空格分隔，名称在前
            parts = line.split()
            if parts:
                name = parts[0]
                name = self._normalize_name(name)
                seq_names.add(name)

        if not seq_names:
            raise AlignmentValidationError(
                f"No sequences found in PHYLIP file: {file_path}"
            )

        return seq_names

    def _detect_alignment_format(self, file_path: Path) -> str:
        """
        检测比对文件格式

        Returns:
            BioPython 支持的格式名称 ('fasta', 'phylip', 'stockholm', 'clustal')
        """
        ext = file_path.suffix.lower()

        # 根据扩展名判断
        format_map = {
            ".fasta": "fasta",
            ".fa": "fasta",
            ".faa": "fasta",
            ".fna": "fasta",
            ".phylip": "phylip",
            ".phy": "phylip",
            ".stockholm": "stockholm",
            ".sto": "stockholm",
            ".clustal": "clustal",
            ".aln": "clustal",
        }

        if ext in format_map:
            return format_map[ext]

        # 根据首行内容判断
        with open(file_path, "r") as f:
            first_line = f.readline().strip()

        if first_line.startswith(">"):
            return "fasta"
        elif first_line.startswith("# STOCKHOLM"):
            return "stockholm"
        elif first_line.upper().startswith("CLUSTAL"):
            return "clustal"
        elif first_line and first_line[0].isdigit():
            return "phylip"

        # 默认使用 FASTA
        self.logger.warning(
            f"Could not detect alignment format for {file_path}, assuming FASTA"
        )
        return "fasta"

    def _suggest_name_corrections(
        self, missing_names: Iterable[str], available_names: Iterable[str]
    ) -> Dict[str, str]:
        """建议可能的名称修正（用于模糊匹配失败时的错误提示）"""
        suggestions = {}
        for missing in missing_names:
            missing_lower = missing.lower().replace("_", "").replace("-", "")
            for available in available_names:
                available_lower = available.lower().replace("_", "").replace("-", "")
                if missing_lower == available_lower:
                    suggestions[missing] = available
                    break
                if self._levenshtein_distance(missing_lower, available_lower) <= 2:
                    suggestions[missing] = available
                    break
        return suggestions

    def _levenshtein_distance(self, s1: str, s2: str) -> int:
        """计算两个字符串之间的 Levenshtein 距离"""
        if len(s1) < len(s2):
            return self._levenshtein_distance(s2, s1)
        if len(s2) == 0:
            return len(s1)

        # 首行是 "编辑距离 = 列号" 的序列，后续行是 list；两者都必须能按下标
        # 读写，所以把首行就算成 list（range 不是一种 list）。
        previous_row: List[int] = list(range(len(s2) + 1))
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row

        return previous_row[-1]

    def _normalize_name(self, name: str) -> str:
        """规范化名称用于比较"""
        # 替换特殊字符
        normalized = re.sub(r"[;:,()\[\]]", "_", name)
        return normalized.lower().strip()

    def full_validation(
        self,
        tree: PhylogeneticTree,
        alignment_path: Optional[Path] = None,
        require_rooted: bool = False,
        method: Optional[str] = None,
        require_time_tree: bool = False,
    ) -> ValidationReport:
        """
        执行完整的验证流程

        Args:
            tree: 系统发育树
            alignment_path: 比对文件路径（可选）
            require_rooted: 是否强制要求有根树（某些方法如 LSD2 不要求）
            method: 本次将要使用的定年方法名（用于超度量性检查的分档，B-18）。
                为 None 时按"假定时间树"的严格档给 warning。
            require_time_tree: 把"输入树不是时间树"从 warning 升级为 error

        Returns:
            验证报告
        """
        errors = []
        warnings = []

        rooting_info = tree.get_rooting_info()

        if not rooting_info["is_rooted"]:
            if require_rooted:
                errors.append(
                    "Tree is not rooted. Please provide a rooted tree or use --rooting option."
                )
            else:
                warnings.append(
                    "Tree is not rooted. Some dating methods may require a rooted tree. "
                    "Consider using --rooting option."
                )
                if rooting_info.get("likely_artifact"):
                    # 使用 .get 守卫，避免 get_rooting_info 未提供这些键时 KeyError
                    warnings.append(
                        "Rooting artifact detected: tree may have been rerooted previously "
                        "(short branch ratio: {shortest_ratio:.3f})".format(
                            shortest_ratio=rooting_info.get("shortest_ratio", 0.0)
                        )
                    )

        # 2. 验证二叉性
        if not self.validate_binary(tree):
            warnings.append("Tree is not strictly binary (contains polytomies).")

        # 3. 验证分支长度
        # validate_branch_lengths 返回 (is_valid, invalid_nodes, errors)，
        # 此处曾误将返回值解包为 (is_valid, errors, warnings)，导致负值/零分支长度
        # 被当作 warnings 而非 errors。修正为正确的三元组解包。
        lengths_valid, invalid_nodes, length_errors = self.validate_branch_lengths(tree)
        if not lengths_valid:
            errors.extend(length_errors)

        # 4. 超度量性 / 时间树判定（B-18）——默认只 warn，不拦下运行
        ultrametric = self.validate_ultrametric(
            tree, method=method, require_time_tree=require_time_tree
        )
        errors.extend(ultrametric.errors)
        warnings.extend(ultrametric.warnings)

        # 5. 验证与比对文件的一致性
        if alignment_path:
            consistency = self.validate_sequence_consistency(alignment_path, tree)
            errors.extend(consistency.errors)
            warnings.extend(consistency.warnings)

        is_valid = len(errors) == 0

        if is_valid:
            self.logger.success("Tree validation passed")
        else:
            self.logger.error(f"Tree validation failed: {errors}")

        return ValidationReport(is_valid, errors, warnings)
