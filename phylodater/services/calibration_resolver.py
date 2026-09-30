"""
CalibrationResolver - MRCA 节点定位服务

五层递进解析策略，最大化系统发育多样性
"""

import math
import re
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

from ..core.exceptions import (
    CalibrationConflictError,
    CalibrationResolutionError,
    MonophylyError,
)
from ..infrastructure.logging import get_logger
from ..models import CalibrationPoint, PhylogeneticTree
from ..models.constraints import AgeConstraint, FixedAgeConstraint, UniformAgeConstraint

if TYPE_CHECKING:  # 只供标注；taxonomy_parser 是可选注入
    from .taxonomy_parser import TaxonomyParser

# 特殊标识符集合（大小写敏感，仅接受大写形式）
SPECIAL_IDENTIFIERS = frozenset({"LUCA", "LBCA", "LACA", "ROOT"})


@dataclass
class ResolutionResult:
    """解析结果"""

    #: 纯节点定位（按叶名/类群找 MRCA）时并不需要一个完整的校准点，
    #: 所以这里按事实是 Optional（旧标注说它必有，与自身的 ``None`` 传入直接矛盾）。
    calibration_point: Optional[CalibrationPoint]
    resolved_taxa: List[str]
    mrca_leaf_pair: Tuple[str, str]
    is_root_node: bool = False


class CalibrationResolver:
    """
    MRCA 节点定位服务

    五层递进解析策略：
    1. 精确名称匹配
    2. 分类学信息匹配
    3. 词边界匹配
    4. 精确分词匹配
    5. 系统发育多样性最大化
    """

    def __init__(
        self,
        tree: PhylogeneticTree,
        taxonomy_parser: Optional["TaxonomyParser"] = None,
    ) -> None:
        """
        初始化解析器

        Args:
            tree: 系统发育树
            taxonomy_parser: 分类学解析器（可选）

        Raises:
            TypeError: 如果 tree 为 None
            ValueError: 如果 tree 没有叶节点
        """
        # 空值验证
        if tree is None:
            raise TypeError("tree cannot be None")

        self.tree = tree
        self.tip_names = tree.tip_names or []

        # 验证树有叶节点
        if not self.tip_names:
            raise ValueError("Tree has no tips (leaf nodes)")

        self.logger = get_logger()
        self.taxonomy_parser = taxonomy_parser

        # 预计算并缓存距离矩阵（用于第五层多样性最大化）。
        # get_distance_matrix 自身有实例级缓存，这里预热一次即可避免
        # 每个解析目标在第五层重复计算 O(n³) 距离矩阵（此前该字段声明后从未使用）。
        try:
            self._distance_matrix = self.tree.get_distance_matrix()
        except Exception as e:
            self.logger.warning(f"Failed to precompute distance matrix: {e}")
            self._distance_matrix = {}

    def resolve(self, target: str) -> ResolutionResult:
        """
        解析校准点目标

        Args:
            target: 目标名称（如 LUCA、LACA、LBCA、Cyanobacteriota）

        Returns:
            ResolutionResult: 解析结果

        Raises:
            TypeError: 如果 target 为 None
            ValueError: 如果 target 为空字符串
            CalibrationResolutionError: 如果无法解析目标或目标非单系群
        """
        # 空值验证
        if target is None:
            raise TypeError("target cannot be None")
        if not isinstance(target, str):
            raise TypeError(f"target must be a string, got {type(target).__name__}")
        if not target.strip():
            raise ValueError("target cannot be empty or whitespace only")

        self.logger.info(f"Resolving calibration target: {target}")

        # 第零层：特殊标识符（LUCA、LACA、LBCA）
        special_resolved = self._resolve_special_identifier(target)
        if special_resolved is not None:
            self.logger.success(
                f"Special identifier '{target}' resolved: {len(special_resolved)} tips"
            )
            # 对 LBCA/LACA（Bacteria/Archaea 作为独立 clade）进行单系性检查：
            # 若 Bacteria/Archaea 在树中并非单系（如根落在 Bacteria 内、或含非细菌
            # 外群），定位到的 MRCA 会错误地包含非目标域的叶节点，必须报错而非
            # 静默定位到错误节点（甚至根）。LUCA 的 MRCA 天然包含真核生物、
            # ROOT 为全树根，不做此严格检查。
            target_upper = target.upper()
            if target_upper in ("LBCA", "LACA"):
                self._check_clade_monophyly(target, special_resolved)
            return self._create_result(target, special_resolved, is_special=True)

        # 第一层：精确名称匹配
        resolved = self._exact_match(target)
        if resolved:
            self.logger.success(f"Exact match found: {len(resolved)} tips")
            return self._resolve_with_monophyly_check(target, resolved)

        # 第二层：分类学信息匹配
        if self.taxonomy_parser:
            resolved = self._taxonomy_match(target)
            if resolved:
                self.logger.success(f"Taxonomy match found: {len(resolved)} tips")
                return self._resolve_with_monophyly_check(target, resolved)

        # 第三层：词边界匹配
        resolved = self._word_boundary_match(target)
        if resolved:
            self.logger.success(f"Word boundary match found: {len(resolved)} tips")
            return self._resolve_with_monophyly_check(target, resolved)

        # 第四层：精确分词匹配
        resolved = self._token_match(target)
        if resolved:
            self.logger.success(f"Token match found: {len(resolved)} tips")
            return self._resolve_with_monophyly_check(target, resolved)

        # 第五层：系统发育多样性最大化
        resolved = self._diversity_maximization(target)
        if resolved:
            self.logger.success(f"Diversity maximization found: {len(resolved)} tips")
            return self._resolve_with_monophyly_check(target, resolved)

        raise CalibrationResolutionError(
            f"Could not resolve calibration target: {target}"
        )

    def resolve_mrca_pair(self, mrca_pair: Tuple[str, str]) -> ResolutionResult:
        """
        通过 mrca_pair 解析校准点

        Args:
            mrca_pair: 两个叶节点名称的元组

        Returns:
            ResolutionResult: 解析结果

        Raises:
            TypeError: 如果 mrca_pair 为 None
            ValueError: 如果 mrca_pair 格式不正确
            CalibrationResolutionError: 如果无法解析
        """
        if mrca_pair is None:
            raise TypeError("mrca_pair cannot be None")
        if not isinstance(mrca_pair, (tuple, list)) or len(mrca_pair) != 2:
            raise ValueError(
                f"mrca_pair must be a tuple/list of 2 elements, got {mrca_pair}"
            )

        leaf1, leaf2 = mrca_pair

        # 验证叶节点存在
        if leaf1 not in self.tip_names:
            raise CalibrationResolutionError(f"Leaf '{leaf1}' not found in tree")
        if leaf2 not in self.tip_names:
            raise CalibrationResolutionError(f"Leaf '{leaf2}' not found in tree")

        # 获取 MRCA 的所有后代叶节点
        mrca_terminals = self.tree.get_mrca_terminals_all([leaf1, leaf2])
        if not mrca_terminals:
            raise CalibrationResolutionError(
                f"Could not find MRCA for {leaf1} and {leaf2}"
            )

        # 取系统发育距离最远的两个叶节点作为 mrca_leaf_pair，提高定位稳定性
        pair1, pair2 = self.tree.get_mrca_terminals([leaf1, leaf2])
        if pair1 is None or pair2 is None:
            # 空叶名对写不出 MRCA 定位；这里必须报错，而不是把
            # (None, None) 当成一对可用叶名传下去。
            raise CalibrationResolutionError(
                f"Could not pick a representative leaf pair for {leaf1} and {leaf2}"
            )

        self.logger.success(
            f"Resolved MRCA pair ({leaf1}, {leaf2}): {len(mrca_terminals)} tips"
        )

        return ResolutionResult(
            calibration_point=None,
            resolved_taxa=mrca_terminals,
            mrca_leaf_pair=(pair1, pair2),
            is_root_node=False,
        )

    def _resolve_special_identifier(self, target: str) -> Optional[List[str]]:
        """
        解析特殊标识符（LUCA、LACA、LBCA）

        LUCA: 所有 Bacteria 和 Archaea 的 MRCA（返回各一个代表叶节点）
        LBCA: 所有 Bacteria 的 MRCA（返回首尾两个叶节点）
        LACA: 所有 Archaea 的 MRCA（返回首尾两个叶节点）

        Args:
            target: 目标名称

        Returns:
            解析到的叶节点列表，非特殊标识符返回 None
        """
        target_upper = target.upper()
        if target_upper not in SPECIAL_IDENTIFIERS:
            return None

        # 特殊标识符大小写敏感：仅接受大写，小写形式直接拒绝
        if target != target_upper:
            raise CalibrationResolutionError(
                f"特殊标识符大小写敏感：'{target}' 应为大写 '{target_upper}'"
            )

        # ROOT：全树根节点（返回根节点下所有叶节点，使 MRCA 定位到根）。
        # 走 PhylogeneticTree 的公开接口 tip_names —— 根节点的后代恒等于全树叶
        # 集合，因此不需要任何树解析后端（ete3/dendropy 皆不必）。
        if target_upper == "ROOT":
            if not self.tip_names:
                raise CalibrationResolutionError("Tree has no root node")
            return list(self.tip_names)

        bacteria = []
        archaea = []

        for name in self.tip_names:
            taxonomy = None
            if self.taxonomy_parser:
                taxonomy = self.taxonomy_parser.parse(name)

            if taxonomy:
                domain = (taxonomy.get("domain") or "").lower()
                if domain == "bacteria":
                    bacteria.append(name)
                elif domain == "archaea":
                    archaea.append(name)
            else:
                name_lower = name.lower()
                if "_d_bacteria" in name_lower:
                    bacteria.append(name)
                elif "_d_archaea" in name_lower:
                    archaea.append(name)
                elif "bacteria" in name_lower:
                    # 注意：Methanobacteriota 是细菌门，必须包含在 Bacteria 中，
                    # 此前错误地以 "methanobacteria" 子串排除会导致甲烷细菌被漏判。
                    bacteria.append(name)
                elif "archaea" in name_lower:
                    archaea.append(name)

        if target_upper == "LUCA":
            if not bacteria or not archaea:
                raise CalibrationResolutionError(
                    f"Cannot resolve LUCA: found {len(bacteria)} Bacteria and "
                    f"{len(archaea)} Archaea. Need at least one of each."
                )
            return list(bacteria) + list(archaea)

        elif target_upper == "LBCA":
            if len(bacteria) < 2:
                raise CalibrationResolutionError(
                    f"Cannot resolve LBCA: found only {len(bacteria)} Bacteria tips. "
                    f"Need at least 2."
                )
            return list(bacteria)

        elif target_upper == "LACA":
            if len(archaea) < 2:
                raise CalibrationResolutionError(
                    f"Cannot resolve LACA: found only {len(archaea)} Archaea tips. "
                    f"Need at least 2."
                )
            return list(archaea)

        return None

    def _check_clade_monophyly(self, target: str, resolved: List[str]) -> None:
        """检查已解析的特殊类群（Bacteria/Archaea）是否为单系 clade。

        Args:
            target: 标识符（LBCA/LACA）
            resolved: 解析到的叶节点列表

        Raises:
            MonophylyError: 若该类群不是单系群（其 MRCA 下包含非该类群的叶节点）
        """
        mrca_terminals = self.tree.get_mrca_terminals_all(resolved)
        mrca_set = set(mrca_terminals)
        resolved_set = set(resolved)
        if mrca_set != resolved_set:
            extra = sorted(mrca_set - resolved_set)
            extra_sample = ", ".join(extra[:5])
            raise MonophylyError(
                f"'{target}' 不是单系群（non-monophyletic）。解析到 {len(resolved)} 个叶节点，"
                f"但其 MRCA 下共有 {len(mrca_set)} 个叶节点，"
                f"含 {len(extra)} 个不属于 '{target}' 的叶节点"
                f"（如: {extra_sample}{'...' if len(extra) > 5 else ''}）。"
                f"请检查树中 {target} 的组成与定根方式"
                f"（例如根是否落在 Bacteria 内或存在非细菌外群）。"
            )

    def _resolve_with_monophyly_check(
        self, target: str, resolved: List[str]
    ) -> ResolutionResult:
        """
        解析目标并检查单系群

        Args:
            target: 目标名称
            resolved: 解析到的叶节点列表

        Returns:
            ResolutionResult

        Raises:
            CalibrationResolutionError: 如果解析到的类群不构成单系群
        """
        if len(resolved) >= 2:
            # 获取 MRCA 下的所有叶节点
            mrca_terminals = self.tree.get_mrca_terminals_all(resolved)
            mrca_set = set(mrca_terminals)
            resolved_set = set(resolved)

            # 单系群判定：MRCA 下的所有叶节点必须恰好等于解析到的叶节点集合
            if mrca_set != resolved_set:
                extra = mrca_set - resolved_set
                missing = resolved_set - mrca_set

                detail_parts = []
                if extra:
                    extra_sample = sorted(extra)[:5]
                    detail_parts.append(
                        f"MRCA 下还有 {len(extra)} 个不属于 '{target}' 的叶节点"
                        f"（如: {', '.join(extra_sample)}{'...' if len(extra) > 5 else ''}）"
                    )
                if missing:
                    missing_sample = sorted(missing)[:5]
                    detail_parts.append(
                        f"有 {len(missing)} 个属于 '{target}' 的叶节点不在 MRCA 下"
                        f"（如: {', '.join(missing_sample)}{'...' if len(missing) > 5 else ''}）"
                    )

                detail_msg = "; ".join(detail_parts)
                raise MonophylyError(
                    f"'{target}' 不是单系群（non-monophyletic），无法通过指定类群名称进行操作。"
                    f"解析到 {len(resolved)} 个叶节点，但 MRCA 下共有 {len(mrca_set)} 个叶节点。"
                    f"{detail_msg}"
                )

            self.logger.info(
                f"Monophyly check passed for '{target}': {len(resolved)} tips"
            )

        return self._create_result(target, resolved)

    def _exact_match(self, target: str) -> List[str]:
        """精确名称匹配"""
        return [name for name in self.tip_names if name == target]

    def _taxonomy_match(self, target: str) -> List[str]:
        """分类学信息匹配"""
        if not self.taxonomy_parser:
            return []

        target_lower = target.lower()
        resolved = []

        for name in self.tip_names:
            taxonomy = self.taxonomy_parser.parse(name)
            if taxonomy:
                # 检查所有分类级别
                for rank, value in taxonomy.items():
                    # TaxonomyMap 的值可以是 None（缺失等级显式入库），先判空。
                    if value and value.lower() == target_lower:
                        resolved.append(name)
                        break

        return resolved

    def _word_boundary_match(self, target: str) -> List[str]:
        """词边界匹配"""
        pattern = rf"(_|^){re.escape(target)}(_|$)"
        return [
            name for name in self.tip_names if re.search(pattern, name, re.IGNORECASE)
        ]

    def _token_match(self, target: str) -> List[str]:
        """精确分词匹配"""
        target_lower = target.lower()

        resolved = []
        for name in self.tip_names:
            # 规范化名称
            normalized = re.sub(r"[;_.\-\s]+", " ", name.lower())
            name_tokens = set(normalized.split())

            # 单词目标：集合成员检测
            if " " not in target_lower:
                if target_lower in name_tokens:
                    resolved.append(name)
            else:
                # 多词目标：短语匹配
                pattern = rf"\b{re.escape(target_lower)}\b"
                if re.search(pattern, normalized):
                    resolved.append(name)

        return resolved

    def _diversity_maximization(self, target: str) -> List[str]:
        """
        系统发育多样性最大化

        选择覆盖面最广的代表性叶节点
        """
        # 获取所有候选节点
        candidates = self._token_match(target)
        if not candidates:
            # 如果分词匹配失败，尝试模糊匹配
            target_lower = target.lower()
            candidates = [
                name for name in self.tip_names if target_lower in name.lower()
            ]

        if len(candidates) <= 1:
            return candidates

        if len(candidates) == 2:
            return candidates

        # 使用预计算并缓存的距离矩阵选择最远端的节点
        try:
            distances = self._distance_matrix

            # 选择距离最远的两个节点作为代表
            # 距离是 float（矩阵里存的就是 float），首值也要按 float 起步，
            # 否则后面的 ``max_dist = dist`` 会被判成"float 赋给 int"。
            max_dist: float = 0.0
            best_pair = (candidates[0], candidates[1])

            for i, tip1 in enumerate(candidates):
                for tip2 in candidates[i + 1 :]:
                    dist = distances.get((tip1, tip2), distances.get((tip2, tip1), 0))
                    if dist > max_dist:
                        max_dist = dist
                        best_pair = (tip1, tip2)

            return list(best_pair)

        except Exception as e:
            self.logger.warning(f"Distance matrix calculation failed: {e}")
            # 回退到首尾选择
            return [candidates[0], candidates[-1]]

    def _create_result(
        self,
        target: str,
        resolved: List[str],
        age_constraint: Optional[AgeConstraint] = None,
        is_special: bool = False,
    ) -> ResolutionResult:
        """创建解析结果

        Args:
            target: 目标名称
            resolved: 解析到的叶节点列表
            age_constraint: 年龄约束（可选，由调用方在后续从原始校准点复制）
            is_special: 是否为特殊标识符（LUCA/LACA/LBCA），跳过单系群检查

        Returns:
            ResolutionResult

        Note:
            age_constraint 默认为 None。调用方（如 CLI 和 Pipeline 的校准解析流程）
            负责从原始 CalibrationPoint 复制 age_constraint 到解析后的校准点上。
            此方法返回的 ResolutionResult.calibration_point 仅包含节点定位信息，
            不包含年龄约束。
        """
        # 确定是否为根节点
        is_root = target.upper() in ["LUCA", "ROOT"]

        # 获取 MRCA 叶节点对：使用系统发育距离最远的一对，确保 PATHd8/treePL/r8s 定位准确
        if len(resolved) >= 2:
            try:
                mrca_pair = self.tree.get_mrca_terminals(resolved)
            except Exception:
                self.logger.warning(
                    f"Failed to compute phylogenetically most-distant pair for '{target}', "
                    f"falling back to first/last taxa"
                )
                mrca_pair = (resolved[0], resolved[-1])
        elif len(resolved) == 1:
            raise CalibrationResolutionError(
                f"Cannot define MRCA for '{target}': only one taxon ({resolved[0]}) resolved. "
                f"Need at least two taxa to define a node."
            )
        else:
            raise CalibrationResolutionError(
                f"Cannot define MRCA for '{target}': no taxa resolved."
            )

        # get_mrca_terminals 只在入参为空时返回 (None, None)，而上面的分支已经保证
        # resolved 至少有两个名字；这里仍然把那个前提核住，免得把 None 当叶名传下去。
        pair1, pair2 = mrca_pair
        if pair1 is None or pair2 is None:
            raise CalibrationResolutionError(
                f"Could not pick a representative leaf pair for '{target}' "
                f"from {resolved[:2]}"
            )
        mrca_pair = (pair1, pair2)

        # 创建 CalibrationPoint（age_constraint 由调用方从原始校准点复制）
        # 抑制 age_constraint=None 的 UserWarning，因为约束由调用方后续设置
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message=".*age_constraint.*", category=UserWarning
            )
            calibration_point = CalibrationPoint(
                name=target,
                age_constraint=age_constraint,
                resolved_taxa=resolved,
                mrca_leaf_pair=mrca_pair,
                is_root_node=is_root,
            )

        return ResolutionResult(
            calibration_point=calibration_point,
            resolved_taxa=resolved,
            mrca_leaf_pair=mrca_pair,
            is_root_node=is_root,
        )

    def validate_ancestry(self, calibration_points: List[CalibrationPoint]) -> bool:
        """
        验证祖先-后代时序约束（"祖先必须老于后代"这一不变量的守卫）

        检查校准点之间是否存在时序冲突。

        审阅项 B-7：每一对校准点都以**两个方向**分别比较（A 是否为 B 的祖先、
        B 是否为 A 的祖先），因此是否漏检不再依赖校准点在 YAML 里的书写顺序。
        审阅项 B-8：祖先关系由 **MRCA 的完整后代叶集合** 决定，而不是用户为求
        MRCA 而指定的代表类群（``resolved_taxa``）——"鸟+鳄鱼"定 Archosauria、
        "三角龙+霸王龙"定其内部节点这类**代表类群互不重叠**的嵌套写法此前会被
        完全放过。
        审阅项 B-8 关联问题：两个校准点定位到**同一节点**（后代叶集合完全相同）
        时同样必须比较，否则"同一节点给了两条不相容的约束"无人拦截。

        Raises:
            CalibrationConflictError: 当检测到时序冲突时
        """
        self.logger.info("Validating ancestry constraints...")

        points = list(calibration_points or [])
        # 每个校准点所在节点的完整后代叶集合（每点只求一次 MRCA，避免 O(n²) 重复查询）
        clades: List[Optional[frozenset]] = [
            self._clade_leaf_set(cal) for cal in points
        ]

        conflicts: List[Tuple] = []
        total_pairs = 0
        unverifiable_pairs = 0

        for i, cal1 in enumerate(points):
            for j in range(i + 1, len(points)):
                total_pairs += 1
                if clades[i] is None or clades[j] is None:
                    unverifiable_pairs += 1
                # B-7：两个方向都查。同一对最多只有一个方向能给出"严格祖先"关系，
                # 因此不会重复报告；区间完全不相交时报告的总是"上界更小的一方"。
                for ancestor_idx, descendant_idx in ((i, j), (j, i)):
                    conflict_info = self._get_pair_conflict_details(
                        points[ancestor_idx],
                        clades[ancestor_idx],
                        points[descendant_idx],
                        clades[descendant_idx],
                    )
                    if conflict_info:
                        conflicts.append(conflict_info)

        if conflicts:
            error_lines = []
            for conflict in conflicts:
                (
                    ancestor,
                    descendant,
                    anc_constraint,
                    desc_constraint,
                    anc_bound,
                    desc_bound,
                    same_node,
                ) = conflict
                anc_str = self._format_constraint_for_error(
                    ancestor.name, anc_constraint, (0.0, anc_bound)
                )
                desc_str = self._format_constraint_for_error(
                    descendant.name, desc_constraint, (desc_bound, float("inf"))
                )
                relation = "is the SAME node as" if same_node else "is ancestral to"
                error_lines.append(
                    f'  - Node "{ancestor.name}" ({anc_str}) {relation} '
                    f'"{descendant.name}" ({desc_str}): the older bound '
                    f"{desc_bound} Ma exceeds the younger bound {anc_bound} Ma"
                )

            error_msg = (
                "Temporal conflicts detected between calibration points "
                "(an ancestor node must be at least as old as its descendants):\n"
                + "\n".join(error_lines)
            )
            raise CalibrationConflictError(error_msg)

        if total_pairs and unverifiable_pairs:
            # 无法判定拓扑时守卫等于没生效——必须说话，不能静默放行（C-8/C-31 同族）
            self.logger.error(
                f"Ancestry validation INCOMPLETE: {unverifiable_pairs}/{total_pairs} "
                f"calibration pairs could not be compared because their node "
                f"topology (MRCA leaf set) could not be resolved. The "
                f"'ancestor older than descendant' invariant is NOT verified for "
                f"those pairs; check that every calibration's taxa/mrca_pair "
                f"exists in the tree and that a tree backend is available."
            )

        self.logger.success("Ancestry validation passed")
        return True

    def _clade_leaf_set(self, cal: CalibrationPoint) -> Optional[frozenset]:
        """校准点所锚定节点的**完整**后代叶集合（拓扑判定 B-8 的唯一依据）。

        Returns:
            frozenset: 该节点下的全部叶节点名（根校准 = 全树叶集合）
            None: 无法判定（既非根、又无代表类群、或树后端求不出 MRCA）
        """
        if cal.is_root_node:
            tips = getattr(self.tree, "tip_names", None) or []
            return frozenset(tips) if tips else None

        taxa = [t for t in (cal.resolved_taxa or [])]
        if not taxa:
            self.logger.warning(
                f"Ancestry check: calibration '{cal.name}' has no resolved taxa; "
                f"its node cannot be located in the tree."
            )
            return None

        try:
            terminals = self.tree.get_mrca(taxa)
            if terminals is None:
                self.logger.warning(
                    f"Ancestry check: could not locate the MRCA of calibration "
                    f"'{cal.name}' (taxa {sorted(set(taxa))[:5]}...); treating its "
                    f"node position as undecidable rather than guessing."
                )
                return None
            return frozenset(terminals)
        except Exception as e:  # 后端不可用（如树解析库缺失）
            self.logger.warning(
                f"Ancestry check: MRCA lookup for '{cal.name}' failed ({e}); "
                f"treating its node position as undecidable."
            )
            return None

    def _get_pair_conflict_details(
        self,
        ancestor: CalibrationPoint,
        ancestor_clade: Optional[frozenset],
        descendant: CalibrationPoint,
        descendant_clade: Optional[frozenset],
    ) -> Optional[Tuple]:
        """在拓扑关系已知的前提下检查一对校准点的时序冲突。

        冲突判据不变：``max_age(祖先) < min_age(后代)``（区间可行性条件）。
        ``ancestor_clade == descendant_clade``（同一节点）时同样比较——两条约束
        必须相交，否则用户给同一个节点写了彼此不相容的校准。

        Returns:
            Tuple: (ancestor, descendant, anc_constraint, desc_constraint,
            max_age_ancestor, min_age_descendant, same_node) 或 None
        """
        if ancestor_clade is None or descendant_clade is None:
            return None

        same_node = ancestor_clade == descendant_clade
        is_ancestor_relation = ancestor_clade > descendant_clade
        if not (same_node or is_ancestor_relation):
            return None

        anc_constraint = ancestor.age_constraint
        desc_constraint = descendant.age_constraint
        if anc_constraint is None or desc_constraint is None:
            return None

        max_age_anc = self._get_max_age(anc_constraint)
        min_age_desc = self._get_min_age(desc_constraint)
        if max_age_anc is None or min_age_desc is None:
            return None
        if not (math.isfinite(max_age_anc) and math.isfinite(min_age_desc)):
            # inf 上界（软下界）与 NaN 都不构成可比较的界；NaN 已在约束层拒绝，
            # 这里只是防止万一有绕过者把比较退化为恒假。
            if math.isnan(max_age_anc) or math.isnan(min_age_desc):
                self.logger.warning(
                    f"Ancestry check: non-finite age bound for '{ancestor.name}' / "
                    f"'{descendant.name}'; the pair cannot be compared."
                )
                return None

        if max_age_anc < min_age_desc:
            return (
                ancestor,
                descendant,
                anc_constraint,
                desc_constraint,
                max_age_anc,
                min_age_desc,
                same_node,
            )
        return None

    def _format_constraint_for_error(
        self,
        node_name: str,
        constraint: Optional[AgeConstraint],
        bound: Tuple[float, float],
    ) -> str:
        """格式化约束信息用于错误提示"""
        from ..models.constraints import (
            GammaPriorConstraint,
            MaximumAgeConstraint,
            SkewNormalConstraint,
            SkewTConstraint,
            SoftBoundsConstraint,
            SoftLowerBoundConstraint,
        )

        if isinstance(constraint, FixedAgeConstraint):
            return f"fixed_age={constraint.fixed_age}Ma"
        elif isinstance(constraint, UniformAgeConstraint):
            return f"range=[{constraint.min_age}-{constraint.max_age}]Ma"
        elif isinstance(constraint, SoftLowerBoundConstraint):
            return f"soft_lower_bound>={constraint.min_age}Ma"
        elif isinstance(constraint, MaximumAgeConstraint):
            return f"upper_bound<={constraint.max_age}Ma"
        elif isinstance(constraint, SoftBoundsConstraint):
            return f"soft_bounds=[{constraint.min_age}-{constraint.max_age}]Ma"
        elif isinstance(constraint, GammaPriorConstraint):
            # 错误提示应反映约束本身的先验统计量，而非借用祖先/后代的冲突边界
            # （bound 为 (max_age1, min_age2)，与本约束无关）。
            # 分位数统一走 _get_min_age/_get_max_age，避免与冲突判据使用两套数值
            # （也避免 scipy 缺失时在此处抛错）。
            _ci_lower = self._get_min_age(constraint)
            _ci_upper = self._get_max_age(constraint)
            _mean = (
                constraint.alpha * constraint.scale / constraint.beta
                + constraint.offset
            )
            return (
                f"gamma_prior(mean={_mean:.1f}Ma, "
                f"95%CI=[{_ci_lower:.1f}-{_ci_upper:.1f}]Ma)"
            )
        elif isinstance(constraint, (SkewNormalConstraint, SkewTConstraint)):
            kind = "normal" if isinstance(constraint, SkewNormalConstraint) else "t"
            _ci_lower = self._get_min_age(constraint)
            _ci_upper = self._get_max_age(constraint)
            if _ci_lower is None or _ci_upper is None:
                return f"skew_{kind}_prior(unquantifiable)"
            return (
                f"skew_{kind}_prior(95%CI=[{_ci_lower:.1f}-{_ci_upper:.1f}]Ma, "
                f"normal approximation)"
            )
        return "unknown_constraint"

    def _get_temporal_conflict_details(
        self, cal1: CalibrationPoint, cal2: CalibrationPoint
    ) -> Optional[Tuple]:
        """获取时序冲突详情（假定 cal1 是 cal2 的祖先）

        保留给外部/历史调用；``validate_ancestry`` 自身走
        :meth:`_get_pair_conflict_details`（一次 MRCA 查询、双向对称比较）。

        Returns:
            Tuple: (cal1, cal2, constraint1, constraint2, max_age1, min_age2) or None
        """
        conflict = self._get_pair_conflict_details(
            cal1, self._clade_leaf_set(cal1), cal2, self._clade_leaf_set(cal2)
        )
        if conflict is None:
            return None
        # 对旧调用方保持 6 元组形状（不带 same_node 标记）
        return conflict[:6]

    def _is_ancestor(self, cal1: CalibrationPoint, cal2: CalibrationPoint) -> bool:
        """检查 cal1 是否是 cal2 的祖先

        审阅项 B-8：判定**只**依据两个节点在树中的拓扑包含关系，即
        ``MRCA(cal1 的代表类群)`` 的完整后代叶集合是否严格包含
        ``MRCA(cal2 的代表类群)`` 的完整后代叶集合。

        此前这里有一道 ``set(cal1.resolved_taxa).issuperset(set(cal2.resolved_taxa))``
        的"类群超集"预过滤，而 ``resolved_taxa`` 是用户为求 MRCA 指定的**代表类群**
        （不是该 MRCA 的全部后代），于是"鸟+鳄鱼"与"两种恐龙"这类真实嵌套写法在
        预过滤处就被判为无关、拓扑判定永远拿不到本该由它裁决的输入（系统性假阴性）。

        对于根节点（``is_root_node=True``），其后代叶集合即全树，天然包含一切。
        无法定位 MRCA 时返回 False（"无法判定"不得声称"是祖先"），由
        ``validate_ancestry`` 汇总成显式的 error 日志。
        """
        clade1 = self._clade_leaf_set(cal1)
        clade2 = self._clade_leaf_set(cal2)
        if clade1 is None or clade2 is None:
            return False
        # 同一节点（叶集合相同）不算"祖先"，由 validate_ancestry 的同节点分支处理
        return clade1 > clade2

    # 95% 等尾区间的正态近似界：mean ± 1.96·sd（Z_{0.975}）
    _NORMAL_975_Z = 1.959963984540054

    def _normal_approximation_bound(
        self,
        constraint: AgeConstraint,
        mean: float,
        sd: float,
        tail: str,
    ) -> float:
        """矩估计 + 正态近似分位数（审阅项 C-3 的显式替代方案）。

        只在精确分位数确实算不出来时调用（scipy 不可用，或分布本身不存在），
        并且**必须**记一条 warning：近似区间比真实的偏态/重尾区间窄，用它做
        祖先-后代可行性判断是偏乐观的，不能伪装成精确值。
        """
        z = self._NORMAL_975_Z if tail == "upper" else -self._NORMAL_975_Z
        bound = mean + z * sd
        self.logger.warning(
            f"Ancestry check: the 2.5%/97.5% quantiles of "
            f"{type(constraint).__name__} are approximated by the symmetric normal "
            f"interval mean ± 1.96·sd (here {tail} bound = {bound:.3f} Ma) because "
            f"the exact quantile could not be evaluated; the true prior is skewed, "
            f"so this bound is optimistic (the conflict test can still miss a "
            f"violation)."
        )
        return bound

    def _get_min_age(self, constraint: AgeConstraint) -> Optional[float]:
        """获取约束的最小年龄（2.5% 分位数）"""
        from ..models.constraints import (
            GammaPriorConstraint,
            MaximumAgeConstraint,
            SkewNormalConstraint,
            SkewTConstraint,
            SoftBoundsConstraint,
            SoftLowerBoundConstraint,
        )

        if isinstance(constraint, FixedAgeConstraint):
            return constraint.fixed_age
        elif isinstance(constraint, UniformAgeConstraint):
            return constraint.min_age
        elif isinstance(constraint, SoftLowerBoundConstraint):
            return constraint.min_age
        elif isinstance(constraint, MaximumAgeConstraint):
            return 0.0
        elif isinstance(constraint, SoftBoundsConstraint):
            return constraint.min_age
        elif isinstance(constraint, GammaPriorConstraint):
            try:
                from scipy import stats as sp_stats

                # scipy.stats 的分位数返回 numpy 标量（无标注）；声明的返回类型
                # 是 Optional[float]，这里显式转成 float。
                return float(
                    sp_stats.gamma.ppf(
                        0.025, constraint.alpha, scale=1.0 / constraint.beta
                    )
                    * constraint.scale
                    + constraint.offset
                )
            except Exception as e:
                # Gamma 的精确矩：E = offset + scale·α/β，SD = scale·√α/β
                mean = (
                    constraint.offset
                    + constraint.scale * constraint.alpha / constraint.beta
                )
                sd = constraint.scale * (constraint.alpha**0.5) / constraint.beta
                self.logger.warning(f"scipy.stats.gamma.ppf unavailable ({e})")
                return self._normal_approximation_bound(constraint, mean, sd, "lower")
        elif isinstance(constraint, (SkewNormalConstraint, SkewTConstraint)):
            if isinstance(constraint, SkewNormalConstraint):
                # skewnorm 是 scipy 真实存在的分布，走精确分位数
                try:
                    from scipy import stats as sp_stats

                    return float(
                        sp_stats.skewnorm.ppf(
                            0.025,
                            constraint.shape,
                            loc=constraint.location,
                            scale=constraint.scale,
                        )
                    )
                except Exception as e:
                    self.logger.warning(
                        f"scipy.stats.skewnorm.ppf failed for SkewNormalConstraint ({e})"
                    )
                    delta = constraint.shape / (1 + constraint.shape**2) ** 0.5
                    mean = (
                        constraint.location
                        + constraint.scale * delta * (2 / math.pi) ** 0.5
                    )
                    sd = constraint.scale * (1 - 2 * delta**2 / math.pi) ** 0.5
                    return self._normal_approximation_bound(
                        constraint, mean, sd, "lower"
                    )

            # SkewTConstraint：scipy 从未提供 skew-t 分布（审阅项 C-3：旧代码
            # `sp_stats.skewt.ppf` 恒抛 AttributeError，被 except 吞掉后静默返回
            # location*0.5 这个凭空中点）。与 models/constraints.py 的
            # SkewTConstraint.moment_approximation() 保持一致：Azzalini 精确前两阶
            # 矩 + 显式声明的正态近似，并且必须发 warning。
            moment_approximation = getattr(constraint, "moment_approximation", None)
            if moment_approximation is None:
                self.logger.error(
                    "Ancestry check: SkewTConstraint exposes no "
                    "moment_approximation(), so no age bound can be derived for "
                    "this calibration; the pair is left un-compared."
                )
                return None
            try:
                mean, sd = moment_approximation()
            except Exception as e:
                self.logger.error(
                    f"Ancestry check: skew-t moment approximation failed ({e}); "
                    f"no age bound can be derived for this calibration."
                )
                return None
            return self._normal_approximation_bound(constraint, mean, sd, "lower")

        return None

    def _get_max_age(self, constraint: AgeConstraint) -> Optional[float]:
        """获取约束的最大年龄（97.5% 分位数）"""
        from ..models.constraints import (
            GammaPriorConstraint,
            MaximumAgeConstraint,
            SkewNormalConstraint,
            SkewTConstraint,
            SoftBoundsConstraint,
            SoftLowerBoundConstraint,
        )

        if isinstance(constraint, FixedAgeConstraint):
            return constraint.fixed_age
        elif isinstance(constraint, UniformAgeConstraint):
            return constraint.max_age
        elif isinstance(constraint, SoftLowerBoundConstraint):
            return float("inf")
        elif isinstance(constraint, MaximumAgeConstraint):
            return constraint.max_age
        elif isinstance(constraint, SoftBoundsConstraint):
            return constraint.max_age
        elif isinstance(constraint, GammaPriorConstraint):
            try:
                from scipy import stats as sp_stats

                return float(
                    sp_stats.gamma.ppf(
                        0.975, constraint.alpha, scale=1.0 / constraint.beta
                    )
                    * constraint.scale
                    + constraint.offset
                )
            except Exception as e:
                mean = (
                    constraint.offset
                    + constraint.scale * constraint.alpha / constraint.beta
                )
                sd = constraint.scale * (constraint.alpha**0.5) / constraint.beta
                self.logger.warning(f"scipy.stats.gamma.ppf unavailable ({e})")
                return self._normal_approximation_bound(constraint, mean, sd, "upper")
        elif isinstance(constraint, (SkewNormalConstraint, SkewTConstraint)):
            if isinstance(constraint, SkewNormalConstraint):
                try:
                    from scipy import stats as sp_stats

                    return float(
                        sp_stats.skewnorm.ppf(
                            0.975,
                            constraint.shape,
                            loc=constraint.location,
                            scale=constraint.scale,
                        )
                    )
                except Exception as e:
                    self.logger.warning(
                        f"scipy.stats.skewnorm.ppf failed for SkewNormalConstraint ({e})"
                    )
                    delta = constraint.shape / (1 + constraint.shape**2) ** 0.5
                    mean = (
                        constraint.location
                        + constraint.scale * delta * (2 / math.pi) ** 0.5
                    )
                    sd = constraint.scale * (1 - 2 * delta**2 / math.pi) ** 0.5
                    return self._normal_approximation_bound(
                        constraint, mean, sd, "upper"
                    )

            moment_approximation = getattr(constraint, "moment_approximation", None)
            if moment_approximation is None:
                self.logger.error(
                    "Ancestry check: SkewTConstraint exposes no "
                    "moment_approximation(), so no age bound can be derived for "
                    "this calibration; the pair is left un-compared."
                )
                return None
            try:
                mean, sd = moment_approximation()
            except Exception as e:
                self.logger.error(
                    f"Ancestry check: skew-t moment approximation failed ({e}); "
                    f"no age bound can be derived for this calibration."
                )
                return None
            return self._normal_approximation_bound(constraint, mean, sd, "upper")

        return None


def resolve_special_nodes(
    tree: PhylogeneticTree, resolver: CalibrationResolver
) -> Dict[str, List[str]]:
    """
    解析特殊节点（LUCA、LBCA、LACA）

    Args:
        tree: 系统发育树
        resolver: 解析器

    Returns:
        特殊节点到叶节点列表的映射
    """
    special_nodes = {}

    # LUCA：细菌和古菌的共同祖先
    bacteria = []
    archaea = []

    for name in tree.tip_names:
        name_lower = name.lower()
        if "_d_bacteria" in name_lower:
            bacteria.append(name)
        elif "_d_archaea" in name_lower:
            archaea.append(name)
        else:
            # 尝试从名称推断
            if "bacteria" in name_lower:
                # Methanobacteriota 是细菌门，必须包含在 Bacteria 中
                bacteria.append(name)
            elif "archaea" in name_lower:
                archaea.append(name)

    if bacteria and archaea:
        special_nodes["LUCA"] = list(bacteria) + list(archaea)

    if len(bacteria) >= 2:
        special_nodes["LBCA"] = list(bacteria)

    if len(archaea) >= 2:
        special_nodes["LACA"] = list(archaea)

    return special_nodes
