"""
ComparisonReporter - 结果比较和报告生成

生成比较表格（TXT / TSV）、执行摘要与比较柱状图。

============================  B-15：节点身份必须按拓扑对齐  ============================
各方法 ``DatingResult.node_ages`` 的**键法互不相同**：

  * mcmctree / pathd8 / treepl / r8s / mdcat：已校准节点用**校准点名**
    （``cal.name``），叶节点用短 taxon 名（ACCESSION 模式 ``GB_GCA_..._1``）；
  * r8s / lsd2 / mcmctree 的未命名内部节点：``internal_node_<前序遍历序号>``；
  * 少数回退路径甚至用 ``node_<id(node)>``（Python 对象 id）。

旧实现把**键字符串本身**当成节点身份，于是同一个分支在两个方法里落成两行、
每行只有一列有值：表格显示"无分歧"，真相是"从未比较"（B-15）。反向同样危险：
两个不同分支可能都叫 ``internal_node_1``，于是被并排放到同一行、显示"两方法
一致"（C-41）。

本模块因此按**拓扑**对齐：节点身份 = 该节点的后代叶集合（规范串 + 短哈希）。
每条键的解析顺序是

  1. 该方法自己的 ``dated_tree_newick`` 里同名的节点标签（结构化证据）；
  2. 校准点名 → 由输入树求其 MRCA 的后代叶集合（需调用方传 ``calibrations``）；
  3. 输入树里的节点标签；
  4. 输入树的叶名（含短名/登录号别名）→ 判为叶节点；
  5. 纯位置序号键（``node_37`` / ``internal_node_1`` / ``node_0x7f..``）且无
     任何结构证据 → **unresolved**：单独成行、显式标注"不可比较"，绝不与另一个
     方法的键并排（fail-closed，同时修掉旧的叶节点过滤器对 ``node_<id>`` 失效）；
  6. 其余（校准点名这类"描述性名字"）→ ``name-only`` 对齐：只在字面相同时并排，
     并在报告里标明"未经拓扑核验"。

行的 ``Topology_ID`` / ``Depth`` 两列现在真的来自拓扑（旧实现只对叶节点计算，
而表里按定义已无叶节点，两列恒为"名字截断"和 "N/A"）。

============================  C-16 / C-15 / C-17 / B-26  ============================
* C-16：TSV 每个方法额外写 ``_<method>_ci_type``，把 HPD95 / CI95 / RANGE / NONE /
  NA（无该估计）区分开，读者与下游脚本可以据此判断区间语义。
* C-15：柱状图按区间数据画误差线，缺失的柱子/无区间的方法数量显式披露。
* C-17：可视化失败至少记 warning，成功语列出实际产出的文件。
* B-26：摘要里增加"校准对账（请求 N / 生效 M）"一节，把适配器已经暴露在
  ``result.metadata`` 里的 ``constraint_reconciliation``、
  ``node_age_interval_kinds``、``chain_interval_conflicts``、区间语义分布，以及
  管道侧可选的 ``calibration_reconciliation`` 汇总出来，使"校准静默缩水"可见。
"""

from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

from ..infrastructure import get_logger
from ..infrastructure.safe_io import safe_writer
from ..models import DatingResult, NodeAgeEstimate

# ==================== 节点名归一化辅助 ====================

_NAME_CHARS_RE = re.compile(r"[^A-Za-z0-9_]")
_ACCESSION_PATTERNS = (
    re.compile(r"([A-Z]{2}_[A-Z]{2,3}_\d+\.\d+)"),  # GTDB: GB_GCA_000252485.1
    re.compile(r"(GC[AF]_\d+\.\d+)"),  # GenBank/Assembly: GCA_000252485.1
)
#: PAML 风格的树文件头（"8 1"），不是 Newick 的一部分
_PAML_HEADER_RE = re.compile(r"^\s*\d+\s+\d+\s*$")

#: 只在"位置序号 / 对象 id"意义上出现的键：它们在另一个方法里不代表同一个分支。
_POSITIONAL_KEY_RES = (
    re.compile(r"^(?:internal|unnamed|int|new)[ _\-]?node[ _\-]?\d+$", re.I),
    re.compile(r"^node[ _\-]?\d+$", re.I),
    re.compile(r"^node[ _\-]?0x[0-9a-f]+$", re.I),
    re.compile(r"^(?:nd|iv|an|tc)[ _\-]?\d+$", re.I),
    re.compile(r"^intnode\d+$", re.I),
)


def _sanitize_name(name: str) -> str:
    """与 ``NameMappingManager``/各适配器一致的短名清洗（非字母数字→``_``，截断 30）。"""
    return _NAME_CHARS_RE.sub("_", name)[:30]


def _accession_of(name: str) -> str:
    """按 ``NameMappingManager._extract_accession`` 的同一规则提取登录号短名。"""
    for pattern in _ACCESSION_PATTERNS:
        match = pattern.search(name)
        if match:
            return match.group(1).replace(".", "_")
    return ""


def _name_variants(name: str) -> Tuple[str, ...]:
    """一个名称在跨方法/跨缩写时应尝试的所有写法。"""
    variants = [name, name.lower()]
    sanitized = _sanitize_name(name)
    variants += [sanitized, sanitized.lower()]
    accession = _accession_of(name)
    if accession:
        variants += [accession, accession.lower()]
    # B-21 的短名冲突后缀（GB_GCA_000252485_1_2）：去掉尾部的 _<n> 再试一次
    stripped = re.sub(r"_\d+$", "", name)
    if stripped and stripped != name:
        variants += [stripped, stripped.lower()]
    seen: List[str] = []
    for variant in variants:
        if variant and variant not in seen:
            seen.append(variant)
    return tuple(seen)


def _is_positional_key(name: str) -> bool:
    """键是否只是"位置序号/对象 id"（不具备跨方法身份含义）。"""
    return any(pattern.match(name) for pattern in _POSITIONAL_KEY_RES)


def _clade_token(leaves: Sequence[str]) -> str:
    """后代叶集合的规范串（拓扑身份）。"""
    return "{" + ",".join(sorted(leaves)) + "}"


def _clade_hash(leaves: Sequence[str]) -> str:
    return hashlib.sha1(_clade_token(leaves).encode("utf-8")).hexdigest()[:8]


# ==================== 一个不依赖 ete3 的最小 Newick 读取器 ====================
#
# ete3 在 Python ≥3.13 上根本无法导入（A-5），而 ``PhylogeneticTree`` 的公开 API
# 里没有"枚举所有分支节点及其后代叶集合"这一件（见本文件末尾的跨文件需求说明）。
# 这里只做只读的拓扑遍历（不需要分支长度、不需要外部依赖），失败时上层会退化为
# name-only 对齐并显式披露，不会给出错误的对齐结论。


class _NwNode:
    """Newick 节点（只读拓扑用）。"""

    __slots__ = ("label", "children", "parent", "_leaf_cache")

    def __init__(self, label: str = "") -> None:
        self.label = label
        self.children: List["_NwNode"] = []
        self.parent: Optional["_NwNode"] = None
        self._leaf_cache: Optional[Tuple[str, ...]] = None

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def descendant_leaf_names(self) -> Tuple[str, ...]:
        if self._leaf_cache is None:
            if self.is_leaf:
                self._leaf_cache = (self.label,) if self.label else ()
            else:
                collected: List[str] = []
                for child in self.children:
                    collected.extend(child.descendant_leaf_names())
                self._leaf_cache = tuple(collected)
        return self._leaf_cache


def _strip_annotations(text: str) -> str:
    """去掉 ``[...]`` 注释（NHX ``[&date=...]``、MD-Cat ``[t=...]``、PAML ``[B(...)]``）。

    单引号内的内容（PAML 的校准标注标签）原样保留。
    """
    out: List[str] = []
    depth = 0
    in_single_quote = False
    for ch in text:
        if ch == "'":
            in_single_quote = not in_single_quote
        if not in_single_quote:
            if ch == "[":
                depth += 1
                continue
            if ch == "]":
                if depth:
                    depth -= 1
                continue
            if depth:
                continue
        out.append(ch)
    return "".join(out)


def _parse_newick(text: Optional[str]) -> Optional[_NwNode]:
    """把 Newick 解析成只读拓扑；无法解析时返回 ``None``（上层退化处理）。"""
    if not text:
        return None
    body = _strip_annotations(text)
    body = "".join(
        line for line in body.splitlines() if not _PAML_HEADER_RE.match(line)
    )
    start = body.find("(")
    if start >= 0:
        body = body[start:]
    else:
        body = body.strip().rstrip(";").strip()
    if not body:
        return None

    pos = 0
    length = len(body)

    def parse_node() -> _NwNode:
        nonlocal pos
        node = _NwNode()
        if pos < length and body[pos] == "(":
            pos += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if pos < length and body[pos] == ",":
                    pos += 1
                    continue
                if pos < length and body[pos] == ")":
                    pos += 1
                    break
                break
        label_parts: List[str] = []
        while pos < length and body[pos] not in "(),;:":
            if body[pos] == "'":
                end = body.find("'", pos + 1)
                if end < 0:
                    label_parts.append(body[pos + 1 :])
                    pos = length
                    break
                label_parts.append(body[pos + 1 : end])
                pos = end + 1
                continue
            end = pos
            while end < length and body[end] not in "'(),;:":
                end += 1
            label_parts.append(body[pos:end])
            pos = end
        node.label = "".join(label_parts).strip()
        if pos < length and body[pos] == ":":
            pos += 1
            while pos < length and body[pos] not in "(),;":
                pos += 1
        return node

    try:
        root = parse_node()
    except RecursionError:
        # 极端深（梯状）树：交给上层退化处理，不冒栈溢出的风险
        return None
    if not root.descendant_leaf_names():
        return None
    return root


@dataclass
class _CladeEntry:
    """树里一个节点的结构事实：后代叶集合、深度、是否叶、原始标签。"""

    leaves: Tuple[str, ...]
    depth: int
    is_leaf: bool
    label: str = ""


class _CladeIndex:
    """按拓扑索引一棵树：``标签 → 节点``、``叶集合 → 节点``、``叶名 → 叶节点``。"""

    def __init__(self, newick: Optional[str]) -> None:
        self.ok = False
        self.entries: List[_CladeEntry] = []
        self.by_label: Dict[str, _CladeEntry] = {}
        self.by_leaves: Dict[Tuple[str, ...], _CladeEntry] = {}
        self.tips: List[str] = []
        self.root_entry: Optional[_CladeEntry] = None

        root = _parse_newick(newick)
        if root is None:
            return
        self.ok = True

        pending: List[Tuple[_NwNode, int]] = [(root, 0)]
        while pending:
            node, depth = pending.pop()
            leaves = tuple(sorted(node.descendant_leaf_names()))
            if not leaves:
                continue
            entry = _CladeEntry(
                leaves=leaves, depth=depth, is_leaf=node.is_leaf, label=node.label
            )
            self.entries.append(entry)
            self.by_leaves.setdefault(leaves, entry)
            if depth == 0 and self.root_entry is None:
                self.root_entry = entry
            if node.is_leaf and node.label:
                self.tips.append(node.label)
            if node.label:
                for variant in _name_variants(node.label):
                    existing = self.by_label.get(variant)
                    # 内部节点标签优先于同名叶标签（叶名里出现 "Node1" 的树极少）
                    if existing is None or (existing.is_leaf and not node.is_leaf):
                        self.by_label[variant] = entry
            for child in node.children:
                pending.append((child, depth + 1))

    def lookup(self, name: str) -> Optional[_CladeEntry]:
        if not self.ok or not name:
            return None
        for variant in _name_variants(name):
            entry = self.by_label.get(variant)
            if entry is not None:
                return entry
        return None

    def entry_for_leaves(self, leaves: Tuple[str, ...]) -> Optional[_CladeEntry]:
        return self.by_leaves.get(tuple(sorted(leaves)))

    def mrca_entry(self, leaf_names: Iterable[str]) -> Optional[_CladeEntry]:
        """含给定叶集合的最小节点（MRCA）。任一叶名找不到则返回 None。"""
        if not self.ok:
            return None
        needed = {leaf for leaf in leaf_names if leaf}
        if not needed:
            return None
        best: Optional[_CladeEntry] = None
        for entry in self.entries:
            if len(entry.leaves) >= len(needed) and needed.issubset(entry.leaves):
                if best is None or len(entry.leaves) < len(best.leaves):
                    best = entry
        return best


class _LeafAlias:
    """短名/登录号 ↔ 输入树叶名 的双向别名表（mcmctree 等会改写叶名）。"""

    def __init__(self, tip_names: Sequence[str]) -> None:
        self.map: Dict[str, str] = {}
        for tip in tip_names:
            if not tip:
                continue
            for variant in _name_variants(tip):
                self.map.setdefault(variant, tip)

    def __bool__(self) -> bool:
        return bool(self.map)

    def canonical(self, name: str) -> str:
        """把（可能是短名的）叶名换回输入树里的叶名；认不出则原样返回。"""
        if not self.map:
            return name
        for variant in _name_variants(name):
            tip = self.map.get(variant)
            if tip is not None:
                return tip
        return name

    def as_tip(self, name: str) -> Optional[str]:
        """名称确实指某个叶节点（含短名/登录号别名）时返回其规范叶名。"""
        if not self.map:
            return None
        for variant in _name_variants(name):
            tip = self.map.get(variant)
            if tip is not None:
                return tip
        return None


# ==================== 对齐结果的数据结构 ====================

ALIGN_TOPOLOGY = "topology"
ALIGN_NAME_ONLY = "name-only"
ALIGN_UNRESOLVED = "unresolved"
ALIGN_LEAF = "leaf"


@dataclass
class _Row:
    """比较表里的一行 = 一个（尽可能）拓扑对齐后的节点。"""

    node_id: str
    display: str
    alignment: str
    leaves: Tuple[str, ...] = ()
    depth: Optional[int] = None
    depth_source: str = "none"
    comparable: bool = True
    values: Dict[str, NodeAgeEstimate] = field(default_factory=dict)
    source_keys: Dict[str, str] = field(default_factory=dict)

    @property
    def n_leaves(self) -> int:
        return len(self.leaves)

    @property
    def topology_id(self) -> str:
        if self.alignment != ALIGN_TOPOLOGY or not self.leaves:
            return "unresolved"
        return f"{len(self.leaves)}L:{_clade_hash(self.leaves)}"

    def methods_present(self) -> List[str]:
        return sorted(self.values)


@dataclass
class _Comparison:
    """一次对齐的完整产物，供 TXT / TSV / 图 / 摘要共用。"""

    methods: List[str]
    rows: List[_Row] = field(default_factory=list)
    leaf_entries: List[Tuple[str, str]] = field(default_factory=list)
    unresolved_entries: List[Tuple[str, str]] = field(default_factory=list)
    alignment_counts: Dict[str, int] = field(default_factory=dict)
    unmapped_leaf_names: int = 0
    tree_available: bool = False
    newick_indexes: int = 0
    ci_type_counts: Dict[str, Dict[str, int]] = field(default_factory=dict)

    @property
    def comparable_rows(self) -> List[_Row]:
        return [r for r in self.rows if r.comparable]

    @property
    def multi_method_rows(self) -> List[_Row]:
        return [r for r in self.comparable_rows if len(r.values) >= 2]

    def rows_missing_by_method(self) -> Dict[str, int]:
        return {
            method: sum(1 for row in self.comparable_rows if method not in row.values)
            for method in self.methods
        }


def _tree_newick(tree: Any) -> Optional[str]:
    """从 ``PhylogeneticTree``（或直接给 Newick 字符串）取 Newick 文本。"""
    if tree is None:
        return None
    if isinstance(tree, str):
        return tree
    newick = getattr(tree, "newick", None)
    if isinstance(newick, str) and newick.strip():
        return newick
    to_string = getattr(tree, "to_nhx_string", None)
    if callable(to_string):
        try:
            value = to_string()
            if isinstance(value, str):
                return value
        except Exception as exc:  # pragma: no cover - 防御
            get_logger().debug(f"tree.to_nhx_string() 失败，改用其它途径取拓扑：{exc}")
    return None


def _tree_tips(tree: Any, index: Optional[_CladeIndex]) -> List[str]:
    tips: List[str] = []
    seen = set()
    sources: List[str] = []
    if index is not None:
        sources.extend(index.tips)
    if tree is not None and not isinstance(tree, str):
        try:
            sources.extend(list(tree.tip_names))  # 公开 API
        except Exception as exc:  # pragma: no cover - 防御
            get_logger().debug(f"tree.tip_names 读取失败：{exc}")
    for name in sources:
        if name and name not in seen:
            seen.add(name)
            tips.append(name)
    return tips


def _ci_type_label(estimate: Optional[NodeAgeEstimate]) -> str:
    """区间语义标签：HPD95 / CI95 / RANGE / NONE / UNSPECIFIED / NA(无该估计)。

    ``UNSPECIFIED`` 专门留给"给了上下界却没声明区间类型"的估计——它既不是"没有
    区间"（NONE），也不是"缺数据"（NA），把它并入任何一类都会骗到下游读者（C-16）。
    """
    if estimate is None:
        return "NA"
    ci_type = estimate.ci_type
    value = getattr(ci_type, "value", ci_type)
    label = "NONE" if value in (None, "") else str(value)
    if label == "NONE" and (
        estimate.ci_lower is not None or estimate.ci_upper is not None
    ):
        return "UNSPECIFIED"
    return label


def _has_interval(estimate: NodeAgeEstimate) -> bool:
    return (
        estimate.ci_lower is not None
        and estimate.ci_upper is not None
        and _ci_type_label(estimate) not in ("NONE", "NA")
    )


def _require_bounds(estimate: NodeAgeEstimate) -> Tuple[float, float]:
    """取成对的置信区间上下界。

    :func:`_has_interval` 已经把"任一界为 None"排除在外（它额外还校了
    ci_type 标签），本函数只负责把那个前提翻译成本地可用的 ``(float, float)``：
    mypy 不会跟进函数调用去推断属性非空。前提被破坏时招错，而不是
    默默把缺失的界当 0 参与比较。
    """
    if estimate.ci_lower is None or estimate.ci_upper is None:
        raise ValueError(
            "NodeAgeEstimate has no paired confidence bounds; guard the call with "
            "_has_interval() first"
        )
    return estimate.ci_lower, estimate.ci_upper


def _format_estimate(estimate: Optional[NodeAgeEstimate]) -> str:
    """与 ``DatingResult.format_age_string`` 同口径：``48.3 [41.2,55.6]_HPD95``。"""
    if estimate is None:
        return "NA"
    if _has_interval(estimate):
        return (
            f"{estimate.mean_age:.1f} "
            f"[{estimate.ci_lower:.1f},{estimate.ci_upper:.1f}]_"
            f"{_ci_type_label(estimate)}"
        )
    return f"{estimate.mean_age:.1f}"


class ComparisonReporter:
    """
    结果比较报告生成器

    生成：
    - comparison_table.txt: 文本比较表格（按拓扑对齐，含区间语义）
    - comparison_table.tsv: 机器可读表格（每方法带 _ci_type / _source_key）
    - summary.txt: 执行摘要（对齐口径 + 校准对账）
    - comparison_plot.png: 比较柱状图（带误差线，来自真实区间）
    """

    #: 摘要里做校准对账时读取的 metadata 键
    RECONCILIATION_KEYS = (
        "calibration_reconciliation",
        "constraint_reconciliation",
        "node_age_interval_kinds",
        "chain_interval_conflicts",
        "clock_tests",
    )

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir
        self.logger = get_logger()
        self._tree = None
        self._last_comparison: Optional[_Comparison] = None
        self._plot_stats: Dict[str, Any] = {}
        self._failures: List[str] = []

    # ------------------------------------------------------------------ 入口

    def generate(
        self,
        results: Dict[str, DatingResult],
        tree: Any = None,
        calibrations: Optional[Sequence[Any]] = None,
    ) -> Dict[str, Any]:
        """生成所有报告。

        Args:
            results: 定年结果字典（``{method_name: DatingResult}``）
            tree: 输入树（``PhylogeneticTree`` 或 Newick 字符串），用于按拓扑对齐
            calibrations: 可选的校准点列表；提供后 ``cal.name`` 这类键也能被解析成
                真实的 clade（而不是只靠字符串相等），从而进一步收紧 B-15 的对齐。

        Returns:
            本次生成的文件清单与告警计数。
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._tree = tree
        self._failures = []

        comparison = self._align(results, tree, calibrations)
        self._last_comparison = comparison

        generated: List[Path] = []
        generated.append(
            self._generate_text_table(results, tree=tree, comparison=comparison)
        )
        generated.append(
            self._generate_tsv_table(results, tree=tree, comparison=comparison)
        )
        generated.append(
            self._generate_summary(results, tree=tree, comparison=comparison)
        )

        # 可视化（matplotlib 可选；C-17：失败必须可见）
        try:
            plot = self.generate_visualization(
                results, tree=tree, comparison=comparison
            )
        except Exception as exc:  # noqa: BLE001 - 图形失败不应吃掉整份报告
            plot = None
            self._failures.append(f"visualization: {exc}")
            self.logger.warning(
                f"Comparison bar chart was NOT generated ({exc}). "
                "The TSV/TXT tables are unaffected; re-run with matplotlib available "
                "if the figure is needed."
            )
        if plot is not None:
            generated.append(plot)

        files = ", ".join(path.name for path in generated) or "(none)"
        self.logger.success(
            f"Generated comparison reports in {self.output_dir}: {files}"
        )
        if self._failures:
            self.logger.warning(
                f"Comparison report completed with {len(self._failures)} degradation(s): "
                + "; ".join(self._failures)
            )
        return {
            "files": [str(path) for path in generated],
            "alignments": dict(comparison.alignment_counts),
            "unresolved_keys": len(comparison.unresolved_entries),
            "plot": self._plot_stats,
            "failures": list(self._failures),
        }

    # ------------------------------------------------ 向后兼容的私有小工具

    def _is_leaf_node(self, node_name: str, tree: Any) -> bool:
        """判断节点键是否对应叶节点（末端类群）——现在对 ``node_<id>`` 同样有效。

        叶节点行的"年龄"其实是根到叶的距离，不能与节点年龄混在一张表里。旧实现
        对 ``node_<id>`` 一律返回 ``False``（B-15(2)），恰好让最需要过滤的那批方法
        漏了进来。现改为走同一套拓扑解析：能确定是叶 → True；是纯位置键且解析不
        出来 → 也按"不可比较"排除（fail-closed），并在对齐结果里留痕。
        """
        if tree is None and not self._last_comparison:
            return False
        comparison = self._last_comparison or self._align({}, tree, None)
        for method, key in comparison.leaf_entries:
            if key == node_name:
                return True
        # 名称层面直接就是输入树的叶名（含短名/登录号别名）
        input_index = _CladeIndex(_tree_newick(tree))
        entry = input_index.lookup(node_name) if input_index.ok else None
        return bool(entry and entry.is_leaf)

    def _get_topology_id_for_node(self, node_name: str, tree: Any) -> str:
        """节点的拓扑标识（``<叶数>L:<叶集合哈希>``），解析不出时显式 "unresolved"。"""
        row = self._find_row(node_name)
        if row is not None:
            return row.topology_id
        input_index = _CladeIndex(_tree_newick(tree))
        if input_index.ok:
            entry = input_index.lookup(node_name)
            if entry is not None and entry.leaves:
                return f"{len(entry.leaves)}L:{_clade_hash(entry.leaves)}"
        return "unresolved"

    def _get_depth_for_node(self, node_name: str, tree: Any) -> str:
        """节点深度（输入树里到根的边数），解析不出时显式 "N/A" 并留痕。"""
        row = self._find_row(node_name)
        if row is not None and row.depth is not None:
            return str(row.depth)
        input_index = _CladeIndex(_tree_newick(tree))
        if input_index.ok:
            entry = input_index.lookup(node_name)
            if entry is not None:
                return str(entry.depth)
        return "N/A"

    def _find_row(self, node_name: str) -> Optional[_Row]:
        if not self._last_comparison:
            return None
        for row in self._last_comparison.rows:
            if row.display == node_name or node_name in row.source_keys.values():
                return row
        return None

    # ------------------------------------------------------------ 拓扑对齐

    def _align(
        self,
        results: Dict[str, DatingResult],
        tree: Any = None,
        calibrations: Optional[Sequence[Any]] = None,
    ) -> _Comparison:
        methods = list(results.keys())
        # “不 ok 就当没有”是这个模块的通用表达，所以这里必须是 Optional。
        # 声明与赋值分开写：``x: Optional[T] = T()`` 不会被收窄。
        input_index: Optional[_CladeIndex]
        input_index = _CladeIndex(_tree_newick(tree))
        if not input_index.ok:
            input_index = None
        alias = _LeafAlias(_tree_tips(tree, input_index))
        cal_clades = self._build_calibration_clades(calibrations, input_index, alias)

        comparison = _Comparison(
            methods=methods,
            tree_available=input_index is not None,
            newick_indexes=0,
        )
        rows: Dict[str, _Row] = {}

        for method, result in results.items():
            dated_index: Optional[_CladeIndex]
            dated_index = _CladeIndex(getattr(result, "dated_tree_newick", None))
            if dated_index.ok:
                comparison.newick_indexes += 1
            else:
                dated_index = None

            for key, estimate in (getattr(result, "node_ages", {}) or {}).items():
                row = self._resolve_and_place(
                    method=method,
                    key=str(key),
                    estimate=estimate,
                    dated_index=dated_index,
                    input_index=input_index,
                    alias=alias,
                    cal_clades=cal_clades,
                    rows=rows,
                    comparison=comparison,
                )
                if row is None:
                    continue

        comparison.rows = self._order_rows(rows.values())
        comparison.alignment_counts = {}
        for row in comparison.rows:
            comparison.alignment_counts[row.alignment] = (
                comparison.alignment_counts.get(row.alignment, 0) + 1
            )
        comparison.ci_type_counts = {
            method: self._ci_type_distribution(
                (getattr(results[method], "node_ages", {}) or {}).values()
            )
            for method in methods
        }

        if comparison.unresolved_entries:
            self.logger.warning(
                f"{len(comparison.unresolved_entries)} node key(s) could not be resolved "
                "to a clade from the tree and are marked uncomparable "
                "(they are NOT aligned with other methods' nodes): "
                + "; ".join(f"{m}:{k}" for m, k in comparison.unresolved_entries[:8])
                + (" …" if len(comparison.unresolved_entries) > 8 else "")
                + ". To close the gap: pass the input tree AND the resolved calibration "
                "list into ComparisonReporter.generate(..., tree=..., calibrations=...), "
                "and/or have adapters keep node labels on the dated tree so positional "
                "keys (node_<n> / internal_node_<n>) become resolvable."
            )
        if comparison.leaf_entries:
            self.logger.info(
                f"Dropped {len(comparison.leaf_entries)} leaf-node entries "
                "(root-to-tip distances, not node ages) from the comparison table."
            )
        if comparison.unmapped_leaf_names:
            self.logger.warning(
                f"{comparison.unmapped_leaf_names} leaf name(s) in the method output "
                "trees could not be mapped back to the input tree; clade identities "
                "built from them may not match other methods."
            )
        return comparison

    def _resolve_and_place(
        self,
        method: str,
        key: str,
        estimate: NodeAgeEstimate,
        dated_index: Optional[_CladeIndex],
        input_index: Optional[_CladeIndex],
        alias: _LeafAlias,
        cal_clades: Dict[str, Tuple[str, ...]],
        rows: Dict[str, _Row],
        comparison: _Comparison,
    ) -> Optional[_Row]:
        leaves: Optional[Tuple[str, ...]] = None
        depth: Optional[int] = None
        in_method_tree = False

        entry = dated_index.lookup(key) if dated_index is not None else None
        if entry is not None:
            leaves = entry.leaves
            depth = entry.depth
            in_method_tree = True
        elif key in cal_clades:
            leaves = cal_clades[key]
        elif input_index is not None:
            input_entry = input_index.lookup(key)
            if input_entry is not None:
                leaves = input_entry.leaves
                depth = input_entry.depth

        # 单叶 = 叶节点：它的"年龄"是根到叶距离，不进比较表（B-15(2) 的过滤器）
        if leaves is not None and len(leaves) == 1:
            comparison.leaf_entries.append((method, key))
            return None
        if leaves is None and alias.as_tip(key) is not None:
            # 键本身就是（可能被缩写的）叶名 —— 同样按叶节点排除，fail-closed
            comparison.leaf_entries.append((method, key))
            return None

        if leaves is None:
            if _is_positional_key(key):
                # 位置序号冒充节点身份：绝不与其它方法的键并排（B-15/C-41）
                comparison.unresolved_entries.append((method, key))
                node_id = f"unresolved:{method}:{key}"
                row = rows.get(node_id)
                if row is None:
                    row = _Row(
                        node_id=node_id,
                        display=key,
                        alignment=ALIGN_UNRESOLVED,
                        comparable=False,
                    )
                    rows[node_id] = row
                self._attach(row, method, key, estimate)
                return row
            # 描述性名字（多数情况下就是校准点名）：只能按名字对齐，未核验拓扑
            node_id = f"name:{key}"
            row = rows.get(node_id)
            if row is None:
                row = _Row(
                    node_id=node_id,
                    display=key,
                    alignment=ALIGN_NAME_ONLY,
                )
                rows[node_id] = row
            self._attach(row, method, key, estimate)
            return row

        # 把方法树里的叶名换成输入树的叶名，使短名方法与全名方法可比
        mapped: List[str] = []
        for leaf in leaves:
            canonical = alias.canonical(leaf)
            if alias.map and leaf not in alias.map:
                comparison.unmapped_leaf_names += 1
            mapped.append(canonical)
        leaves = tuple(sorted(set(mapped)))

        if len(leaves) <= 1:
            comparison.leaf_entries.append((method, key))
            return None

        # 输入树里存在同一 clade → 深度按输入树（规范口径）计
        depth_source = "method-tree" if in_method_tree else "inferred"
        if input_index is not None:
            input_entry = input_index.entry_for_leaves(leaves)
            if input_entry is not None:
                depth, depth_source = input_entry.depth, "input-tree"

        node_id = "clade:" + _clade_hash(leaves)
        row = rows.get(node_id)
        if row is None:
            row = _Row(
                node_id=node_id,
                display=key,
                alignment=ALIGN_TOPOLOGY,
                leaves=leaves,
                depth=depth,
                depth_source=depth_source,
            )
            rows[node_id] = row
        else:
            # 后到的、来自输入树的深度是更权威的口径，覆盖方法树推断值
            if depth_source == "input-tree" and row.depth_source != "input-tree":
                row.depth, row.depth_source = depth, "input-tree"
            elif row.depth is None and depth is not None:
                row.depth, row.depth_source = depth, depth_source
        self._attach(row, method, key, estimate)
        return row

    @staticmethod
    def _attach(row: _Row, method: str, key: str, estimate: NodeAgeEstimate) -> None:
        previous_key = row.source_keys.get(method)
        if previous_key is not None and previous_key != key:
            # 同一方法有两个键落进同一个 clade：要么树里有两个同名/等价节点，
            # 要么该方法自己的键法内部就不自洽。静默取舍会伪装成"只有一个估计"。
            get_logger().warning(
                f"{method} contributed two node keys to the same clade "
                f"'{row.topology_id}': '{previous_key}' (kept) and '{key}' "
                f"(values {row.values[method].mean_age:.3f} vs "
                f"{estimate.mean_age:.3f} Ma). The later one is kept; check this "
                "method's node-key naming."
            )
        row.values[method] = estimate
        row.source_keys[method] = key

    @staticmethod
    def _order_rows(rows: Iterable[_Row]) -> List[_Row]:
        rank = {ALIGN_TOPOLOGY: 0, ALIGN_NAME_ONLY: 1, ALIGN_UNRESOLVED: 2}
        return sorted(
            rows,
            key=lambda r: (
                rank.get(r.alignment, 3),
                r.depth if r.depth is not None else 9999,
                -r.n_leaves,
                r.display.lower(),
            ),
        )

    @staticmethod
    def _build_calibration_clades(
        calibrations: Optional[Sequence[Any]],
        input_index: Optional[_CladeIndex],
        alias: _LeafAlias,
    ) -> Dict[str, Tuple[str, ...]]:
        """``cal.name → 该校准节点的后代叶集合``（用输入树的真实拓扑求 MRCA）。"""
        clades: Dict[str, Tuple[str, ...]] = {}
        if not calibrations or input_index is None:
            return clades
        for cal in calibrations:
            name = getattr(cal, "name", None)
            if not name:
                continue
            if getattr(cal, "is_root_node", False):
                if input_index.root_entry is not None:
                    clades[name] = input_index.root_entry.leaves
                continue
            pair = getattr(cal, "mrca_leaf_pair", None)
            taxa = getattr(cal, "resolved_taxa", None)
            for candidates in (pair, taxa):
                if not candidates:
                    continue
                leaves = [alias.canonical(str(t)) for t in candidates if t]
                entry = input_index.mrca_entry(leaves)
                if entry is not None and not entry.is_leaf:
                    clades[name] = entry.leaves
                    break
        return clades

    @staticmethod
    def _ci_type_distribution(estimates: Iterable[NodeAgeEstimate]) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for estimate in estimates:
            label = _ci_type_label(estimate)
            counts[label] = counts.get(label, 0) + 1
        return counts

    # ------------------------------------------------------------ TXT 表格

    def _generate_text_table(
        self,
        results: Dict[str, DatingResult],
        tree: Any = None,
        comparison: Optional[_Comparison] = None,
    ) -> Path:
        """生成文本比较表格（行按拓扑对齐，两列拓扑信息是真的）。"""
        output_file = self.output_dir / "comparison_table.txt"
        comparison = comparison or self._align(results, tree)
        rows = comparison.comparable_rows
        methods = comparison.methods

        lines: List[str] = []
        lines.append("=" * 150)
        lines.append("PhyloDater Comparison Table")
        lines.append("=" * 150)
        lines.append("")
        lines.append(
            "Rows are aligned across methods by NODE TOPOLOGY (descendant leaf set), "
            "not by the raw node-key string each method emits (review B-15)."
        )
        if not comparison.tree_available:
            lines.append(
                "WARNING: no input tree was supplied — clade identities were derived from "
                "each method's own dated tree where possible; other rows are name-only."
            )
        lines.append("")

        header = f"{'Node':<26} | {'Topology_ID':<18} | {'Depth':<6} | {'Algn':<10}"
        for method in methods:
            header += f" | {method:<26}"
        lines.append(header)
        lines.append("-" * len(header))

        for row in rows:
            depth = str(row.depth) if row.depth is not None else "N/A"
            line = (
                f"{row.display:<26} | {row.topology_id:<18} | {depth:<6} | "
                f"{row.alignment:<10}"
            )
            for method in methods:
                line += f" | {_format_estimate(row.values.get(method)):<26}"
            lines.append(line)

        if not rows:
            lines.append(
                "(No internal/calibration nodes to display after filtering leaf nodes)"
            )

        lines.append("")
        lines.append("=" * 150)
        lines.append("Node identity / alignment legend")
        lines.append("=" * 150)
        lines.append(
            "  Topology_ID = <#descendant leaves>L:<sha1-8 of sorted leaf set>; "
            "Depth = edges from the root; Algn = how the row was aligned."
        )
        lines.append(
            "  topology   : aligned by descendant leaf set (safe to compare side by side)."
        )
        lines.append(
            "  name-only  : keys carry no structural evidence; aligned by identical "
            "name only — NOT topology-verified."
        )
        lines.append(
            "  unresolved : positional key (node_37 / internal_node_1 / node_<id>) that "
            "could not be resolved to a clade — UNCOMPARABLE, never merged."
        )
        lines.append("")
        lines.append("Per-row source keys (which raw key each method contributed):")
        for row in rows:
            keys = " | ".join(
                f"{method}='{row.source_keys.get(method, '—')}'" for method in methods
            )
            clade = (
                ",".join(row.leaves[:6]) + ("…" if len(row.leaves) > 6 else "")
                if row.leaves
                else "(unresolved)"
            )
            lines.append(f"  {row.display}: leaves={clade}; {keys}")
        lines.append("")
        lines.extend(self._alignment_summary_lines(comparison))
        lines.append("")

        uncomparable = [row for row in comparison.rows if not row.comparable]
        if uncomparable:
            lines.append("=" * 150)
            lines.append(
                "Uncomparable entries (kept out of the table above; shown so no data is "
                "silently lost). These node keys carry only positional meaning, so they "
                "can NEVER be placed next to another method's node."
            )
            lines.append("-" * 150)
            for row in uncomparable:
                for method in row.methods_present():
                    lines.append(
                        f"  {method}:'{row.source_keys[method]}' = "
                        f"{_format_estimate(row.values[method])}"
                    )
            lines.append("")

        lines.append("=" * 150)

        with safe_writer(output_file, encoding="utf-8") as handle:
            handle.write("\n".join(lines))
        self.logger.info(f"Generated text table: {output_file}")
        return output_file

    def _alignment_summary_lines(self, comparison: _Comparison) -> List[str]:
        rows = comparison.comparable_rows
        multi = comparison.multi_method_rows
        counts = comparison.alignment_counts
        lines = ["Alignment summary (B-15):"]
        lines.append(
            f"  aligned rows shown: {len(rows)} "
            f"(topology={counts.get(ALIGN_TOPOLOGY, 0)}, "
            f"name-only={counts.get(ALIGN_NAME_ONLY, 0)}); "
            f"uncomparable rows: {counts.get(ALIGN_UNRESOLVED, 0)}"
        )
        lines.append(
            f"  rows carrying estimates from >=2 methods (genuinely compared): "
            f"{len(multi)}; rows with a single method (no comparison possible): "
            f"{len(rows) - len(multi)}"
        )
        lines.append(
            f"  leaf-node entries dropped as root-to-tip distances: "
            f"{len(comparison.leaf_entries)}"
        )
        if comparison.unresolved_entries:
            lines.append(
                "  unresolved keys (listed per method, never aligned): "
                + "; ".join(f"{m}:{k}" for m, k in comparison.unresolved_entries[:20])
            )
        for method, missing in sorted(comparison.rows_missing_by_method().items()):
            if missing:
                lines.append(
                    f"  {method}: no estimate for {missing}/{len(rows)} aligned row(s) "
                    "(blank cells below, not 'agreement')"
                )
        disagreement = self._disagreement_lines(comparison)
        if disagreement:
            lines.append("")
            lines.extend(disagreement)
        return lines

    def _disagreement_lines(self, comparison: _Comparison) -> List[str]:
        rows = comparison.multi_method_rows
        if not rows:
            return []
        worst: Optional[Tuple[float, str, str, float, float]] = None
        non_overlapping = 0
        mixed_semantics = 0
        for row in rows:
            present = sorted(row.values)
            for i, method_a in enumerate(present):
                for method_b in present[i + 1 :]:
                    est_a = row.values[method_a]
                    est_b = row.values[method_b]
                    delta = abs(est_a.mean_age - est_b.mean_age)
                    if worst is None or delta > worst[0]:
                        worst = (
                            delta,
                            row.display,
                            f"{method_a} vs {method_b}",
                            est_a.mean_age,
                            est_b.mean_age,
                        )
                    if _has_interval(est_a) and _has_interval(est_b):
                        lower_a, upper_a = _require_bounds(est_a)
                        lower_b, upper_b = _require_bounds(est_b)
                        overlap = not (upper_a < lower_b or upper_b < lower_a)
                        if not overlap:
                            non_overlapping += 1
                        if _ci_type_label(est_a) != _ci_type_label(est_b):
                            mixed_semantics += 1
        lines = [
            f"Cross-method disagreement over {len(rows)} compared row(s) "
            "(topology-aligned only):"
        ]
        if worst is not None:
            lines.append(
                f"  largest mean difference: {worst[0]:.1f} Ma at '{worst[1]}' "
                f"({worst[2]}: {worst[3]:.1f} vs {worst[4]:.1f} Ma)"
            )
        lines.append(f"  pairs with non-overlapping intervals: {non_overlapping}")
        if mixed_semantics:
            lines.append(
                f"  NOTE: {mixed_semantics} of those comparisons mix interval kinds "
                "(HPD95 vs CI95 vs RANGE); they are not statistically equivalent (C-16)."
            )
        return lines

    # ------------------------------------------------------------ TSV 表格

    def _generate_tsv_table(
        self,
        results: Dict[str, DatingResult],
        tree: Any = None,
        comparison: Optional[_Comparison] = None,
    ) -> Path:
        """生成 TSV：区间语义逐列保留（C-16），对齐口径逐行保留（B-15）。"""
        output_file = self.output_dir / "comparison_table.tsv"
        comparison = comparison or self._align(results, tree)
        rows = comparison.comparable_rows + [
            row for row in comparison.rows if not row.comparable
        ]
        methods = comparison.methods

        with safe_writer(output_file, newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            header = [
                "Node",
                "Node_ID",
                "Topology_ID",
                "Alignment",
                "Comparable",
                "Depth",
                "Depth_Source",
                "N_Leaves",
                "Clade_Leaves",
                "Source_Keys",
            ]
            for method in methods:
                # C-16: mean / 区间端点之外，必须带上区间类型与该方法的原始键
                header.extend(
                    [
                        f"{method}_mean",
                        f"{method}_ci_lower",
                        f"{method}_ci_upper",
                        f"{method}_ci_type",
                        f"{method}_source_key",
                    ]
                )
            writer.writerow(header)

            for row in rows:
                source_keys = ";".join(
                    f"{method}={row.source_keys[method]}"
                    for method in methods
                    if method in row.source_keys
                )
                line: List[Union[int, str, float]] = [
                    row.display,
                    row.node_id,
                    row.topology_id,
                    row.alignment,
                    "TRUE" if row.comparable else "FALSE",
                    row.depth if row.depth is not None else "NA",
                    row.depth_source,
                    row.n_leaves if row.leaves else "NA",
                    "|".join(row.leaves) if row.leaves else "NA",
                    source_keys or "NA",
                ]
                for method in methods:
                    estimate = row.values.get(method)
                    if estimate is None:
                        line.extend(["NA", "NA", "NA", "NA", "NA"])
                    else:
                        line.extend(
                            [
                                estimate.mean_age,
                                (
                                    estimate.ci_lower
                                    if estimate.ci_lower is not None
                                    else "NA"
                                ),
                                (
                                    estimate.ci_upper
                                    if estimate.ci_upper is not None
                                    else "NA"
                                ),
                                _ci_type_label(estimate),
                                row.source_keys.get(method, "NA"),
                            ]
                        )
                writer.writerow(line)

        self.logger.info(f"Generated TSV table: {output_file}")
        return output_file

    # ------------------------------------------------------------ 摘要

    def _generate_summary(
        self,
        results: Dict[str, DatingResult],
        tree: Any = None,
        comparison: Optional[_Comparison] = None,
    ) -> Path:
        output_file = self.output_dir / "summary.txt"
        comparison = comparison or self._align(results, tree)

        unique_clades = {row.node_id for row in comparison.rows if row.leaves}
        lines: List[str] = []
        lines.append("=" * 80)
        lines.append("PhyloDater Execution Summary")
        lines.append("=" * 80)
        lines.append("")

        lines.append(f"Total methods executed: {len(results)}")
        # 去重后的节点数：旧实现按"方法 × 键"计数，多方法会把同一支系重复计入
        lines.append(
            f"Distinct nodes dated (topology-deduplicated): {len(unique_clades)}"
        )
        lines.append(
            f"Distinct nodes comparable across methods: "
            f"{len(comparison.comparable_rows)}"
        )
        lines.append(
            f"Node keys recorded by all methods (raw, may double-count one branch): "
            f"{sum(len(r.node_ages) for r in results.values())}"
        )
        if comparison.unresolved_entries:
            lines.append(
                f"Uncomparable node keys (not resolvable to a clade): "
                f"{len(comparison.unresolved_entries)}"
            )
        lines.append("")

        lines.extend(self._alignment_summary_lines(comparison))
        lines.append("")

        lines.extend(self._reconciliation_lines(results, comparison))
        lines.append("")

        lines.append("Method Statistics:")
        lines.append("-" * 80)
        for method_name, result in results.items():
            raw_keys = list((getattr(result, "node_ages", {}) or {}).keys())
            attached = [
                row for row in comparison.rows if method_name in row.source_keys
            ]
            internal = [row for row in attached if row.comparable]
            dropped_leaf = sum(
                1 for m, _ in comparison.leaf_entries if m == method_name
            )
            unresolved = sum(
                1 for m, _ in comparison.unresolved_entries if m == method_name
            )
            lines.append(f"\n{method_name}:")
            lines.append(f"  Execution time: {result.execution_seconds:.1f}s")
            lines.append(f"  Node keys emitted: {len(raw_keys)}")
            lines.append(f"  Nodes shown in comparison table: {len(internal)}")
            lines.append(
                f"  Entries dropped as leaf (root-to-tip) nodes: {dropped_leaf}"
            )
            lines.append(f"  Entries unresolvable to a clade: {unresolved}")
            counts = comparison.ci_type_counts.get(method_name, {})
            if counts:
                lines.append(
                    "  Interval semantics of emitted nodes: "
                    + ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
                )

            if result.is_converged is not None:
                status = "Converged" if result.is_converged else "Not converged"
                lines.append(f"  Convergence: {status}")

            if result.warnings:
                lines.append(f"  Warnings: {len(result.warnings)}")
                for warning in result.warnings:
                    lines.append(f"    - {warning}")

            other_metadata = {
                key: value
                for key, value in (getattr(result, "metadata", None) or {}).items()
                if key not in self.RECONCILIATION_KEYS
            }
            if other_metadata:
                lines.append("  Metadata:")
                for key, value in other_metadata.items():
                    lines.append(f"    {key}: {value}")

        lines.append("")
        lines.append("=" * 80)

        with safe_writer(output_file, encoding="utf-8") as handle:
            handle.write("\n".join(lines))
        self.logger.info(f"Generated summary: {output_file}")
        return output_file

    def _reconciliation_lines(
        self,
        results: Dict[str, DatingResult],
        comparison: _Comparison,
    ) -> List[str]:
        """B-26：把"请求 N / 生效 M"的校准对账摊到报告里。"""
        lines = [
            "Calibration reconciliation (requested N / effective M, B-26):",
            "  Adapters that lose or soften constraints during input writing now report",
            "  it in result.metadata; without this section a shrunk calibration set is invisible.",
        ]
        anything = False
        for method, result in results.items():
            metadata = getattr(result, "metadata", None) or {}
            buckets: List[str] = []

            recon = metadata.get("calibration_reconciliation") or metadata.get(
                "constraint_reconciliation"
            )
            if isinstance(recon, dict):
                requested = recon.get("requested", recon.get("calibrations_supplied"))
                effective = recon.get(
                    "effective", recon.get("constrained_nodes_written")
                )
                if requested is not None or effective is not None:
                    buckets.append(
                        f"calibrations requested={requested} / "
                        f"effective={effective}"
                    )
                    anything = True
                not_written = recon.get("calibrations_not_written")
                if not_written:
                    buckets.append(f"not written to input file: {list(not_written)}")
                written = recon.get("directives_written")
                if written:
                    buckets.append(f"directives written: {dict(written)}")

                # 请求数与真实节点数不一致也要能看见
                if isinstance(requested, int) and requested:
                    emitted = len(getattr(result, "node_ages", {}) or {})
                    if emitted < requested:
                        buckets.append(
                            f"only {emitted} node age(s) parsed for {requested} "
                            "requested calibration(s)"
                        )
                        anything = True

            interval_kinds = metadata.get("node_age_interval_kinds")
            if isinstance(interval_kinds, dict) and interval_kinds:
                grouped: Dict[str, int] = {}
                for description in interval_kinds.values():
                    grouped[str(description)] = grouped.get(str(description), 0) + 1
                buckets.append(
                    "node interval kinds: "
                    + "; ".join(f"{k} ×{v}" for k, v in sorted(grouped.items()))
                )
                anything = True

            conflicts = metadata.get("chain_interval_conflicts")
            if isinstance(conflicts, dict) and conflicts:
                buckets.append(
                    f"chain interval conflicts on {len(conflicts)} node(s): "
                    + ", ".join(sorted(str(k) for k in conflicts)[:10])
                )
                anything = True

            if isinstance(metadata.get("clock_tests"), dict):
                clock = metadata["clock_tests"]
                buckets.append(
                    f"clock tests accepted={clock.get('accepted')}/{clock.get('total')}"
                )
                anything = True

            counts = comparison.ci_type_counts.get(method, {})
            mixed = {k for k in counts if k not in ("NONE", "NA")}
            if len(mixed) > 1:
                buckets.append(
                    "MIXED interval semantics across nodes: "
                    + ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
                    + " — do not aggregate these intervals arithmetically (C-16)."
                )
                anything = True
            if "RANGE" in counts:
                buckets.append(
                    f"{counts['RANGE']} node(s) carry RANGE (mean±SD / envelope) "
                    "intervals labelled as if they were credible intervals — see "
                    "metadata['node_age_interval_kinds']."
                )
                anything = True

            if buckets:
                anything = True
                lines.append(f"  {method}:")
                for bucket in buckets:
                    lines.append(f"    - {bucket}")
            else:
                lines.append(
                    f"  {method}: no reconciliation metadata published by this adapter "
                    "(requested vs effective calibration counts NOT verifiable from the "
                    "report; the adapter should expose metadata['calibration_reconciliation'])"
                )
        if not anything:
            lines.append("  (no calibration reconciliation information available)")
        return lines

    # ------------------------------------------------------------ 可视化

    def generate_visualization(
        self,
        results: Dict[str, DatingResult],
        tree: Any = None,
        comparison: Optional[_Comparison] = None,
    ) -> Optional[Path]:
        """生成比较柱状图：误差线来自真实区间，缺失值显式披露（C-15）。"""
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt
        import numpy as np

        self._configure_offline_fonts()

        comparison = comparison or self._align(results, tree)
        rows = comparison.comparable_rows
        methods = comparison.methods

        stats: Dict[str, Any] = {
            "rows": len(rows),
            "unresolved_rows": len(comparison.rows) - len(comparison.comparable_rows),
            "leaf_entries_dropped": len(comparison.leaf_entries),
            "bars_per_method": {},
            "error_bars_per_method": {},
            "missing_cells": 0,
            "cells_without_interval": 0,
            "figure": None,
        }
        self._plot_stats = stats

        if not rows:
            self.logger.info(
                "No internal/calibration nodes to visualize after filtering leaves"
            )
            return None

        figure, axes = plt.subplots(figsize=(max(12, 1.6 * len(rows)), 6))
        x = np.arange(len(rows))
        width = 0.8 / max(len(methods), 1)

        for i, method in enumerate(methods):
            means = []
            low_err: List[float] = []
            high_err: List[float] = []
            bars = 0
            with_error = 0
            for row in rows:
                estimate = row.values.get(method)
                if estimate is None:
                    # NaN 高度 → 该柱不画；数量在下面显式披露，不留"看起来像 0"
                    means.append(float("nan"))
                    low_err.append(0.0)
                    high_err.append(0.0)
                    stats["missing_cells"] += 1
                    continue
                bars += 1
                means.append(estimate.mean_age)
                if _has_interval(estimate):
                    lower, upper = _require_bounds(estimate)
                    low = max(0.0, estimate.mean_age - lower)
                    high = max(0.0, upper - estimate.mean_age)
                    low_err.append(low)
                    high_err.append(high)
                    with_error += 1
                else:
                    stats["cells_without_interval"] += 1
                    low_err.append(0.0)
                    high_err.append(0.0)
            stats["bars_per_method"][method] = bars
            stats["error_bars_per_method"][method] = with_error
            ci_labels = sorted(
                {
                    _ci_type_label(row.values[method])
                    for row in rows
                    if row.values.get(method) is not None
                }
            )
            suffix = (
                f" (bars {bars}/{len(rows)}"
                + (f", err±{'/'.join(ci_labels)} {with_error}" if with_error else "")
                + ")"
            )
            axes.bar(
                x + i * width,
                means,
                width,
                label=f"{method}{suffix}",
                yerr=np.array([low_err, high_err]),
                capsize=2.5,
                error_kw={
                    "ecolor": "black",
                    "elinewidth": 0.9,
                    "capthick": 0.9,
                    "alpha": 0.75,
                },
            )

        axes.set_xlabel("Node (aligned by descendant leaf set; Algn shown)")
        axes.set_ylabel("Age (Ma)")
        axes.set_title(
            "PhyloDater Comparison — error bars are the reported intervals "
            "(HPD95 / CI95 / RANGE are NOT equivalent, see *_ci_type in the TSV)"
        )
        axes.set_xticks(x + width * (len(methods) - 1) / 2)
        axes.set_xticklabels(
            [f"{row.display}\n[{row.alignment[:4]}]" for row in rows],
            rotation=45,
            ha="right",
        )
        axes.legend(fontsize=8)

        notes = (
            f"{stats['missing_cells']} method×node cell(s) have no estimate "
            f"(bar omitted, NOT age 0); {stats['cells_without_interval']} cell(s) "
            f"report no interval (plotted without error bar); "
            f"{stats['unresolved_rows']} unresolvable node key(s) omitted "
            f"(see 'uncomparable' in the tables); "
            f"{stats['leaf_entries_dropped']} leaf entry(ies) filtered."
        )
        figure.text(0.01, 0.012, "Disclosure: " + notes, fontsize=7)
        self.logger.info("Comparison chart disclosure: " + notes)
        if stats["missing_cells"] or stats["cells_without_interval"]:
            self.logger.warning(
                "Comparison chart has omitted bars / bars without error bars: " + notes
            )

        figure.tight_layout(rect=(0, 0.035, 1, 1))
        output_file = self.output_dir / "comparison_plot.png"
        figure.savefig(output_file, dpi=300)
        plt.close(figure)

        stats["figure"] = str(output_file)
        self.logger.info(f"Generated visualization: {output_file}")
        return output_file

    def _configure_offline_fonts(self) -> None:
        """配置离线模式字体支持

        在无网络环境下（如HPC集群计算节点），尝试使用系统可用字体，
        避免因缺少外部字体而崩溃。
        """
        # 禁用字体下载（防止尝试从网络获取字体）。
        # 'font.download' 只存在于 matplotlib < 3.8（3.8 起该机制被移除），
        # 而 matplotlib 的类型标注只承认"当前版本存在的键名"，走 RcParams 的
        # __setitem__ 会被判成非法键。这里改为先查键在不在、再通过 RcParams
        # 继承自 dict 的那套接口写入：键不存在时根本不写，也不再靠 except 兜。
        import matplotlib
        import matplotlib.font_manager as fm

        legacy_key = "font.download"
        if legacy_key in dict.keys(matplotlib.rcParams):
            dict.__setitem__(matplotlib.rcParams, legacy_key, False)

        available_fonts = [f.name for f in fm.fontManager.ttflist]

        preferred_fonts = [
            "DejaVu Sans",  # 大多数Linux系统默认
            "Liberation Sans",  # RHEL/CentOS
            "Arial",  # macOS/Windows
            "Helvetica",  # macOS
            "sans-serif",  # 通用回退
        ]

        selected_font = None
        for font in preferred_fonts:
            if font in available_fonts or font == "sans-serif":
                selected_font = font
                break

        if selected_font:
            matplotlib.rcParams["font.family"] = selected_font
            self.logger.debug(f"Using font: {selected_font}")

        # 处理Unicode字符（如希腊字母）
        matplotlib.rcParams["axes.unicode_minus"] = False
