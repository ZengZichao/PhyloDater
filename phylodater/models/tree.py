"""
系统发育树类

不可变值对象设计，多后端只读适配器（DendroPy / BioPython / ETE3-可选）。

## 树后端策略（审阅项 A-5 / B-10 的修复）

ete3 3.1.x 无条件导入标准库 ``cgi``，而 ``cgi`` 已在 Python 3.13 中被移除
（PEP 594），因此在 3.13+ 的解释器上 **ete3 根本无法导入**。此前本类的
``is_monophyletic`` / ``get_mrca`` / ``get_mrca_terminals_all`` /
``with_calibration_annotations`` / ``with_paml48_calibration_annotations``
全部直接依赖 ``self._ete3_cache``，没有任何回退路径，导致：

* ``pip install phylodater``（ete3 是 optional extra，默认不装）的用户在
  MCMCTree 的 ``prepare_inputs`` 上撞 ``ModuleNotFoundError: No module named 'cgi'``；
* 全部 5 个校准标注单元测试失败。

现在的分工是：

* **拓扑运算（MRCA、单系性、祖先判定、定根、校准标注写树）走本文件自带的
  零依赖纯 Python Newick 游走器**（``_NwNode`` / ``_NwTree``）——不要求任何
  第三方树库可导入；
* 距离矩阵等数值运算继续用 DendroPy（硬依赖，已在 ``pyproject.toml`` 声明）；
* ``_ete3_cache`` / ``_ete3_tree`` / ``_ete4_tree`` **仍然保留**：ete3 能导入时
  它们是可视化层需要的"ete 形状"树后端；不能导入时抛出带解释的 ``ImportError``
  （viz 层 :func:`phylodater.viz.tree_plot` 的后端解析逻辑正是按 ``ImportError``
  设计的），绝不让 ``import phylodater`` 因此失败。

本类不使用 frozen dataclass（需兼容遗留调用点以 _canonical_newick= 传入），
而是在实例上通过 functools.cached_property 缓存各后端树对象
（_nw_tree / _dendropy_tree / _biopython_tree / _ete3_cache），并用实例字段
_distance_matrix_cache 缓存距离矩阵，避免重复解析整棵树。
"""

import hashlib
import math
import re
from functools import cached_property
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)

# --------------------------------------------------------------------- #
# 第三方树后端的对象别名
#
# ete3 / Bio.Phylo / DendroPy 都不随包发布类型标注（py.typed），它们的
# 树/节点对象在本仓里只能按 ``Any`` 传递。下面三个别名把"真实类型来源"
# 写进名字与注释，而不是让读者误以 PhyloDater 自己的类型是 Any。
# --------------------------------------------------------------------- #
Ete3Tree = Any  # ete3.Tree 实例（无类型标注可用）
BioPhyloTree = Any  # Bio.Phylo.Newick.Tree 实例（无类型标注可用）
DendropyTree = Any  # dendropy.Tree 实例（无类型标注可用）

# 本模块的约束/校准类型。constraints.py 与 calibration.py 都不反向依赖本文件
# （calibration 只 import constraints），所以这里可以模块级导入而不成环；
# 旧代码里的函数内局部 import 只是当时的保守写法。
from .calibration import CalibrationPoint  # noqa: E402
from .constraints import AgeConstraint  # noqa: E402

# ======================================================================= #
# ete3 可用性：惰性探测（A-5/B-10）
# ======================================================================= #
#
# 模块导入时**不**尝试 import ete3：本机解释器上那次尝试会真的执行 ete3 的
# 包初始化（拉起 web 插件、matplotlib、tkinter 等），既慢又可能在别的环境里
# 抛非 ImportError 异常。这里用三态惰性标记，只有真正需要 ete 形状后端的
# 调用点（可视化层）才会触发探测。

_ETE_UNTESTED = 0
_ETE_AVAILABLE = 1
_ETE_UNAVAILABLE = 2

_ete_import_state: int = _ETE_UNTESTED
_ete_tree_class: Optional[Callable[..., Ete3Tree]] = None
_ete_import_error: Optional[BaseException] = None


def ete3_available() -> bool:
    """当前解释器能否导入 ete3（惰性探测，结果按进程缓存）。"""
    global _ete_import_state, _ete_tree_class, _ete_import_error
    if _ete_import_state == _ETE_UNTESTED:
        try:
            from ete3 import Tree as _Tree  # noqa: WPS433 (惰性可选依赖探测)

            _ete_tree_class = _Tree
            _ete_import_state = _ETE_AVAILABLE
        except BaseException as exc:  # pragma: no cover - 取决于解释器
            # ete3 的包初始化在 3.13+ 抛 ModuleNotFoundError（ImportError 子类），
            # 但第三方包的初始化也可能抛别的异常；这里一律按"不可用"处理，
            # 因为本项目的核心功能已经不再需要 ete3。
            _ete_import_error = exc
            _ete_import_state = _ETE_UNAVAILABLE
    return _ete_import_state == _ETE_AVAILABLE


def _require_ete3_tree_class() -> Callable[..., Ete3Tree]:
    """返回 ete3 的 ``Tree`` 类；不可用时抛出带解释的 ``ImportError``。"""
    tree_class = _ete_tree_class
    if not ete3_available() or tree_class is None:  # pragma: no branch - 取决于解释器
        cause = f" ({_ete_import_error})" if _ete_import_error else ""
        raise ImportError(
            "ete3 在当前解释器上不可用" + cause + "。ete3 3.1.x 会导入标准库 "
            "``cgi``，而 ``cgi`` 已在 Python 3.13 被移除（PEP 594），因此 ete3 只"
            "能在 Python <= 3.12 上使用。PhyloDater 的核心功能（MRCA、单系性、"
            "定根、MCMCTree 校准标注写树）已不依赖 ete3；只有可视化层需要 ete "
            "形状的树后端。请安装 ete4，或改用 "
            "``pip install 'phylodater[ete3]'``（需 Python <= 3.12）。"
        )
    return tree_class


# ======================================================================= #
# 零依赖纯 Python Newick 模型（A-5/B-10 的回退实现）
# ======================================================================= #


class _NwParseError(ValueError):
    """Newick 文本无法解析。"""


class _NwNode:
    """Newick 节点（可变，仅在本模块内部的重根过程中被构造/改写）。"""

    __slots__ = ("label", "children", "parent", "length", "_leafset")

    def __init__(
        self,
        label: str = "",
        length: Optional[float] = None,
        parent: Optional["_NwNode"] = None,
    ) -> None:
        self.label = label
        self.length = length
        self.parent = parent
        self.children: List["_NwNode"] = []
        self._leafset: Optional[frozenset] = None

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def iter_preorder(self) -> Iterator["_NwNode"]:
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def iter_postorder(self) -> Iterator["_NwNode"]:
        """后序（子先于父），迭代实现，避免深树上的递归溢出。"""
        out: List["_NwNode"] = []
        stack = [(self, False)]
        while stack:
            node, visited = stack.pop()
            if visited:
                out.append(node)
                continue
            stack.append((node, True))
            stack.extend((child, False) for child in reversed(node.children))
        return iter(out)

    def leaf_labels(self) -> Tuple[str, ...]:
        return tuple(
            node.label for node in self.iter_preorder() if node.is_leaf and node.label
        )

    def ancestor_chain(self) -> List["_NwNode"]:
        chain = [self]
        node = self
        while node.parent is not None:
            node = node.parent
            chain.append(node)
        return chain


_WHITESPACE = " \t\r\n\v\f"
_QUOTE_CHARS = "'\""


def _nw_skip_comment(text: str, index: int) -> int:
    """``text[index] == '['``：返回配对的 ``]`` 之后的下标（引号感知）。"""
    depth = 0
    length = len(text)
    quote: Optional[str] = None
    while index < length:
        ch = text[index]
        if quote is not None:
            if ch == quote:
                if index + 1 < length and text[index + 1] == quote:
                    index += 2
                    continue
                quote = None
        elif ch in _QUOTE_CHARS:
            quote = ch
        elif ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0:
                return index + 1
        index += 1
    return length


def _nw_read_label(text: str, index: int) -> Tuple[str, int]:
    """读取节点标签（支持 ``'...'`` / ``'"...'`` 引号与双写转义）。"""
    length = len(text)
    parts: List[str] = []
    while index < length and text[index] not in "(),;:[]":
        ch = text[index]
        if ch in _QUOTE_CHARS:
            quote = ch
            index += 1
            buf: List[str] = []
            while index < length:
                if text[index] == quote:
                    if index + 1 < length and text[index + 1] == quote:
                        buf.append(quote)
                        index += 2
                        continue
                    index += 1
                    break
                buf.append(text[index])
                index += 1
            parts.append("".join(buf))
            continue
        parts.append(ch)
        index += 1
    return "".join(parts).strip(), index


def _nw_read_length(text: str, index: int) -> Tuple[Optional[float], int]:
    """``text[index] == ':'``：读取分支长度。非法数值按"缺失"处理。"""
    length = len(text)
    end = index + 1
    while end < length and text[end] not in "(),;[]":
        end += 1
    token = text[index + 1 : end].strip()
    try:
        value = float(token)
    except ValueError:
        value = None
    if value is not None and not math.isfinite(value):
        # inf/nan 不是合法分支长度；交给 validate_branch_lengths 去披露。
        value = None
    return value, end


def _nw_parse_node(
    text: str, index: int, parent: Optional[_NwNode]
) -> Tuple[_NwNode, int]:
    length = len(text)
    while index < length and text[index] in _WHITESPACE:
        index += 1
    # 注释可以出现在一个节点的**起始**位置：BEAST / IQ-TREE / DendroPy 会把有根标记
    # （``[&R]`` / ``[&&R]`` / ``[&U]``）写在树字符串最前面，NHX 节点注释也可能紧跟在
    # ``(`` 之前。旧实现在这里只读标签，遇到 ``[`` 就停下，于是 ``[&R] (A:1,B:1);``
    # 被解析成一个无叶子的根，MCMCTree / treePL 通路与祖先关系校验直接报
    # "Cannot parse Newick for tree operations"。
    while index < length and text[index] == "[":
        index = _nw_skip_comment(text, index)
        while index < length and text[index] in _WHITESPACE:
            index += 1
    node = _NwNode(parent=parent)
    if index < length and text[index] == "(":
        index += 1
        while True:
            child, index = _nw_parse_node(text, index, node)
            node.children.append(child)
            while index < length and text[index] in _WHITESPACE:
                index += 1
            if index < length and text[index] == ",":
                index += 1
                continue
            if index < length and text[index] == ")":
                index += 1
                break
            raise _NwParseError(
                f"Malformed Newick near position {index}: {text[index:index+20]!r}"
            )
    node.label, index = _nw_read_label(text, index)
    while index < length:
        ch = text[index]
        if ch in _WHITESPACE:
            index += 1
            continue
        if ch == "[":
            index = _nw_skip_comment(text, index)
            continue
        if ch == ":":
            node.length, index = _nw_read_length(text, index)
            continue
        break
    return node, index


_PHYLIP_HEADER_RE = re.compile(r"^\s*\d+\s+\d+\s*$")

_LEADING_COMMENT_RE = re.compile(r"^\s*(?:\[[^\]]*\]\s*)+")


def strip_leading_newick_comments(newick: str) -> str:
    """去掉树字符串最前面的 ``[&R]`` / ``[&&U]`` / NHX 注释，返回可直接交给
    ete3 的写法。

    ete3 的新ick 解析器**不接受开头的有根标记注释**（BEAST、IQ-TREE、DendroPy、
    RAxML-NG 等主流工具输出的树常常以 ``[&R]`` 开头），传入会抛
    ``NewickError: Unexpected newick format``。本项目的零依赖游走器已能容忍该注释，
    但适配器在需要 ete3 时必须先用本函数预处理。

    只处理**开头**的注释：节点后面的 ``[&date=…]`` 等注释仍由各自的解析器处理。
    """
    if not newick:
        return newick or ""
    return _LEADING_COMMENT_RE.sub("", newick, count=1)


def _nw_parse(newick: Optional[str]) -> Optional[_NwNode]:
    """把 Newick 解析成 ``_NwNode`` 根节点；解析不出任何叶时返回 ``None``。

    容忍：PHYLIP 头（``"4 1"``）、NHX/``[&R]`` 注释、引号标签、缺失分支长度、
    多树文件（只取第一棵）。
    """
    if not newick or not newick.strip():
        return None
    text = newick.strip()
    # PAML 4.8 的树文件带 "n_species 1" 头部，必须跳过才能解析
    lines = text.split("\n")
    if len(lines) > 1 and _PHYLIP_HEADER_RE.match(lines[0]):
        text = "\n".join(lines[1:])
    try:
        root, _ = _nw_parse_node(text, 0, None)
    except (_NwParseError, RecursionError):
        return None
    if not root.leaf_labels():
        return None
    return root


def _nw_needs_quoting(label: str) -> bool:
    return label == "" or bool(re.search(r"[,;()\[\]:\"'\s]", label))


def _nw_quote(label: str) -> str:
    if _nw_needs_quoting(label):
        return "'" + label.replace("'", "''") + "'"
    return label


def _nw_fmt_length(value: float) -> str:
    """分支长度字符串：去掉浮点噪声，且绝不产出科学计数法。"""
    text = f"{value:.10g}"
    if "e" in text or "E" in text:
        text = f"{value:.12f}".rstrip("0")
        if text.endswith("."):
            text += "0"
    return text


def _nw_rooted_digest(node: _NwNode) -> str:
    """有根树的同构指纹：**兄弟顺序无关**，但对拓扑与枝长敏感。

    用于如实回答"重根到底改变了什么"：Newick 里两个子树谁写在前面无含义，
    因此不能拿字符串相等当判据；而根位置真的移动时（哪怕根分裂不变、只是沿着
    同一条枝挪了位置），枝长组合一定会变，这个指纹能抓住它。
    """
    digests: Dict[int, str] = {}
    for node in node.iter_postorder():
        children = "|".join(sorted(digests[id(child)] for child in node.children))
        length = "NA" if node.length is None else f"{round(node.length, 9):.9f}"
        digests[id(node)] = f"[{children}]{length}"
    return digests[id(node)]


def _nw_write(
    root: _NwNode,
    *,
    branch_lengths: bool = True,
    internal_labels: bool = True,
    annotations: Optional[Dict[int, str]] = None,
    terminate: bool = True,
) -> str:
    """序列化为一棵 Newick 串（迭代式后序，不依赖递归深度）。

    ``annotations`` 以 ``id(node)`` 为键：命中的内部节点在右括号后写入
    ``'<约束串>'``（PAML/MCMCTree 的校准标注约定，见审阅报告 §六.1）。
    """
    rendered: Dict[int, str] = {}
    for node in root.iter_postorder():
        if node.is_leaf:
            text = _nw_quote(node.label)
        else:
            inner = ",".join(rendered[id(child)] for child in node.children)
            text = f"({inner})"
            annotation = annotations.get(id(node)) if annotations else None
            if annotation:
                # 校准标注占据节点名的位置：MCMCTree 把引号内的 token 读进
                # .annotation 再按 B()/U()/L() 匹配，因此不再输出内部节点名。
                text += f"'{annotation}'"
            elif internal_labels and node.label:
                text += _nw_quote(node.label)
        if branch_lengths and node.length is not None:
            text += f":{_nw_fmt_length(node.length)}"
        rendered[id(node)] = text
    out = rendered[id(root)]
    return out + ";" if terminate else out


def _nw_adjacency(root: _NwNode) -> Dict[int, List[Tuple[_NwNode, Optional[float]]]]:
    """无根视角的邻接表：``id(node) -> [(邻居, 边长)]``（根自身的入枝被丢弃）。"""
    adj: Dict[int, List[Tuple[_NwNode, Optional[float]]]] = {}
    for node in root.iter_preorder():
        for child in node.children:
            adj.setdefault(id(node), []).append((child, child.length))
            adj.setdefault(id(child), []).append((node, child.length))
    return adj


def _nw_edge_length(value: Optional[float]) -> float:
    """缺失分支长度按 0 处理（与 ete3/dendropy 的 "0 枝长" 约定一致）。"""
    return 0.0 if value is None else float(value)


def _nw_expand(
    parent: _NwNode,
    start: _NwNode,
    edge_length: Optional[float],
    adjacency: Dict[int, List[Tuple[_NwNode, Optional[float]]]],
    visited: Set[int],
) -> _NwNode:
    """把 ``start`` 及其全部未访问邻居作为 ``parent`` 的子树挂上去（保持兄弟顺序）。"""
    visited.add(id(start))
    node = _NwNode(label=start.label, length=edge_length, parent=parent)
    parent.children.append(node)
    stack: List[Tuple[_NwNode, _NwNode]] = [(node, start)]
    while stack:
        current, origin = stack.pop()
        neighbours = adjacency.get(id(origin), [])
        pending = []
        for neighbour, length in neighbours:
            if id(neighbour) in visited:
                continue
            visited.add(id(neighbour))
            pending.append((neighbour, length))
        created = []
        for neighbour, length in pending:  # 前序入列，保持兄弟顺序
            child = _NwNode(label=neighbour.label, length=length, parent=current)
            current.children.append(child)
            created.append((child, neighbour))
        for item in reversed(created):  # 逆序入栈，弹出时仍是原顺序
            stack.append(item)
    return node


def _nw_suppress_unifurcations(root: _NwNode) -> _NwNode:
    """删除"只有一个孩子"的内部节点并把两段枝长相加（永不丢叶）。

    重根后原根（度 2）会退化成度 1 的链节点，必须合并；这与 ete3
    ``set_outgroup`` / dendropy ``suppress_unifurcations`` 的行为一致。
    """
    for node in root.iter_postorder():
        if node is root:
            continue
        if len(node.children) != 1:
            continue
        parent = node.parent
        if parent is None:
            continue
        child = node.children[0]
        merged = node.length
        if child.length is not None:
            merged = _nw_edge_length(merged) + child.length
        child.length = merged
        child.parent = parent
        try:
            position = parent.children.index(node)
        except ValueError:  # pragma: no cover - 后序遍历下不应发生
            continue
        parent.children[position] = child
    # 根本身退化为度 1 时（例如对单孩子树重根），把孩子提为根
    while len(root.children) == 1:
        root = root.children[0]
        root.parent = None
        root.length = None
    return root


def _nw_reroot_on_edge(
    root: _NwNode, child: _NwNode, distance_from_child: float
) -> _NwNode:
    """在 ``child`` 的入枝上、距该枝**子端** ``distance_from_child`` 处插入一个新的二元根。

    该枝原长 ``L`` 被劈成两段：靠 ``child`` 一侧长 ``x``，靠 ``child.parent`` 一侧长
    ``L - x``。新根恒有 2 个子节点，且**所有 tip-to-tip 距离保持不变**（这是重根
    操作唯一允许的不变量）。两个端点各有用途：

    * ``x = 0`` → 根贴在 ``child`` 上（"以 child 这个节点为根"，根到各叶的距离
      就等于从 child 量出的距离，最小方差定根要的正是这个位置）；
    * ``x = L`` → 根贴在 ``child.parent`` 上（"在 child 所在枝上、其父端定根"，
      若树本来就 rooted 在这一枝上则结果恒等，外群定根要的正是这个形态）；
    * ``0 < x < L`` → 把该枝劈成两段（中点定根落在枝内部时使用）。
    """
    parent = child.parent
    if parent is None:
        # 目标节点已经是根：重根即恒等变换（新根放在它与任一子节点之间不影响拓扑，
        # 直接复制一份返回，避免在"无入枝"处凭空造 0 长枝）。
        return _nw_copy_tree(root)
    total = _nw_edge_length(child.length)
    offset = min(max(distance_from_child, 0.0), total)
    adjacency = _nw_adjacency(root)
    new_root = _NwNode()
    visited: Set[int] = {id(child), id(parent)}
    _nw_expand(new_root, child, offset, adjacency, visited)
    _nw_expand(new_root, parent, total - offset, adjacency, visited)
    return _nw_suppress_unifurcations(new_root)


def _nw_copy_tree(root: _NwNode) -> _NwNode:
    """结构复制（枝长/标签/拓扑一致），用于"在原根上重根 = 原样返回"。"""
    adjacency = _nw_adjacency(root)
    new_root = _NwNode(label=root.label)
    visited: Set[int] = {id(root)}
    for neighbour, length in adjacency.get(id(root), []):
        _nw_expand(new_root, neighbour, length, adjacency, visited)
    new_root.length = root.length
    return new_root


def _nw_distances_from(
    root: _NwNode,
    start: _NwNode,
    adjacency: Optional[Dict[int, List[Tuple[_NwNode, Optional[float]]]]] = None,
) -> Tuple[Dict[int, float], Dict[int, Optional[_NwNode]]]:
    """树上两点路径唯一：一次 DFS 即得 ``start`` 到全部节点的距离与前驱。

    枝长为负时也给出"沿树路径的真实距离"（不像 Dijkstra 那样在负权上静默错乱）。
    ``adjacency`` 可传入预先建好的邻接表，让多候选扫描保持 O(n) 建表 + O(n)/候选。
    """
    if adjacency is None:
        adjacency = _nw_adjacency(root)
    distances: Dict[int, float] = {id(start): 0.0}
    previous: Dict[int, Optional[_NwNode]] = {id(start): None}
    stack: List[_NwNode] = [start]
    while stack:
        node = stack.pop()
        base = distances[id(node)]
        for neighbour, length in adjacency.get(id(node), []):
            if id(neighbour) in distances:
                continue
            distances[id(neighbour)] = base + _nw_edge_length(length)
            previous[id(neighbour)] = node
            stack.append(neighbour)
    return distances, previous


def _nw_find_node_by_id(root: _NwNode, node_id: int) -> Optional[_NwNode]:
    for node in root.iter_preorder():
        if id(node) == node_id:
            return node
    return None


def _nw_diameter_tips(root: _NwNode) -> Optional[Tuple[_NwNode, _NwNode, float]]:
    """树中最长 tip-to-tip 路径（"直径"，中点定根的依据）。

    树上"离任一节点最远的点必是直径端点"，因此两次 DFS 即可拿到直径两端。
    """
    leaves = [node for node in root.iter_preorder() if node.is_leaf]
    if len(leaves) < 2:
        return None
    adjacency = _nw_adjacency(root)
    distances, _ = _nw_distances_from(root, leaves[0], adjacency)
    first_id = max(
        (id(node) for node in leaves), key=lambda nid: distances.get(nid, float("-inf"))
    )
    first = _nw_find_node_by_id(root, first_id)
    if first is None:  # pragma: no cover
        return None
    distances, _ = _nw_distances_from(root, first, adjacency)
    last_id = max(
        (id(node) for node in leaves), key=lambda nid: distances.get(nid, float("-inf"))
    )
    total = distances.get(last_id)
    last = _nw_find_node_by_id(root, last_id)
    if total is None or first is None or last is None:  # pragma: no cover
        return None
    return first, last, total


def _nw_path(from_node: _NwNode, to_node: _NwNode) -> List[_NwNode]:
    """树上的唯一路径（含两端）。"""
    parent_of: Dict[int, _NwNode] = {}
    seen: Set[int] = {id(from_node)}
    stack: List[_NwNode] = [from_node]
    while stack:
        node = stack.pop()
        if id(node) == id(to_node):
            break
        for child in node.children:
            if id(child) not in seen:
                seen.add(id(child))
                parent_of[id(child)] = node
                stack.append(child)
        if node.parent is not None and id(node.parent) not in seen:
            seen.add(id(node.parent))
            parent_of[id(node.parent)] = node
            stack.append(node.parent)
    if id(to_node) not in seen:  # pragma: no cover - 同一棵树内不会发生
        return []
    path = [to_node]
    while path[-1] is not from_node:
        path.append(parent_of[id(path[-1])])
    path.reverse()
    return path


class _NwTree:
    """只读拓扑视图：叶索引、MRCA、后代叶集合、重根。"""

    def __init__(self, root: _NwNode) -> None:
        self.root = root
        self._leaf_index: Dict[str, _NwNode] = {}
        self._duplicate_labels: Set[str] = set()
        for node in root.iter_preorder():
            if not node.is_leaf or not node.label:
                continue
            if node.label in self._leaf_index:
                self._duplicate_labels.add(node.label)
                continue
            self._leaf_index[node.label] = node

    # ---- 基本信息 ---------------------------------------------------- #

    @property
    def ok(self) -> bool:
        return bool(self._leaf_index)

    def tip_labels(self) -> List[str]:
        return [
            node.label
            for node in self.root.iter_preorder()
            if node.is_leaf and node.label
        ]

    def leaf_nodes(self) -> List[_NwNode]:
        return [node for node in self.root.iter_preorder() if node.is_leaf]

    def node_leaf_set(self, node: _NwNode) -> frozenset:
        """``node`` 的**全部**后代叶名（缓存于节点自身）。"""
        if node._leafset is None:
            node._leafset = frozenset(node.leaf_labels())
        return node._leafset

    # ---- 拓扑查询 ---------------------------------------------------- #

    def lookup_leaf(self, name: Optional[str]) -> Optional[_NwNode]:
        if name is None:
            return None
        return self._leaf_index.get(str(name))

    def mrca(self, names: Sequence) -> Optional[_NwNode]:
        """``names``（叶名）的最近共同祖先；任一名不在叶集合中时返回 ``None``。"""
        wanted = [str(name) for name in names if name is not None and str(name) != ""]
        if not wanted:
            return None
        nodes = []
        for name in wanted:
            leaf = self.lookup_leaf(name)
            if leaf is None:
                return None
            nodes.append(leaf)
        chains = [node.ancestor_chain() for node in nodes]
        common: Set[int] = {id(node) for node in chains[0]}
        for chain in chains[1:]:
            common &= {id(node) for node in chain}
        for node in chains[0]:  # 从最深的一条链往上找，第一个公共祖先即 MRCA
            if id(node) in common:
                return node
        return None  # pragma: no cover - 至少共享根节点

    def mrca_of_nodes(self, names: Sequence) -> Optional[frozenset]:
        node = self.mrca(names)
        return None if node is None else self.node_leaf_set(node)

    # ---- 重根（返回一棵**新**树，本对象不被改写） --------------------- #

    def rooted_with_outgroup(
        self, names: Sequence
    ) -> Optional[Tuple[_NwNode, frozenset]]:
        """以外群所在枝为根：返回 ``(新根, 外群叶集合)``；外群无法定位时返回 ``None``。

        新根放在"外群枝"上、紧贴该枝的父端，因此外群整体作为根的两个子节点之一
        并保留其原有枝长；若树本来就 rooted 在这一枝上，结果是恒等的（所有枝长不变）。
        """
        node = self.mrca(names)
        if node is None:
            return None
        leaves = self.node_leaf_set(node)
        if node is self.root or node.parent is None:
            # 外群覆盖了整棵树：无法据此定根（否则会得到单孩子根）
            return None
        new_root = _nw_reroot_on_edge(self.root, node, _nw_edge_length(node.length))
        return new_root, leaves

    def midpoint_rooted(self) -> Optional[_NwNode]:
        """中点定根：最长 tip-to-tip 路径（树的"直径"）的中点置根，必要时劈开枝条。

        沿直径从一端走向另一端，累计距离首次达到 ``total/2`` 的那条枝即为新根所在枝；
        中点恰好落在节点上时，该节点成为根的一侧端点（枝长被夹到 ``[0, L]``），
        因此根**恒为二元根**，且所有 tip-to-tip 距离保持不变。
        """
        diameter = _nw_diameter_tips(self.root)
        if diameter is None:
            return None
        start, end, total = diameter
        if not math.isfinite(total) or total <= 0.0:
            return None
        target = total / 2.0
        path = _nw_path(start, end)
        if len(path) < 2:
            return None
        tolerance = 1e-9 * max(1.0, abs(total))
        travelled = 0.0
        for index in range(len(path) - 1):
            near, far = path[index], path[index + 1]
            if far.parent is near:
                # 顺流而下：far 是该枝的子端，走在 travelled + edge 处
                child, edge, child_at_far_end = far, _nw_edge_length(far.length), True
            elif near.parent is far:
                # 逆流而上：near 才是该枝的子端，走在 travelled 处
                child, edge, child_at_far_end = (
                    near,
                    _nw_edge_length(near.length),
                    False,
                )
            else:  # pragma: no cover - _nw_path 保证相邻节点间必有亲缘关系
                return None
            if target <= travelled + edge + tolerance:
                from_near = min(max(target - travelled, 0.0), edge)
                offset = edge - from_near if child_at_far_end else from_near
                return _nw_reroot_on_edge(self.root, child, offset)
            travelled += edge
        last = path[-1]
        if last.parent is None:  # pragma: no cover - 直径终点不会是根则不触发
            return _nw_copy_tree(self.root)
        return _nw_reroot_on_edge(self.root, last, _nw_edge_length(last.length))

    def min_variance_rooted(self) -> Optional[Tuple[_NwNode, float]]:
        """根到叶距离方差最小的位置定根（候选为全部现有节点，与原 ete3 实现同口径）。

        判据是"以该节点为根时根到各叶距离的总体方差"，即 ``numpy.var`` 的口径；
        在树上两点路径唯一，因此一次 DFS 就能拿到某个候选位置到全部叶的距离，
        邻接表只建一次（原 ete3 实现对每个候选都就地改拓扑并全树重算）。
        """
        leaves = self.leaf_nodes()
        if len(leaves) < 2:
            return None
        adjacency = _nw_adjacency(self.root)
        best_node: Optional[_NwNode] = None
        best_variance = float("inf")
        for node in self.root.iter_preorder():
            distances, _ = _nw_distances_from(
                root=self.root, start=node, adjacency=adjacency
            )
            values = [distances.get(id(leaf)) for leaf in leaves]
            if any(
                value is None for value in values
            ):  # pragma: no cover - 连通树不会发生
                continue
            variance = _nw_population_variance(
                [float(value) for value in values if value is not None]
            )
            if variance < best_variance:
                best_variance = variance
                best_node = node
        if best_node is None:  # pragma: no cover
            return None
        # x = 0：新根贴在该候选节点上，于是"根到各叶的距离"就是上面算出的那组距离。
        return (
            _nw_reroot_on_edge(self.root, best_node, 0.0),
            best_variance,
        )


def _nw_population_variance(values: Sequence[float]) -> float:
    """总体方差（与 ``numpy.var`` 默认口径一致），两遍算法以保证数值稳定。"""
    count = len(values)
    if count == 0:
        return 0.0
    mean = math.fsum(values) / count
    return math.fsum((value - mean) ** 2 for value in values) / count


# ======================================================================= #
# 校准约束串：一律经 constraints.py 的公开渲染器出口（审阅项 B-1/A-3/#4）
# ======================================================================= #

_MCMCTREE_B_RE = re.compile(
    r"^B\(\s*(?P<lo>[+-]?[0-9][0-9.eE+-]*)\s*,\s*(?P<hi>[+-]?[0-9][0-9.eE+-]*)\s*(?:,|$)"
)


def _mcmctree_ga_interval(
    constraint: Optional[AgeConstraint],
) -> Optional[Tuple[float, float]]:
    """从约束的 MCMCTree 渲染串中取出 ``(下界 Ga, 上界 Ga)``。

    不重复实现 Ma→Ga 换算与格式化：直接解析 :meth:`AgeConstraint.to_mcmctree_calib_string`
    的产物。那才是"MCMCTree 实际会拿到的界"，也是同节点校准求交时唯一正确的输入。
    非 ``B()`` 形态（``U()``/``L()``/``G()``/``SN()``/``ST()``，即单边或密度先验）
    无法安全求交，返回 ``None`` 由调用方显式报错。
    """
    if constraint is None:
        return None
    try:
        rendered = constraint.to_mcmctree_calib_string()
    except Exception:
        return None
    if not rendered:
        return None
    match = _MCMCTREE_B_RE.match(rendered.strip())
    if not match:
        return None
    try:
        lower = float(match.group("lo"))
        upper = float(match.group("hi"))
    except ValueError:  # pragma: no cover - 正则已保证是数字
        return None
    if not math.isfinite(lower) or not math.isfinite(upper) or upper <= lower:
        return None
    return lower, upper


def _render_merged_ga_interval(lower_ga: float, upper_ga: float) -> str:
    """把求交后的 Ga 区间交回 constraints.py 渲染成 MCMCTree ``B()`` 串。"""
    from .constraints import UniformAgeConstraint

    return UniformAgeConstraint(
        min_age=lower_ga * 1000.0, max_age=upper_ga * 1000.0
    ).to_mcmctree_calib_string()


def _render_ga_value(age_ma: float, source: str) -> str:
    """Ma → Ga 字符串（含地球年龄上限校验），复用 constraints.py 的公开渲染路径。

    此前这里是手写 ``f"{age_ma/1000:.4f}"``：任何 < 5 Ma 的年龄都会被压成
    ``0.0000``（PAML 读到零宽度/零上界）。constraints.py 的 ``_format_ga`` 会在
    4 位小数不够时自动扩位，而它唯一的公开出口是各约束的
    ``to_mcmctree_calib_string()``，故此处借 ``U()``/``L()`` 渲染器取同一个数。
    """
    from .constraints import MaximumAgeConstraint

    if not math.isfinite(age_ma):
        raise ValueError(f"{source}: age ({age_ma} Ma) is not a finite number")
    rendered = MaximumAgeConstraint(max_age=float(age_ma)).to_mcmctree_calib_string()
    match = re.match(r"^U\(\s*([0-9][0-9.eE+-]*)", rendered.strip())
    if not match:  # pragma: no cover - U() 渲染器形态固定
        raise ValueError(f"{source}: unexpected MCMCTree rendering {rendered!r}")
    return match.group(1)


def _constraint_to_b_string(
    constraint: Optional[AgeConstraint], prefix: str = "B("
) -> str:
    """把约束转成 PAML/MCMCTree 串（Ga 单位）。

    ``B(`` → 直接委托 constraints.py 的集中渲染器 :meth:`to_mcmctree_calib_string`，
    因此 14 种约束类型（含 SoftBounds/Maximum/Gamma/SkewNormal/SkewT）全部可用，
    且尾部概率槽位数与 PAML 4.10 的上游槽位一致（审阅项 B-1）。

    ``<`` / ``>`` → PAML 4.8 根节点标注形态（``'<value'`` / ``'>value'``）：
    只借 constraints.py 的 Ga 渲染取数，不在此重复 ``.4f`` 换算（审阅项 A-3）。
    """
    if constraint is None:
        return ""

    if prefix == "<":
        upper_ma = _constraint_upper_bound_ma(constraint)
        if upper_ma is None:
            return ""
        return f"<{_render_ga_value(upper_ma, 'root calibration')}"
    if prefix == ">":
        lower_ma = _constraint_lower_bound_ma(constraint)
        if lower_ma is None:
            return ""
        return f">{_render_ga_value(lower_ma, 'root calibration')}"

    return constraint.to_mcmctree_calib_string()


def _constraint_upper_bound_ma(
    constraint: Optional[AgeConstraint],
) -> Optional[float]:
    """约束显式上界（Ma）；无上界语义（如软下界、密度先验）时返回 ``None``。"""
    for attribute in ("max_age", "fixed_age"):
        value = getattr(constraint, attribute, None)
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if math.isfinite(float(value)):
                return float(value)
    return None


def _constraint_lower_bound_ma(
    constraint: Optional[AgeConstraint],
) -> Optional[float]:
    """约束显式下界（Ma）；无下界语义（如上限约束）时返回 ``None``。"""
    for attribute in ("min_age", "fixed_age"):
        value = getattr(constraint, attribute, None)
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if math.isfinite(float(value)):
                return float(value)
    return None


class PhylogeneticTree:
    """系统发育树 (不可变值对象)"""

    def __init__(self, newick: Optional[str] = None, **_kw: Any) -> None:
        # Public attribute (ETE3/adapters unified API)
        # 兼容遗留调用点以 _canonical_newick= 传入
        resolved_newick: str = (
            newick if newick is not None else _kw.get("_canonical_newick", "")
        )
        self.newick = resolved_newick.strip()
        # 实例级缓存：距离矩阵（树不可变，解析一次后复用，避免 O(n³) 重复计算）
        self._distance_matrix_cache: Optional[Dict] = None

    # ------------------------------------------------------------------ #
    # 后端
    # ------------------------------------------------------------------ #

    @cached_property
    def _nw_tree(self) -> Optional[_NwTree]:
        """本类的**默认拓扑后端**：零依赖纯 Python Newick 游走器（A-5/B-10）。

        不可变：需要改写拓扑（重根）的调用点一律走 :meth:`_fresh_nw_root`，
        以免污染这份共享缓存。
        """
        root = _nw_parse(self.newick)
        if root is None:
            return None
        tree = _NwTree(root)
        return tree if tree.ok else None

    def _fresh_nw_root(self) -> _NwNode:
        """重新解析一棵**可安全改写**的树（重根过程会就地构造新拓扑）。"""
        root = _nw_parse(self.newick)
        if root is None:
            from ..core.exceptions import TreeValidationError

            raise TreeValidationError(
                f"Cannot parse Newick for tree operations: {self.newick[:120]!r}"
            )
        return root

    # ------------------------------------------------------------------ #
    # 重根原语（供 services/tree_rooting.py 使用；不依赖 ete3）
    #
    # 约定：外群定根把新根放在"外群所在枝"上、紧贴外群一侧，于是外群仍是根的两个
    # 子节点之一并保留其原枝长（等价于 ete3 ``set_outgroup``）；中点定根把新根放在
    # 最长 tip-to-tip 路径的正中（必要时劈枝）；最小方差定根在**现有节点**中寻找
    # 根到叶距离方差最小的位置（候选含当前根，与原 ete3 实现口径一致）。
    # 三者都会合并重根后原根退化出的度-1 链节点，并且**永不增删叶节点**。
    # ------------------------------------------------------------------ #

    def _rooted_digest(self) -> str:
        """本树的有根同构指纹（兄弟顺序无关、对枝长敏感）；无法解析时返回 ``""``。

        定根服务用它判断"用户请求的重根到底有没有改变树"——Newick 字符串相等
        不是可用的判据（两棵子树写在根的哪一侧没有含义）。
        """
        root = _nw_parse(self.newick)
        return "" if root is None else _nw_rooted_digest(root)

    def _rerooted(self, new_root: _NwNode) -> "PhylogeneticTree":
        """序列化重根结果，并核验"重根没有增删任何叶节点"。

        重根是纯拓扑操作：叶集合必须逐字不变。任何偏差都意味着实现错误（这类
        "输出树看上去完全正常、实际少了一个分类单元"的失效形态正是审阅项 B-21
        警告的那一类），因此宁可报错也不写出一棵悄悄丢叶的树。
        """
        from ..core.exceptions import TreeValidationError

        before = sorted(self.tip_names)
        after = sorted(
            node.label
            for node in new_root.iter_preorder()
            if node.is_leaf and node.label
        )
        if before != after:
            missing = sorted(set(before) - set(after))
            extra = sorted(set(after) - set(before))
            raise TreeValidationError(
                f"重根后叶集合发生变化（缺失 {missing[:5]}，多出 {extra[:5]}），"
                f"拒绝输出丢叶树。输入树：{self.newick[:120]!r}"
            )
        return PhylogeneticTree(_nw_write(new_root))

    def rerooted_with_outgroup(
        self, outgroup_taxa: Sequence
    ) -> Tuple[Optional["PhylogeneticTree"], frozenset]:
        """以外群定根，返回 ``(新树, 外群后代叶集合)``；无法定根时返回 ``(None, 空集)``。

        第一个元素是 Optional：无法定根是预期内的业务分支（外群覆盖了整棵树），
        调用方 ``tree_rooting`` 已经按此写入跳过理由。旧标注写成必定返回一棵树，
        与本文档字符串里的 ``(None, 空集)`` 直接矛盾。
        """
        tree = _NwTree(self._fresh_nw_root())
        result = tree.rooted_with_outgroup(outgroup_taxa)
        if result is None:
            return None, frozenset()
        new_root, leaves = result
        return self._rerooted(new_root), leaves

    def rerooted_at_midpoint(self) -> Optional["PhylogeneticTree"]:
        """中点定根；无法确定直径（叶少于 2 或全部枝长为 0）时返回 ``None``。"""
        tree = _NwTree(self._fresh_nw_root())
        new_root = tree.midpoint_rooted()
        return None if new_root is None else self._rerooted(new_root)

    def rerooted_at_min_variance(self) -> Optional[Tuple["PhylogeneticTree", float]]:
        """最小方差定根，返回 ``(新树, 该位置的根到叶距离方差)``；无法计算时 ``None``。"""
        tree = _NwTree(self._fresh_nw_root())
        result = tree.min_variance_rooted()
        if result is None:
            return None
        new_root, variance = result
        return self._rerooted(new_root), variance

    @cached_property
    def _ete3_cache(self) -> Ete3Tree:
        """ETE3 树缓存 (向后兼容旧版接口)。

        使用 format=1(保留内部节点名)与 quoted_node_names=True,以支持
        含 ``.`` ``-`` 等特殊字符的 GTDB 嵌入式节点名。

        通过 functools.cached_property 真正缓存解析结果（newick 不可变，
        重复访问不再每次重新解析整棵树；此前每次访问都重建 Tree 属于性能 bug）。

        **仅供可视化层取"ete 形状"的树后端**：本类的拓扑/标注运算已不再依赖
        ete3（审阅项 A-5/B-10）。ete3 不可导入（Python >= 3.13 移除了 ``cgi``）
        时抛 :class:`ImportError`，viz 层正是按 ``ImportError`` 设计的降级路径。
        """
        return _require_ete3_tree_class()(self.newick, format=1, quoted_node_names=True)

    @property
    def _ete3_tree(self) -> Ete3Tree:
        """ETE3 树 (viz 模块期望的属性名,等同于 _ete3_cache)"""
        return self._ete3_cache

    @property
    def _ete4_tree(self) -> Ete3Tree:
        """ETE4 树 (viz 模块期望的属性名,通过 ete3 兼容返回)"""
        return self._ete3_cache

    @cached_property
    def _biopython_tree(self) -> BioPhyloTree:
        """BioPython 树（缓存解析结果）"""
        from io import StringIO

        from Bio import Phylo as BioPhylo

        return BioPhylo.read(StringIO(self.newick), "newick")

    @cached_property
    def _dendropy_tree(self) -> DendropyTree:
        """DendroPy 树（数值后端的解析结果，缓存）"""
        from dendropy import Tree

        return Tree.get(
            data=self.newick,
            schema="newick",
            preserve_underscores=True,
            rooting="force-rooted",
        )

    # ------------------------------------------------------------------ #
    # 读写
    # ------------------------------------------------------------------ #

    def write(self, file_path: Union[str, Path], format: str = "newick") -> None:
        """写入树到文件。

        Args:
            file_path: 输出路径
            format: 输出格式，支持 "newick"（默认）、"nhx"（与 newick 等价，
                NHX 注解已内嵌于 newick 字符串中）、"nexus"。对于 "nexus"，
                使用 BioPython 做最佳努力转换；转换失败则回退原始 newick 并告警，
                避免扩展名与内容不符（此前该参数被静默忽略）。
        """
        fmt = (format or "newick").lower()
        text = self.newick
        if fmt == "nexus":
            try:
                from io import StringIO

                from Bio import Phylo as BioPhylo

                tree = BioPhylo.read(StringIO(self.newick), "newick")
                buf = StringIO()
                BioPhylo.write(tree, buf, "nexus")
                text = buf.getvalue()
            except Exception as e:
                from ..infrastructure.logging import get_logger

                get_logger().warning(
                    f"Failed to convert tree to NEXUS, writing newick: {e}"
                )
                text = self.newick

        from ..infrastructure.safe_io import safe_writer

        with safe_writer(file_path, encoding="utf-8", newline="\n") as f:
            f.write(text)
            if not text.endswith("\n"):
                f.write("\n")

    # ------------------------------------------------------------------ #
    # 基本信息
    # ------------------------------------------------------------------ #

    @property
    def num_tips(self) -> int:
        """叶节点数量（DendroPy 优先，失败时用零依赖游走器计数）。"""
        try:
            return len(self._dendropy_tree.leaf_nodes())
        except Exception:
            tree = self._nw_tree
            if tree is not None:
                return len(tree.tip_labels())
            # 两个后端都失败时的保守回退
            return self.newick.count(",") + 1

    @property
    def tip_names(self) -> List[str]:
        """叶节点名称列表（按 Newick 中出现的顺序）。"""
        try:
            tree = self._dendropy_tree
            return [str(leaf.taxon.label) for leaf in tree.leaf_node_iter()]
        except Exception:
            # DendroPy 不可用/解析失败时不再返回空清单（那会让上游把
            # "拿不到叶名"误判成"树上没有叶"），改用零依赖游走器。
            tree = self._nw_tree
            return list(tree.tip_labels()) if tree is not None else []

    def without_branch_lengths(self) -> "PhylogeneticTree":
        """返回无分支长度的 Newick 字符串"""
        import re

        return PhylogeneticTree(re.sub(r":[0-9.eE+-]+", "", self.newick))

    def without_internal_labels(self) -> "PhylogeneticTree":
        """返回无内部节点标签的树"""
        import re

        return PhylogeneticTree(re.sub(r"\)([A-Za-z_][A-Za-z0-9_]*)", ")", self.newick))

    def is_binary(self) -> bool:
        """检查是否为严格二叉树"""
        try:
            tree = self._dendropy_tree
            return all(
                len(node.child_nodes()) in (0, 2) for node in tree.postorder_node_iter()
            )
        except Exception:
            root = self._nw_tree
            if root is None:
                return True
            return all(
                len(node.children) in (0, 2) for node in root.root.iter_preorder()
            )

    def validate_branch_lengths(self) -> Tuple[bool, List[str], List[str]]:
        """验证分支长度是否有效（检测负值与解析错误）。

        返回: (is_valid, invalid_nodes, errors)
        """
        invalid_nodes = []
        errors = []
        try:
            tree = self._dendropy_tree
        except Exception as e:
            return False, [], [f"无法解析 Newick: {e}"]

        for node in tree.postorder_node_iter():
            length = node.edge.length
            if length is not None and length < 0:
                label = str(node.taxon) if node.taxon else (node.label or "internal")
                invalid_nodes.append(label)
                errors.append(f"负分支长度 {length} 在节点 {label}")

        return len(errors) == 0, invalid_nodes, errors

    def get_rooting_info(self) -> dict:
        """返回根节点信息（真正检测有根性，而非总是返回 True）。

        检测启发式：DendroPy 以 force-rooted 解析 Newick。若根节点恰好有 2 个
        子节点，视为已定根（真实根放置）；若根节点有 ≥3 个子节点（典型的三叉
        根），则极可能是未定根树（unrooted）用多元根（multifurcation）表示，
        应标记为未定根，交由 tree_validator 按 require_rooted 报错或告警——这是
        分子定年必需的校验（此前永远返回 is_rooted=True 使该分支成为死代码）。

        同时计算 likely_artifact / shortest_ratio 供 tree_validator 使用：
        - likely_artifact: 根分支极短（可能为人为定根假象）时为 True；
        - shortest_ratio: 最短内部分支长度 / 最长内部分支长度（∈(0,1]），
          用于量化"短枝假象"的严重程度（根分支极短则其值趋近于 0）。
        """
        try:
            tree = self._dendropy_tree
            root = tree.seed_node
            children = root.child_nodes()
            n_children = len(children)
            # 二叉根 = 已定根；三叉及以上根 = 未定根（常态表示）
            is_rooted = n_children == 2

            # 外类群检测：已定根且某一子树仅含单一叶节点
            has_outgroup = False
            outgroup: List[str] = []
            if is_rooted:
                for ch in children:
                    leaves = [
                        nd.taxon.label for nd in ch.leaf_iter() if nd.taxon is not None
                    ]
                    if len(leaves) == 1:
                        has_outgroup = True
                        outgroup = [leaves[0]]
                        break

            # 内部分支长度统计（含根分支），用于 artifact 判定
            internal_lengths: List[float] = []
            root_edge = root.edge.length
            if root_edge is not None:
                internal_lengths.append(root_edge)
            for nd in root.preorder_iter():
                if nd is root:
                    continue
                if not nd.is_leaf() and nd.edge.length is not None:
                    internal_lengths.append(nd.edge.length)

            if internal_lengths:
                shortest = min(internal_lengths)
                longest = max(internal_lengths)
                shortest_ratio = shortest / longest if longest > 0 else 0.0
            else:
                shortest_ratio = 0.0

            # 根分支极短（相对内部分支）视为人为定根假象
            likely_artifact = (
                is_rooted
                and root_edge is not None
                and (root_edge < 1e-6 or (longest > 0 and root_edge < 0.01 * longest))
            )

            return {
                "is_rooted": is_rooted,
                "root_children": n_children,
                "root_depth": 0,
                "has_outgroup": has_outgroup,
                "outgroup": outgroup,
                "likely_artifact": likely_artifact,
                "shortest_ratio": shortest_ratio,
            }
        except Exception:
            # 解析失败时保守回退：视为未定根，便于上层按 require_rooted 处理
            return {
                "is_rooted": False,
                "root_children": 0,
                "root_depth": 0,
                "has_outgroup": False,
                "outgroup": [],
                "likely_artifact": False,
                "shortest_ratio": 0.0,
            }

    @property
    def is_rooted(self) -> bool:
        """树是否已定根(便捷属性,源于 get_rooting_info)"""
        return bool(self.get_rooting_info().get("is_rooted", True))

    def get_taxonomy_filename(self, filepath: str) -> str:
        return Path(filepath).name

    def strip_annotations(self) -> "PhylogeneticTree":
        """返回去除 NHX 标注 [xxx] 的树"""
        import re

        cleaned = re.sub(r"\[.*?\]", "", self.newick)
        return PhylogeneticTree(cleaned)

    def with_renamed_leaves(self, mapping: Mapping[str, str]) -> "PhylogeneticTree":
        """批量替换叶节点名称，同时保留内部节点标注（如 MCMCTree 校准串）。

        旧实现直接对 Newick 字符串做 ``str.replace``，当叶名（如 ``B``）与内部节点
        的 MCMCTree 校准标注（如 ``'B(0.0100, ...)'``）中的字符冲突时，会破坏标注。
        此处改为在解析后的树结构上仅替换叶节点标签，再重新序列化。
        """
        root = _nw_parse(self.newick)
        if root is None:
            return self
        for node in root.iter_preorder():
            if node.is_leaf and node.label in mapping:
                node.label = mapping[node.label]
        return PhylogeneticTree(
            _nw_write(
                root,
                branch_lengths=True,
                internal_labels=True,
                annotations=None,
                terminate=True,
            )
        )

    @classmethod
    def from_newick(cls, newick: str) -> "PhylogeneticTree":
        """从 Newick 字符串创建"""
        return cls(newick)

    @classmethod
    def from_file(
        cls,
        file_path: Union[str, Path],
        format: str = "newick",
        **kw: Any,
    ) -> "PhylogeneticTree":
        """从文件加载树

        编码:优先 UTF-8,失败时回退 GBK/GB18030(兼容中文 Windows 环境下的文件)。
        空文件:内容为空或纯空白时 raise ValueError。
        """
        raw = Path(file_path).read_bytes()
        if not raw:
            raise ValueError(f"Tree file is empty: {file_path}")
        text = None
        for enc in ("utf-8", "gbk", "gb18030"):
            try:
                text = raw.decode(enc)
                break
            except UnicodeDecodeError:
                continue
        if text is None:
            raise ValueError(
                f"Cannot decode tree file (tried utf-8/gbk/gb18030): {file_path}"
            )
        text = text.strip()
        if not text:
            raise ValueError(f"Tree file is empty: {file_path}")
        return cls(text)

    # ------------------------------------------------------------------ #
    # 拓扑查询（零 ete3 依赖）
    # ------------------------------------------------------------------ #

    def is_monophyletic(self, tip_names: Sequence[str]) -> bool:
        """检查一组叶节点是否单系（真实拓扑判定，不依赖 ete3）。

        判定与 ete3 ``check_monophyly(..., unrooted=False)`` 一致：目标叶集合的
        MRCA 所包含的叶节点恰好等于目标集合（不多不少）时为单系。

        含 ``.`` ``-`` 等特殊字符的 GTDB 嵌入式节点名按精确叶名匹配，
        绝不做子串匹配（``"A"`` 不会命中 ``"AB"``）。
        """
        if not tip_names:
            return True
        target = {
            str(name) for name in tip_names if name is not None and str(name) != ""
        }
        if not target:
            return True
        tree = self._nw_tree
        if tree is None:
            # 解析失败时退化为精确叶名集合包含判定（保守返回）。
            return target.issubset(set(self.tip_names))
        node = tree.mrca(sorted(target))
        if node is None:
            return target.issubset(set(self.tip_names))
        return tree.node_leaf_set(node) == frozenset(target)

    def get_mrca(self, tip_names: Sequence[str]) -> Optional[List[str]]:
        """返回 MRCA 包含的叶节点集合（真实拓扑，不依赖 ete3）。

        异常时返回 None（而非把输入当 MRCA 返回），避免调用方（如
        CalibrationResolver._is_ancestor）在出错时得到错误结果。
        """
        if not tip_names:
            return []
        tree = self._nw_tree
        if tree is None:
            return None
        node = tree.mrca(list(tip_names))
        if node is None:
            return None
        return sorted(tree.node_leaf_set(node))

    def get_mrca_terminals(
        self, tip_names: Sequence[str]
    ) -> Tuple[Optional[str], Optional[str]]:
        """返回遗传距离最远的代表性叶节点对（系统发育多样性最大化）。

        基于 DendroPy 的真实距离矩阵，在候选叶节点中选择两两系统发育
        距离最大的一对，用于 MRCA 定位时覆盖最广的进化分支。

        入参为空时没有可锁定的代表叶，返回 ``(None, None)``（调用方
        CalibrationResolver 就是按"拿到两个名字才能写 MRCA"处理的）。
        """
        if not tip_names:
            return None, None
        if len(tip_names) <= 2:
            last = tip_names[-1]
            return (tip_names[0], last)
        try:
            distances = self.get_distance_matrix()
            candidates = list(tip_names)
            max_dist = -1.0
            best_pair = (candidates[0], candidates[1])
            for i, tip1 in enumerate(candidates):
                for tip2 in candidates[i + 1 :]:
                    dist = distances.get((tip1, tip2), distances.get((tip2, tip1), 0.0))
                    if dist > max_dist:
                        max_dist = dist
                        best_pair = (tip1, tip2)
            return best_pair
        except Exception:
            return tip_names[0], tip_names[-1]

    def get_mrca_terminals_all(self, tip_names: Sequence[str]) -> List[str]:
        """返回 MRCA 下的所有后代叶节点（真实拓扑，非简单排序回显）。

        此前该方法仅返回排序后的输入，导致 CalibrationResolver 的
        单系性检查（_resolve_with_monophyly_check）形同虚设。现改为基于
        零依赖 Newick 游走器计算 MRCA 的真实叶节点集合。解析失败时回退到
        排序后的输入。
        """
        if not tip_names:
            return []
        leaves = self.get_mrca(tip_names)
        if leaves is None:
            return sorted(tip_names)
        return leaves

    def get_distance_matrix(self) -> Dict[Tuple[str, str], float]:
        """返回所有叶节点对的系统发育距离矩阵。

        矩阵以 (tip1, tip2) -> distance 的字典形式返回，包含对称条目。
        解析失败时返回空字典。结果按实例缓存（树不可变），避免 O(n³) 重复计算
        （此前每次调用都重新解析整棵树并两两求距）。
        """
        if self._distance_matrix_cache is not None:
            return self._distance_matrix_cache
        matrix: Dict = {}
        try:
            tree = self._dendropy_tree
            pdm = tree.phylogenetic_distance_matrix()
            taxa = list(tree.taxon_namespace)
            for i, t1 in enumerate(taxa):
                for t2 in taxa[i + 1 :]:
                    try:
                        dist = float(pdm.distance(t1, t2))
                        matrix[(str(t1.label), str(t2.label))] = dist
                        matrix[(str(t2.label), str(t1.label))] = dist
                    except Exception:
                        continue
        except Exception:
            matrix = {}
        self._distance_matrix_cache = matrix
        return matrix

    # ------------------------------------------------------------------ #
    # 公开的"分支/单系群"枚举 API（跨模块需求：比较层可据此退役本地解析器）
    # ------------------------------------------------------------------ #

    def iter_clades(self) -> List[Dict]:
        """枚举树上**每一个节点**的结构事实（前序：父先于子，深度沿路径递增）。

        每个条目是 ``{"label", "leaves", "is_leaf", "depth", "topology_id"}``：

        * ``leaves``：该节点的**全部**后代叶名（``tuple``，Newick 出现顺序）；
        * ``depth``：该节点到根的**边数**（根为 0，叶节点为其祖先边数）；
        * ``topology_id``：``"<叶数>L:<sha1(排序叶集合)[:8]>"`` —— 同一棵树里
          唯一的拓扑身份，也可跨树比对（与 ``comparison_reporter`` 的口径一致）；
        * ``label``：原始节点名（叶节点即其叶名；无名内部节点为 ``""``）。

        不依赖 ete3/DendroPy。解析失败时返回 ``[]``。
        """
        tree = self._nw_tree
        if tree is None:
            return []
        root = tree.root
        depth_of: Dict[int, int] = {id(root): 0}
        entries: List[Dict] = []
        for node in root.iter_preorder():
            depth = depth_of[id(node)]
            for child in node.children:
                depth_of[id(child)] = depth + 1
            leaves = node.leaf_labels()
            entries.append(
                {
                    "label": node.label,
                    "leaves": leaves,
                    "is_leaf": node.is_leaf,
                    "depth": depth,
                    "topology_id": _topology_id(leaves),
                }
            )
        return entries

    def get_descendant_leaves(self, tip_names: Sequence[str]) -> Optional[List[str]]:
        """``tip_names`` 所锚定节点的**全部**后代叶名（排序）；无法定位时 ``None``。

        与 :meth:`get_mrca` 同义，但名字直接说明"返回的是后代叶集合"，供
        比较层/报告层按拓扑对齐节点时使用。
        """
        return self.get_mrca(tip_names)

    def get_node_topology_id(self, tip_names: Sequence[str]) -> str:
        """目标节点的稳定拓扑标识符：由其 MRCA 的**完整**后代叶集合导出。

        此前实现对入参 ``tip_names`` 排序拼接，因此"用哪两个代表类群去锚定"
        会决定身份，同一节点上的两条校准拿不到同一个 ID。现在锚定同一节点的
        不同叶对必然得到同一 ID（``"<叶数>L:<sha1-8>"``）。

        无法定位 MRCA 时退回输入自身的排序拼接（``"?"`` 前缀），使"没能按拓扑
        定位"这件事在结果里可见，而不是伪装成一个正常 ID。
        """
        names = [str(name) for name in (tip_names or []) if name]
        leaves = self.get_mrca(names) if names else None
        if leaves is None:
            return "?" + ",".join(sorted(names))
        return _topology_id(sorted(leaves))

    def get_node_depth(self, tip_names: Sequence[str]) -> int:
        """目标节点（MRCA）到根节点的深度（以边数计）。

        解析失败或无法定位 MRCA 时返回 0。
        """
        names = [str(name) for name in (tip_names or []) if name]
        if not names:
            return 0
        tree = self._nw_tree
        if tree is None:
            return 0
        node = tree.mrca(names)
        if node is None:
            return 0
        depth = 0
        current = node
        while current.parent is not None:
            depth += 1
            current = current.parent
        return depth

    def _get_tip_by_name(self, name: str) -> Optional[Any]:
        try:
            tree = self._dendropy_tree
            for leaf in tree.leaf_node_iter():
                if leaf.taxon and str(leaf.taxon.label) == name:
                    return leaf
            return None
        except Exception:
            nw = self._nw_tree
            return None if nw is None else nw.lookup_leaf(name)

    def _get_distance_matrix_biopython(self) -> Dict[Tuple[str, str], float]:
        return {}

    def iter_tip_pairs(self) -> Iterator[Tuple[str, str]]:
        names = self.tip_names
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                yield names[i], names[j]

    def get_distance_matrix_chunked(
        self, chunk_size: int = 1000
    ) -> Dict[Tuple[str, str], float]:
        return {}

    def to_nhx_string(self) -> str:
        """返回 NHX 格式字符串。

        注意：当前实现返回原始 Newick（不含 NHX 标注）。
        如需完整 NHX 标注支持，请使用 ete3 的 ``write(format=1, features=...)`` 方法。
        """
        return self.newick

    # ------------------------------------------------------------------ #
    # 校准标注写树（MCMCTree / PAML）
    # ------------------------------------------------------------------ #

    def with_calibration_annotations(
        self,
        calibrations: Sequence[CalibrationPoint],
        include_root: bool = False,
        default_root_constraint: Optional[AgeConstraint] = None,
    ) -> "PhylogeneticTree":
        """
        PAML 4.9+ / 4.10.8+ 格式: 返回无分支长度、带 MCMCTree 约束标注的 Newick。

        约束通过 `to_mcmctree_calib_string()` 写入内部节点名，单位转换为 Ga。
        默认情况下根节点约束由控制文件中的 `RootAge` 处理，因此不在树中标注。
        对 PAML 4.10+，``RootAge`` 已被弃用，可将 ``include_root=True`` 把根约束
        直接写入树根；未提供根校准时可用 ``default_root_constraint`` 写入宽松默认界。
        PHYLIP 头 "n_species 1" 由调用方（mcmctree_method._add_phylip_header_to_tree）
        单独追加，避免重复。
        """
        root = self._fresh_nw_root()
        tree = _NwTree(root)
        annotations = {
            id(node): text
            for node, text in self._resolve_annotations_on(tree, calibrations)
        }

        if include_root:
            root_cal = next(
                (cal for cal in calibrations if getattr(cal, "is_root_node", False)),
                None,
            )
            constraint = getattr(root_cal, "age_constraint", None)
            if constraint is None:
                constraint = default_root_constraint
            if constraint is not None:
                rendered = constraint.to_mcmctree_calib_string()
                if rendered:
                    annotations[id(tree.root)] = rendered

        return PhylogeneticTree(
            _nw_write(
                tree.root,
                branch_lengths=False,
                internal_labels=False,
                annotations=annotations,
                terminate=False,
            )
        )

    def _resolve_annotations_on(
        self, tree: _NwTree, calibrations: Sequence[CalibrationPoint]
    ) -> List[Tuple[_NwNode, str]]:
        """把校准点解析到树上的节点并渲染约束串（供两条写树路径共用）。

        返回 ``[(节点, 约束串)]``。两条曾经静默失败的路径（审阅项 B-6）现在是
        显式行为：

        1. 非根校准缺少可用的 ``mrca_leaf_pair``、或其叶名在树中定位不到、
           或定位到的是叶节点而非内部节点 → :class:`CalibrationResolutionError`
           （此前 ``continue`` 直接丢弃，MCMCTree 用比用户指定更少的校准跑完）；
        2. 多条校准解析到**同一个节点**（以"完整后代叶集合"判同，因此叶对不同的
           嵌套写法也会被抓住）→ 区间求交（下界取各下界的最大、上界取各上界
           的最小）并 ``logger.warning`` 披露；交集为空 →
           :class:`CalibrationConflictError`（此前用 ``id(mrca)`` 作字典键，后者
           覆盖前者，且不留任何痕迹）。
        """
        from ..core.exceptions import (
            CalibrationConflictError,
            CalibrationResolutionError,
        )
        from ..infrastructure.logging import get_logger

        logger = get_logger()
        grouped: Dict[Tuple, List] = {}
        order: List[Tuple] = []
        for cal in calibrations:
            if getattr(cal, "is_root_node", False):
                continue
            pair = getattr(cal, "mrca_leaf_pair", None)
            usable = tuple(pair)[:2] if pair else ()
            usable = tuple(
                str(name) for name in usable if name is not None and str(name) != ""
            )
            if len(usable) < 2:
                raise CalibrationResolutionError(
                    f"校准点 '{getattr(cal, 'name', '?')}' 是内部节点校准，却没有可用的 "
                    f"mrca_leaf_pair（当前为 {pair!r}），无法定位到树上的节点。"
                    f"写树时丢弃它会让 MCMCTree 用比用户指定更少的校准完成定年，"
                    f"因此此处直接报错。请先用 CalibrationResolver 解析出代表叶节点对，"
                    f"或在 YAML 中显式提供 mrca_pair。"
                )
            node = tree.mrca(list(usable))
            if node is None:
                raise CalibrationResolutionError(
                    f"校准点 '{getattr(cal, 'name', '?')}' 的 mrca_leaf_pair "
                    f"{usable} 无法在树中定位（叶名不存在或树已改动）。"
                    f"拒绝静默丢弃该校准点。"
                )
            leaves = tree.node_leaf_set(node)
            if len(leaves) < 2:
                raise CalibrationResolutionError(
                    f"校准点 '{getattr(cal, 'name', '?')}' 的 mrca_leaf_pair "
                    f"{usable} 定位到了叶节点 {sorted(leaves)} 而不是内部节点。"
                    f"MCMCTree 的 B() 标注只能挂在内部节点上；两个代表类群必须"
                    f"分属不同分支。"
                )
            key = tuple(sorted(leaves))
            if key not in grouped:
                grouped[key] = [node]
                order.append(key)
            grouped[key].append(cal)

        annotations: List[Tuple[_NwNode, str]] = []
        for key in order:
            node = grouped[key][0]
            entries = grouped[key][1:]
            if not entries:
                continue
            if len(entries) == 1:
                constraint = getattr(entries[0], "age_constraint", None)
                if constraint is None:
                    logger.warning(
                        f"Calibration '{getattr(entries[0], 'name', '?')}' has no "
                        f"age_constraint: nothing can be written for its node "
                        f"{list(key)[:3]}...; it is NOT silently dropped but simply "
                        f"carries no prior to render."
                    )
                    continue
                rendered = constraint.to_mcmctree_calib_string()
                if rendered:
                    annotations.append((node, rendered))
                continue

            intervals = []
            for cal in entries:
                interval = _mcmctree_ga_interval(getattr(cal, "age_constraint", None))
                if interval is None:
                    raise CalibrationConflictError(
                        f"节点（后代叶集合 {sorted(key)[:5]}...）上有多条校准"
                        f"（{[getattr(c, 'name', '?') for c in entries]}），其中至少一条"
                        f"不是可求交的 B(min,max) 区间，无法安全合并。拒绝按"
                        f"'后者覆盖前者'静默丢弃其余校准（审阅项 B-6）。"
                        f"请人工把它们合并为一条区间校准。"
                    )
                intervals.append(interval)
            lower = max(item[0] for item in intervals)
            upper = min(item[1] for item in intervals)
            names = [getattr(cal, "name", "?") for cal in entries]
            if not upper > lower:
                raise CalibrationConflictError(
                    f"节点（后代叶集合 {sorted(key)[:5]}...）上的 {len(entries)} 条校准"
                    f"（{names}）区间交集为空：合并后下界 {lower:.6g} Ga 不小于上界 "
                    f"{upper:.6g} Ga。这些先验互相矛盾，MCMCTree 会以 "
                    f"'fossil bounds in tree incorrect' 中止。请放宽其中一条。"
                )
            merged = _render_merged_ga_interval(lower, upper)
            logger.warning(
                f"{len(entries)} calibrations {names} resolve to the SAME node "
                f"(descendant leaves: {sorted(key)[:5]}"
                f"{'...' if len(key) > 5 else ''}); their MCMCTree intervals were "
                f"merged by intersection into {merged} (lower = max of lowers, "
                f"upper = min of uppers). Previously the last one silently "
                f"overwrote the others (review item B-6)."
            )
            annotations.append((node, merged))
        return annotations

    def with_paml48_calibration_annotations(
        self, calibrations: List
    ) -> "PhylogeneticTree":
        """
        PAML 4.8 专用: 返回满足 PAML 4.8 约束格式 + 无分支长度的树文件

        PAML 4.8 wiki 格式:
          ((human,chimp) 'B(14.8,15.2,0.01)', (dog,(mouse,rat) 'B(11.0,13.0)')'<60<');

        注: 约束以 Ga 为单位（1000 Ma = 1 Ga）。B() 的槽位数由 constraints.py 的
        集中渲染器决定（PAML 4.8/4.10 的 ``npfossils[BOUND_F]`` 均为 4 槽，见
        审阅项 B-1），根节点上界走 ``'<value'`` 后缀。
        """
        root = self._fresh_nw_root()
        tree = _NwTree(root)
        annotations = {
            id(node): text
            for node, text in self._resolve_annotations_on(tree, calibrations)
        }
        newick = _nw_write(
            tree.root,
            branch_lengths=False,
            internal_labels=False,
            annotations=annotations,
            terminate=False,
        )

        # 追加根节点约束 (如果有)：一律经 constraints.py 的 Ga 渲染出口
        root_cal = next(
            (cal for cal in calibrations if getattr(cal, "is_root_node", False)), None
        )
        if root_cal is not None:
            s = _constraint_to_b_string(
                getattr(root_cal, "age_constraint", None), prefix="<"
            )
            if s:
                newick += f"'{s}'"

        # PAML 4.8 要求树文件带头部: "n_species 1"
        # n_species = 叶节点数量
        n_leaves = len(tree.tip_labels())
        return PhylogeneticTree(f"{n_leaves} 1\n{newick}")


def _topology_id(leaves: Sequence[str]) -> str:
    """``"<叶数>L:<sha1-8 of 排序叶集合>"``：跨方法/跨树对齐节点的稳定身份。"""
    token = "{" + ",".join(sorted(leaves)) + "}"
    digest = hashlib.sha1(token.encode("utf-8")).hexdigest()[:8]
    return f"{len(set(leaves))}L:{digest}"
