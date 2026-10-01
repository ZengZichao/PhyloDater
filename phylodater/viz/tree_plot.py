"""
Tree plotting module for PhyloDater visualization.

Based on ChronoPhylo, adapted for PhyloDater's DatingResult.

时间轴约定（审阅项 B-13 / B-14 / C-13）
------------------------------------
所有节点（**包括叶节点**）的横坐标都是"距今 Ma"的年龄：

* 节点自带 ``age`` 注解时优先使用（定年结果、相对时间树都走这条）；
* 否则按同一套公式 ``age = root_age - dist(root, node)`` 计算——叶节点不再
  被无条件钉在 0（B-14：灭绝类群与带采样日期的尖端会被画到"现在"）；
* 叶节点若带有采样日期 / 化石年龄注解（``A|-0.045``、``[&date=-0.045]``、
  ``[&&NHX:age=55]``、``age``/``tip_age`` 特征），使用其真实年龄，并在日志里
  报告赋值数量；
* 距离计算失败一律上抛，不再静默画成"该节点发生在现在"（C-13）；只有"既没有
  年龄注解、到根路径上也没有任何枝长信息"时才退回 0 Ma，并且必定告警。

``root_age`` 缺省时按真实推断处理（B-13）：根节点注解 → root-to-tip 最长分支
路径 → 已有节点年龄注解的最大值；全部不可用时显式报错，绝不静默取 0。

本模块**不导入 ete3/ete4**（A-5/B-10）：树一律按鸭子类型访问，
:class:`~phylodater.models.tree.PhylogeneticTree` 则经 :func:`as_ete_tree`
走模型层暴露的后端入口。
"""

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection

from ..infrastructure.logging import get_logger

logger = get_logger()

#: root-to-tip 距离相对极差超过该值即判为"非超度量（不是时间树）"。
#: 口径与 ``services/deep_validator.ULTRAMETRIC_DEFAULT_TOLERANCE`` 一致（B-18），
#: 此处独立定义以避免 viz → services 的反向依赖。
ULTRAMETRIC_TOLERANCE = 1e-6

#: 显式承载"节点年龄（距今 Ma）"的树内特征键，按优先级排列。
NODE_AGE_FEATURE_KEYS: Tuple[str, ...] = (
    "age",
    "node_age",
    "tip_age",
    "sampling_age",
)

#: 显式承载"采样日期"的特征键。取值轴与 ``adapters/lsd2_method.py`` 一致：
#: 负值 = 过去（``-0.045`` 即 45 ka 前采样），单位 Ma。
NODE_DATE_FEATURE_KEYS: Tuple[str, ...] = (
    "date",
    "tip_date",
    "sampling_date",
    "collection_date",
    "collect_date",
)

#: 叶标签尾随 ``|`` 的采样日期后缀（BEAST / IQ-TREE ``--date TAXNAME`` 约定）
_LABEL_DATE_SUFFIX_RE = re.compile(r"\|([^|]*)$")
#: Newick 标签里的注解块：``A[&date=-0.045]`` / ``A[&&DATE=2020-03-01]`` /
#: ``A[&&NHX:age=55:tip_age=55]``（ETE 的 NHX 块用 ``:`` 分隔，phylocom 风格用 ``;``）
_ANNOTATION_BLOCK_RE = re.compile(r"\[(?:&&|&)?\s*(?:NHX:)?([^\]]*)\]", re.IGNORECASE)
_ANNOTATION_PAIR_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+?)\s*$")
_ANNOTATION_TOKEN_RE = re.compile(r"[;:]")
_NUMERIC_TOKEN_RE = re.compile(r"^[+-]?\d+(?:\.\d+)?$")
_BARE_YEAR_RE = re.compile(r"^[+-]?\d{4}$")

#: 逐节点告警的去重登记表：同一类问题只报一次（附上首个 offenders），
#: 避免上千个尖端把日志刷满。
_WARNED_KEYS: Dict[str, str] = {}


def _warn_once(key: str, message: str) -> None:
    """同一 ``key`` 的告警只发一次（消息里带首个触发对象名）。"""
    if key in _WARNED_KEYS:
        return
    _WARNED_KEYS[key] = message
    logger.warning(message)


@dataclass
class NodeCoordinate:
    """Data class for storing node coordinates."""

    x: float
    y: float
    name: str
    is_leaf: bool
    is_root: bool
    age: float
    props: Dict[str, Any]


# --------------------------------------------------------------------------- #
# ete3 / ete4 通用的小工具（不硬依赖 ete3：全部按鸭子类型访问）
# --------------------------------------------------------------------------- #


def is_ete4_node(node: Any) -> bool:
    """ete4 节点带 ``props`` 映射，ete3 节点没有（用 getattr 读特征）。"""
    return isinstance(getattr(node, "props", None), dict)


def get_tree_node_attr(node: Any, attr: str, default: Any = None) -> Any:
    """读取节点特征/注解（兼容 ete3 的 feature 与 ete4 的 props）。"""
    props = getattr(node, "props", None)
    if isinstance(props, dict):
        return props.get(attr, default)
    return getattr(node, attr, default)


def _is_leaf(node: Any) -> bool:
    """ete3 的 is_leaf 是方法，ete4 的是属性。"""
    checker = getattr(node, "is_leaf", None)
    if callable(checker):
        return bool(checker())
    if checker is not None:
        return bool(checker)
    return not list(getattr(node, "children", []) or [])


def _iter_nodes(root_node: Any) -> List[Any]:
    """后序/前序皆可：只需要遍历全树。"""
    traverse = getattr(root_node, "traverse", None)
    if callable(traverse):
        return list(traverse())
    nodes = [root_node]
    stack = list(getattr(root_node, "children", []) or [])
    while stack:
        current = stack.pop()
        nodes.append(current)
        stack.extend(getattr(current, "children", []) or [])
    return nodes


def _branch_length_to_parent(node: Any) -> Optional[float]:
    """节点到父节点枝长。返回 None 表示该节点根本没有枝长信息。"""
    raw = getattr(node, "dist", None)
    if raw is None:
        raw = getattr(node, "branch_length", None)
    if raw is None:
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        raise ValueError(
            f"节点 '{getattr(node, 'name', '') or '(未命名)'}' 的分支长度为 {raw!r}，"
            "非有限值无法用于定年绘图。"
        )
    return value


def _sum_up_chain(root: Any, node: Any) -> float:
    """沿 ``up`` 链将枝长累加到根；链在到达根之前断开则报错。"""
    total = 0.0
    seen = {id(node)}
    current = node
    while current is not root:
        step = _branch_length_to_parent(current)
        total += 0.0 if step is None else step
        parent = getattr(current, "up", None)
        if parent is None:
            raise ValueError(
                f"从节点 '{getattr(node, 'name', '') or '(未命名)'}' 回溯到根节点失败："
                "父链在到达根节点之前断开（树不完整或节点不属于该树）。"
            )
        if id(parent) in seen:
            raise ValueError("树中存在循环的父节点引用，无法计算根到叶距离。")
        seen.add(id(parent))
        current = parent
    return total


def distance_from_root(root: Any, node: Any) -> float:
    """根到 ``node`` 的分支长度距离。

    审阅项 C-13：此前任何失败都被 ``except Exception: return 0.0`` 画成
    "该节点发生在现在"，图上"年龄为 0"与"计算失败"长得一模一样。现在一律上抛，
    由调用方中止出图。

    Raises
    ------
    ValueError
        距离无法计算（后端 API 漂移、树不完整、枝长为 NaN/inf 等）。
    """
    failures: List[str] = []

    getter = getattr(root, "get_distance", None)
    if callable(getter):
        try:
            value = float(getter(node))
        except Exception as exc:  # 后端实现差异：记下原因后换一条路径
            failures.append(f"root.get_distance() 失败: {exc}")
        else:
            if math.isfinite(value):
                return value
            failures.append(f"root.get_distance() 返回非有限值 {value!r}")
    else:
        failures.append("该树后端不提供 get_distance()")

    try:
        value = _sum_up_chain(root, node)
    except ValueError as exc:
        failures.append(str(exc))
    else:
        if math.isfinite(value):
            return value
        failures.append(f"枝长累加得到非有限值 {value!r}")

    name = getattr(node, "name", "") or "(未命名)"
    raise ValueError(
        f"无法计算节点 '{name}' 到根节点的距离，因此其年龄无法确定；"
        "可视化层不会把它画成'现在'。原因: " + "; ".join(failures)
    )


# --------------------------------------------------------------------------- #
# 尖端 / 节点年龄注解（B-14）
# --------------------------------------------------------------------------- #


def _sampling_axis_to_age_ma(token: Any, *, where: str) -> Optional[float]:
    """把"采样日期"取值解释为距今 Ma。

    沿用 PhyloDater 既有轴向约定（见 ``adapters/lsd2_method.py`` 的
    ``_extract_tip_sampling_dates``）：负值 = 过去，单位为 Ma。
    刻意**不做**单位猜测：日历式日期（``2019-12-31``）与裸四位数（``2019``，
    既可能是公元年份也可能是 2019 Ma）一律不使用，只告警。
    """
    text = str(token).strip().strip('"').strip("'")
    if not text:
        return None

    value = _parse_finite_age(text) if _NUMERIC_TOKEN_RE.match(text) else None
    if value is None:
        _warn_once(
            "sampling-date-not-ma-axis",
            f"{where}: 采样日期 '{text}' 无法解释为距今 Ma（日历式日期需要先确定"
            "'现在'是哪一年，可视化层不做该假设），改用分支长度推算该节点年龄。",
        )
        return None

    if value < 0:
        return -value  # 负值 = 过去（PhyloDater / LSD2-IQ-TREE 约定）
    if value == 0:
        return 0.0
    if _BARE_YEAR_RE.match(text):
        _warn_once(
            "sampling-date-bare-year",
            f"{where}: 采样日期 '{text}' 是裸四位数，无法区分公元年份与 Ma，"
            "忽略该注解；请写成 '-2019'（Ma 轴）或 '2019.0'。",
        )
        return None
    return value


def _parse_finite_age(value: Any) -> Optional[float]:
    """把年龄取值解析成有限浮点数；不是数值或为 NaN/inf 时返回 None。"""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _label_annotations(label: str) -> List[Tuple[str, str]]:
    """从 Newick 标签里抽出 ``key=value`` 注解对（键统一小写）。

    支持 ETE 的 NHX 块（``[&&NHX:age=55:tip_age=55]``，``:`` 分隔）与
    phylocom/BEAST 风格（``[&date=-0.045]``、``[&&DATE=2020-03-01]``，``;`` 分隔）。
    """
    pairs: List[Tuple[str, str]] = []
    for block in _ANNOTATION_BLOCK_RE.findall(label):
        for token in _ANNOTATION_TOKEN_RE.split(block):
            match = _ANNOTATION_PAIR_RE.match(token)
            if match:
                value = match.group(2).strip().strip('"').strip("'")
                pairs.append((match.group(1).lower(), value))
    return pairs


def get_recorded_node_age(node: Any, *, node_kind: str = "节点") -> Optional[float]:
    """读取节点上**显式记录**的年龄（距今 Ma）；没有则返回 None。

    优先级：年龄类特征（``age``/``tip_age``/…）> 日期类特征（``date``/…）>
    叶标签里的 ``[&age=…]`` / ``[&&NHX:age=…]`` 注解 > ``|日期`` 后缀。

    负年龄不会被折成正数（同审阅项 C-55 的取向：那是方向/单位写错的信号，
    折正会把不可能结果伪装成合理值），只告警并退回按枝长推算。
    """
    name = getattr(node, "name", "") or "(未命名)"
    where = f"{node_kind} '{name}'"

    for key in NODE_AGE_FEATURE_KEYS:
        raw = get_tree_node_attr(node, key, None)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            continue
        value = _parse_finite_age(raw)
        if value is None:
            _warn_once(
                "recorded-age-not-numeric",
                f"{where}: 特征 '{key}'={raw!r} 不是有限的数值年龄，已忽略。",
            )
            continue
        if value < 0:
            _warn_once(
                "recorded-negative-age",
                f"{where}: 记录的年龄 {value} Ma 为负（年龄轴上不可能），"
                "已忽略并改用分支长度推算；请检查该注解的单位或方向。",
            )
            continue
        return value

    for key in NODE_DATE_FEATURE_KEYS:
        raw = get_tree_node_attr(node, key, None)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            continue
        age = _sampling_axis_to_age_ma(raw, where=f"{where} 的特征 '{key}'")
        if age is not None:
            return age

    label = getattr(node, "name", "") or ""
    if label:
        annotations = _label_annotations(label)
        # 注解里的年龄键优先于日期键（与上面的特征顺序一致）
        for key, raw in annotations:
            if key not in NODE_AGE_FEATURE_KEYS:
                continue
            value = _parse_finite_age(raw)
            if value is None:
                _warn_once(
                    "recorded-age-not-numeric",
                    f"{where}: Newick 注解 '{key}={raw}' 不是有限的数值年龄，已忽略。",
                )
                continue
            if value < 0:
                _warn_once(
                    "recorded-negative-age",
                    f"{where}: Newick 注解 '{key}={raw}' 为负（年龄轴上不可能），"
                    "已忽略并改用分支长度推算；请检查该注解的单位或方向。",
                )
                continue
            return value
        for key, raw in annotations:
            if key in NODE_DATE_FEATURE_KEYS:
                age = _sampling_axis_to_age_ma(
                    raw, where=f"{where} 的 Newick 注解 '{key}'"
                )
                if age is not None:
                    return age
        suffix = _LABEL_DATE_SUFFIX_RE.search(label)
        if suffix:
            age = _sampling_axis_to_age_ma(
                suffix.group(1), where=f"{where} 的标签日期后缀"
            )
            if age is not None:
                return age

    return None


# --------------------------------------------------------------------------- #
# 根年龄推断（B-13）
# --------------------------------------------------------------------------- #


def _root_to_tip_heights(root_node: Any) -> List[float]:
    """所有可计算出的 root-to-tip 距离（距今 Ma / 相对单位）。"""
    heights: List[float] = []
    for leaf in (n for n in _iter_nodes(root_node) if _is_leaf(n)):
        try:
            heights.append(distance_from_root(root_node, leaf))
        except ValueError:
            # 单条路径算不出来不改变"用最长路径"这一口径，
            # 该节点稍后在坐标计算里会再次报错。
            continue
    return heights


def _max_root_to_tip_height(root_node: Any) -> Optional[float]:
    """树自身枝长隐含的根高度；无任何枝长信息时返回 None。"""
    heights = _root_to_tip_heights(root_node)
    if not heights:
        return None
    max_height = max(heights)
    return max_height if max_height > 0 else None


def infer_root_age_from_branch_lengths(root_node: Any) -> Optional[float]:
    """从分支长度推断根年龄 = 最大的 root-to-tip 路径长度（距今 Ma）。

    对超度量（时间）树这就是根年龄；对非超度量树取最长路径，同时告警并说明
    所用口径（审阅项 B-18：本层此前完全不检查输入树是否为时间树）。

    Returns
    -------
    float or None
        推断出的根年龄；树没有任何可用枝长（纯拓扑树）时返回 None。
    """
    heights = _root_to_tip_heights(root_node)
    if not heights:
        return None

    max_height = max(heights)
    if max_height <= 0:
        return None

    min_height = min(heights)
    spread = (max_height - min_height) / max_height
    if spread > ULTRAMETRIC_TOLERANCE:
        logger.warning(
            f"输入树不是时间树（非超度量）：root-to-tip 距离在 {min_height:.6g} 与 "
            f"{max_height:.6g} 之间，相对极差 {spread:.3g} > {ULTRAMETRIC_TOLERANCE:g}。"
            f"根年龄按**最长** root-to-tip 路径 = {max_height:.6g} 取值，"
            "各尖端按自身枝长落到隐含的采样时间点上（B-14）。"
            "尖端定年 / 含灭绝类群的树出现该情形是正常的，可忽略本提示；"
            "否则请先用定年方法（或 services.deep_validator.check_ultrametricity）核验。"
        )
    return max_height


def _warn_if_root_age_inconsistent_with_branch_lengths(
    root_node: Any,
    root_age: float,
    *,
    context: str,
) -> None:
    """显式 ``root_age`` 与树自身枝长隐含的根高度不一致时告警。

    节点年龄一律按 ``root_age - dist(root, node)`` 计算（B-14 之后叶节点也走
    同一公式）。若用户给的根年龄和枝长隐含高度不同，尖端就不会落在 0 Ma，
    而是落在 ``root_age - 枝长隐含高度`` —— 这是"枝长是相对单位、没有按
    root_age 归一"的典型症状，必须说出来而不是画出来。
    """
    if root_age <= 0:
        return
    tree_height = _max_root_to_tip_height(root_node)
    if tree_height is None:
        return
    if abs(tree_height - root_age) / root_age <= ULTRAMETRIC_TOLERANCE:
        return
    _warn_once(
        "root-age-vs-branch-lengths",
        f"{context}: 显式给出的 root_age={root_age:g} Ma 与树自身分支长度隐含的"
        f"根高度 {tree_height:g} 不一致。图上所有节点按 root_age - 根到该节点距离"
        "定位，因此尖端不会落在 0 Ma，而是落在各自枝长对应的时间点上。"
        "如果这些分支长度是相对单位（替代数/超时），请先按 root_age/根高度 归一化，"
        "或直接用定年后的时间树出图。",
    )


def resolve_root_age(
    root_node: Any,
    root_age: Optional[float] = None,
    *,
    context: str = "DatingFigure.add_tree",
) -> float:
    """确定用于绘图的根年龄（距今 Ma）。

    顺序：调用方显式传入 > 根节点自身的 ``age`` 注解 > **从分支长度真实推断**
    （B-13）> 树上已有节点年龄的最大值。全部不可用时显式报错——此前会静默取 0，
    导致所有节点塌到 x=0 且没有任何 warning。
    """
    if root_age is not None:
        value = _parse_finite_age(root_age)
        if value is None:
            raise ValueError(
                f"{context}: root_age={root_age!r} 不是有限的数值年龄（距今 Ma），"
                "无法绘制时间轴。"
            )
        if value < 0:
            raise ValueError(
                f"{context}: root_age={root_age!r} 为负，年龄轴上不可能；"
                "请给出正的根年龄（距今 Ma）。"
            )
        _warn_if_root_age_inconsistent_with_branch_lengths(
            root_node, value, context=context
        )
        return value

    annotated = get_tree_node_attr(root_node, "age", None)
    if annotated is not None:
        value = _parse_finite_age(annotated)
        if value is not None and value >= 0:
            logger.debug(f"{context}: 根年龄取自根节点注解 age={value}")
            return value
        logger.warning(
            f"{context}: 根节点的 age 注解 {annotated!r} 不是有效的非负有限年龄，"
            "改为按分支长度推断。"
        )

    inferred = infer_root_age_from_branch_lengths(root_node)
    if inferred is not None:
        logger.info(
            f"{context}: 未提供 root_age，按 root-to-tip 最长路径长度推断根年龄 = "
            f"{inferred:.6g} Ma。如需覆盖该值请显式传入 root_age。"
        )
        return inferred

    recorded = [
        age
        for node in _iter_nodes(root_node)
        if (age := get_recorded_node_age(node, node_kind="节点")) is not None
    ]
    if recorded:
        value = max(recorded)
        logger.warning(
            f"{context}: 树没有分支长度、也没有根年龄注解，按树上已有节点年龄的"
            f"最大值推断根年龄 = {value:.6g} Ma（共 {len(recorded)} 个节点带年龄注解）。"
        )
        return value

    raise ValueError(
        f"{context}: 无法确定根年龄——既未显式传入 root_age，树上也没有根节点年龄注解，"
        "而且所有 root-to-tip 分支长度都为 0（纯拓扑树/无枝长）。"
        "请显式传入 root_age（距今 Ma），或先做定年再出图；"
        "把根年龄当成 0 会让所有节点落在 x=0，产出一张没有信息量却看似正常的图。"
    )


def has_plot_tree_shape(tree: Any) -> bool:
    """该对象是否具备本层绘图所需的"ete 形状"。

    可视化层按鸭子类型访问节点，要求的最低契约是：
    ``traverse()``、``children``/``root``、``name``、``dist``、``up``、
    ``is_leaf``、``is_root``（``get_distance``/``ladderize`` 缺失时有回退路径）。
    """
    if not callable(getattr(tree, "traverse", None)):
        return False
    return hasattr(tree, "children") or hasattr(tree, "root")


def as_ete_tree(tree: Any) -> Any:
    """把 :class:`~phylodater.models.tree.PhylogeneticTree` 换成绘图用的树后端。

    按名字顺序探测模型层暴露的后端入口（公开名优先，向后兼容私有名），
    这样 ``models/tree.py`` 把默认后端从 ete3 换成 ete4/dendropy 时本层无需改动；
    可视化层**不**自己 import 任何树后端（审阅项 A-5/B-10）。

    Raises
    ------
    TypeError
        传入的对象既不是绘图可用的树后端，也不带任何已知的后端入口；
        或后端入口返回的对象不具备 :func:`has_plot_tree_shape` 所要求的形状
        （典型情形：返回一棵 DendroPy ``Tree``——它没有 ``traverse()``）。
    ImportError
        对象提供了后端入口，但所有入口都无法返回一棵树（典型情形：ete3 在
        Python >=3.13 上因 stdlib ``cgi`` 被移除而无法导入，见审阅项 A-5）。
    """
    if has_plot_tree_shape(tree):
        # 已经是 ete3/ete4 形状的树后端（或测试替身），直接使用
        return tree

    failures: List[str] = []
    for accessor_name in ("to_ete_tree", "ete_tree", "_ete4_tree", "_ete3_tree"):
        try:
            accessor = getattr(tree, accessor_name, None)
            backend = accessor() if callable(accessor) else accessor
        except ImportError as e:
            failures.append(f"{accessor_name}: {e}")
            continue
        if backend is None:
            continue
        if not has_plot_tree_shape(backend):
            raise TypeError(
                f"{type(tree).__name__}.{accessor_name}() 返回的 "
                f"{type(backend).__name__} 不具备绘图所需的树接口："
                "需要 traverse()/children、name、dist、up、is_leaf、is_root"
                "（ete3/ete4 的形状）。"
                "请在模型层返回 ete 形状的树，或直接把 ete 树对象传给可视化层。"
            )
        return backend

    if failures:
        raise ImportError(
            f"无法从 {type(tree).__name__} 取得可视化用的树后端："
            + "; ".join(failures)
            + "。可视化层不自行导入 ete，请安装 ete4（推荐；ete3 自 Python 3.13 起"
            "因 stdlib 'cgi' 被移除而无法导入），或让模型层提供无 ete 依赖的树后端。"
        )

    raise TypeError(
        f"无法为绘图取得 {type(tree).__name__} 的树后端："
        "该对象既不是 ete3/ete4 树，也没有 to_ete_tree()/ete_tree 入口。"
        "请直接传入 ete 树对象。"
    )


def compute_node_coordinates(
    tree: Any,
    root_age: Optional[float] = None,
) -> Dict[Any, NodeCoordinate]:
    """
    Compute x and y coordinates for all nodes in the tree.

    X coordinates represent absolute age (time before present).
    Y coordinates are assigned based on leaf ordering after ladderization.

    Parameters
    ----------
    tree : ete4.Tree or ete3.Tree
        The tree to compute coordinates for.
    root_age : float, optional
        The root age in Ma. When omitted it is resolved by
        :func:`resolve_root_age`: root-node annotation -> **longest root-to-tip
        branch path** (B-13) -> largest recorded node age. If none of these is
        available a ``ValueError`` is raised instead of silently using 0
        (which collapsed every node onto x = 0).

    Returns
    -------
    Dict[Any, NodeCoordinate]
        Dictionary mapping nodes to their coordinates.

    Raises
    ------
    ValueError
        If the root age cannot be determined, or if the distance from the root
        to some node cannot be computed (C-13: such nodes used to be drawn at
        "the present" without any warning).
    """
    ladderize = getattr(tree, "ladderize", None)
    if callable(ladderize):
        ladderize()

    # ete3的树本身就是根节点，ete4的树有root属性；特征读取统一走
    # get_tree_node_attr（ete4 在 props 里，ete3 直接是实例属性）
    root_node = tree.root if hasattr(tree, "root") else tree

    root_age = resolve_root_age(
        root_node, root_age, context="tree_plot.compute_node_coordinates"
    )

    coords: Dict[Any, NodeCoordinate] = {}
    y_counter = [0]
    age_sources: Dict[int, str] = {}

    def is_collapsed(node: Any) -> bool:
        """Check if a node should be treated as collapsed."""
        return bool(get_tree_node_attr(node, "collapsed", False))

    def assign_y(node: Any) -> float:
        """Recursively assign Y coordinates."""
        # ete3/ete4 兼容性: is_leaf和is_root在ete3中是方法，ete4中是属性
        is_leaf_node = node.is_leaf() if callable(node.is_leaf) else node.is_leaf
        is_root_node = node.is_root() if callable(node.is_root) else node.is_root

        if is_leaf_node or is_collapsed(node):
            y_counter[0] += 1
            y = float(y_counter[0])
        else:
            child_ys = [assign_y(child) for child in node.children]
            # np.mean 返回 numpy 标量，与上方分支的 float 类型不一致；
            # 显式转换以满足 mypy 的 assignment 检查。
            y = float(np.mean(child_ys))

        age, source = _node_age_with_source(node, root_node, root_age)
        x = float(age)
        age_sources[id(node)] = source

        if is_ete4_node(node):
            props = dict(node.props)
        else:
            # ete3的features是set类型，需要特殊处理
            try:
                props = {f: getattr(node, f, None) for f in node.features}
            except Exception:
                props = {}

        coords[node] = NodeCoordinate(
            x=x,
            y=y,
            name=node.name or "",
            is_leaf=is_leaf_node,
            is_root=is_root_node,
            age=age,
            props=props,
        )

        return y

    assign_y(root_node)

    _check_age_axis_consistency(coords, age_sources, root_node, root_age)

    _log_tip_age_accounting(coords, age_sources)

    return coords


def _check_age_axis_consistency(
    coords: Dict[Any, NodeCoordinate],
    age_sources: Dict[int, str],
    root_node: Any,
    root_age: float,
    tolerance: float = 1e-9,
) -> None:
    """B-14 的后置守卫：时间轴自相矛盾时必须说出来，而不是只画出来。

    两种形态都只告警（不阻断出图），因为它们的常见成因是"部分节点用注解年龄、
    其余用 ``root_age - 分支长度``"这一**混合口径**，用户往往仍希望看图：

    1. 某节点比其父节点还老 —— 即审阅项 B-14 §3 说的"枝干在时间轴上倒挂"；
    2. 根节点的注解年龄与坐标口径 ``root_age`` 不一致 —— ``add_result`` 只改写
       结果里出现的那些节点，其余节点仍按旧尺度定位。
    """
    inversions: List[str] = []
    for node, coord in coords.items():
        parent = getattr(node, "up", None)
        if parent is None:
            continue
        parent_coord = coords.get(parent)
        if parent_coord is None:
            continue
        if coord.age > parent_coord.age + tolerance:
            child_label = coord.name or "(未命名)"
            parent_label = parent_coord.name or "(未命名)"
            inversions.append(
                f"'{child_label}' {coord.age:g} Ma 比其父节点"
                f" '{parent_label}' {parent_coord.age:g} Ma 更老"
            )

    if inversions:
        preview = "; ".join(inversions[:3])
        suffix = "…" if len(inversions) > 3 else ""
        _warn_once(
            "age-axis-inversion",
            f"TreePlotter: 时间轴倒挂——{len(inversions)} 个节点比其父节点更老"
            f"（{preview}{suffix}）。这通常意味着部分节点的年龄取自注解、其余取自"
            " root_age - 分支长度，两套口径不一致；请检查注解的单位/方向，"
            "或用完整定年的树出图。",
        )

    root_coord = coords.get(root_node)
    if (
        root_coord is not None
        and age_sources.get(id(root_node)) == "recorded"
        and abs(root_coord.age - root_age) > tolerance
    ):
        _warn_once(
            "root-age-vs-coordinate-scale",
            f"TreePlotter: 根节点记录的年龄 {root_coord.age:g} Ma 与坐标口径"
            f" root_age={root_age:g} Ma 不同。没有年龄注解的节点仍按"
            " root_age - 分支长度 定位，因此图上会同时出现两套时间尺度。"
            "请以定年结果为准重新出图，或显式传入与注解一致的 root_age。",
        )


def _log_tip_age_accounting(
    coords: Dict[Any, NodeCoordinate],
    age_sources: Dict[int, str],
) -> None:
    """报告有多少尖端用了树内记录的真实年龄（B-14 要求"报告有多少叶被赋值"）。"""
    tips = [
        (node, coord)
        for node, coord in coords.items()
        if coord.is_leaf and not coord.is_root
    ]
    if not tips:
        return

    recorded = [
        coord.name for node, coord in tips if age_sources.get(id(node)) == "recorded"
    ]
    unknown = [
        coord.name for node, coord in tips if age_sources.get(id(node)) == "unknown"
    ]

    if recorded:
        preview = ", ".join(name or "(未命名)" for name in recorded[:3])
        logger.info(
            f"TreePlotter: {len(recorded)}/{len(tips)} 个尖端按树内记录的采样日期/"
            f"化石年龄绘制（如 {preview}）。"
        )
    if unknown:
        preview = ", ".join(name or "(未命名)" for name in unknown[:5])
        suffix = "..." if len(unknown) > 5 else ""
        logger.warning(
            f"TreePlotter: {len(unknown)} 个节点（共 {len(tips)} 个尖端）既没有年龄"
            "注解、其到根路径上也没有任何分支长度信息，只能按'现在'(0 Ma) 绘制："
            f"{preview}{suffix}。若其中含灭绝类群或带采样日期的序列，它们会被错误地"
            "画到今天——请在树里提供 tip age（age 特征或 '|-0.045' 形式的采样日期）。"
        )


def _node_age_with_source(
    node: Any,
    root: Any,
    root_age: float,
) -> Tuple[float, str]:
    """节点年龄（距今 Ma）及其来源标记。

    Returns
    -------
    (age, source)
        ``source`` 为 ``"recorded"``（树里显式记录的年龄）/ ``"root"`` /
        ``"branch-length"`` / ``"unknown"``（无任何年龄信息，只能按 0 绘制）。
    """
    # 显式记录的年龄优先——含 add_result 之后覆盖根节点年龄的情形。
    is_leaf = _is_leaf(node)
    recorded = get_recorded_node_age(node, node_kind="尖端" if is_leaf else "节点")
    if recorded is not None:
        return max(0.0, recorded), "recorded"

    if node is root:
        return float(root_age), "root"

    if not _path_has_branch_length_info(root, node):
        return 0.0, "unknown"

    # 叶节点与内部节点用同一套公式（B-14）：超度量时间树上尖端自然回到 0，
    # 含灭绝类群 / 非超度量输入下则落在各自枝长隐含的时间点上。
    distance = distance_from_root(root, node)
    age = root_age - distance
    if not math.isfinite(age):
        raise ValueError(
            f"节点 '{getattr(node, 'name', '') or '(未命名)'}' 的年龄 "
            f"root_age({root_age}) - dist({distance}) = {age} 不是有限值。"
        )
    return max(0.0, age), "branch-length"


def _path_has_branch_length_info(root: Any, node: Any) -> bool:
    """该节点到根的路径上是否存在任何枝长记录。

    ete3 解析无枝长的 Newick 时把 ``dist`` 置为 0.0（视为"有信息"，由
    :func:`resolve_root_age` 的"全 0 即无法推断"守卫兜住）；只有当路径上
    所有 ``dist`` 都缺失（None / 非数值）时才判为"没有任何年龄信息"。
    """
    current = node
    seen = {id(current)}
    while current is not root:
        parent = getattr(current, "up", None)
        if parent is None or id(parent) in seen:
            # 父链不完整/成环：交给 distance_from_root 去报明确错误
            return True
        if _branch_length_to_parent(current) is not None:
            return True
        seen.add(id(parent))
        current = parent
    return False


def _compute_node_age(node: Any, root: Any, root_age: float) -> float:
    """
    Compute the age of a node in Ma (time before present).

    叶节点与内部节点共用 ``root_age - dist(root, node)``；节点自身带有采样
    日期 / 化石年龄注解时优先用注解（B-14）。

    Parameters
    ----------
    node : ete4.TreeNode or ete3.TreeNode
        The node to compute age for.
    root : ete4.TreeNode or ete3.TreeNode
        The root of the tree.
    root_age : float
        The age of the root node.

    Returns
    -------
    float
        The computed age.

    Raises
    ------
    ValueError
        If the root-to-node distance cannot be computed. Silently returning
        0.0 used to draw "failed" identically to "happened today" (C-13).
    """
    return _node_age_with_source(node, root, root_age)[0]


def get_tree_segments(
    tree: Any,
    coords: Optional[Dict[Any, NodeCoordinate]] = None,
    root_age: Optional[float] = None,
) -> Tuple[List[List[Tuple[float, float]]], List[Dict[str, Any]]]:
    """
    Generate line segments for drawing the tree.

    Creates L-shaped connections (rectangular tree style) for all branches.

    Parameters
    ----------
    tree : ete4.Tree
        The tree to generate segments for.
    coords : Dict[Any, NodeCoordinate], optional
        Pre-computed node coordinates. If not provided, will be computed.
    root_age : float, optional
        The root age for coordinate computation.

    Returns
    -------
    Tuple[List[List[Tuple[float, float]]], List[Dict[str, Any]]]
        A tuple of (segments, segment_properties) where each segment is a
        list of (x, y) coordinate tuples, and properties contain metadata.
    """
    if coords is None:
        coords = compute_node_coordinates(tree, root_age)

    segments = []
    properties = []

    for node in tree.traverse():
        # ete3/ete4 兼容性
        is_root_node = node.is_root() if callable(node.is_root) else node.is_root
        is_leaf_node = node.is_leaf() if callable(node.is_leaf) else node.is_leaf

        if is_root_node:
            continue

        node_coord = coords.get(node)
        parent_coord = coords.get(node.up)

        if node_coord is None or parent_coord is None:
            continue

        x_child = node_coord.x
        y_child = node_coord.y
        x_parent = parent_coord.x

        horizontal_segment = [(x_parent, y_child), (x_child, y_child)]
        segments.append(horizontal_segment)

        props = dict(node_coord.props) if node_coord.props else {}
        props["is_horizontal"] = True
        props["node_name"] = node_coord.name
        props["is_leaf"] = is_leaf_node
        properties.append(props)

    for node in tree.traverse():
        is_leaf_node = node.is_leaf() if callable(node.is_leaf) else node.is_leaf
        if is_leaf_node:
            continue

        node_coord = coords.get(node)
        if node_coord is None:
            continue

        child_ys = []
        for child in node.children:
            child_coord = coords.get(child)
            if child_coord is not None:
                child_ys.append(child_coord.y)

        if len(child_ys) >= 1:
            x_parent = node_coord.x
            # 垂直线段：从父节点y坐标到子节点y坐标的范围
            min_y = min(min(child_ys), node_coord.y)
            max_y = max(max(child_ys), node_coord.y)
            vertical_segment = [(x_parent, min_y), (x_parent, max_y)]
            segments.append(vertical_segment)

            props = dict(node_coord.props) if node_coord.props else {}
            props["is_horizontal"] = False
            props["is_vertical"] = True
            props["node_name"] = node_coord.name
            properties.append(props)

    return segments, properties


class TreePlotter:
    """
    Class for plotting phylogenetic trees on matplotlib axes.

    This class handles the visualization of tree topology, branch styling,
    node markers, and labels.
    """

    def __init__(
        self,
        branch_color: str = "#2C3E50",
        branch_linewidth: float = 2.0,
        leaf_color: str = "#E74C3C",
        leaf_size: float = 8.0,
        internal_node_size: float = 4.0,
        show_internal_labels: bool = False,
        show_leaf_labels: bool = True,
        label_offset: float = 1.5,
        label_fontsize: int = 10,
    ) -> None:
        """
        Initialize the TreePlotter.

        Parameters
        ----------
        branch_color : str
            Default color for branches.
        branch_linewidth : float
            Default linewidth for branches.
        leaf_color : str
            Default color for leaf nodes.
        leaf_size : float
            Size of leaf node markers.
        internal_node_size : float
            Size of internal node markers.
        show_internal_labels : bool
            Whether to show labels for internal nodes.
        show_leaf_labels : bool
            Whether to show labels for leaf nodes.
        label_offset : float
            Horizontal offset for labels.
        label_fontsize : int
            Font size for labels.
        """
        self.branch_color = branch_color
        self.branch_linewidth = branch_linewidth
        self.leaf_color = leaf_color
        self.leaf_size = leaf_size
        self.internal_node_size = internal_node_size
        self.show_internal_labels = show_internal_labels
        self.show_leaf_labels = show_leaf_labels
        self.label_offset = label_offset
        self.label_fontsize = label_fontsize

        self._color_mapping: Dict[str, Any] = {}
        # ``map_size`` 往这里同时写“列名（str）”与“区间（Tuple[float, float]）”，
        # 旧标注 Dict[str, float] 与两个写入都不符。
        self._size_mapping: Dict[str, Any] = {}

    def plot(
        self,
        tree: Any,
        ax: Optional[Axes] = None,
        output: Optional[Union[str, Path]] = None,
        coords: Optional[Dict[Any, NodeCoordinate]] = None,
        root_age: Optional[float] = None,
        color_column: Optional[str] = None,
        size_column: Optional[str] = None,
        cmap: str = "viridis",
        palette: Optional[List[str]] = None,
        size_range: Tuple[float, float] = (2.0, 8.0),
        figsize: Tuple[float, float] = (10.0, 8.0),
    ) -> Dict[Any, NodeCoordinate]:
        """
        Plot the tree on the given axes or save to a file.

        Parameters
        ----------
        tree : ete4.Tree
            The tree to plot.
        ax : matplotlib.axes.Axes, optional
            The axes to plot on. If None and ``output`` is provided, a new
            figure and axes will be created internally.
        output : str or pathlib.Path, optional
            File path to save the plot. When given, a new figure is created
            and saved via ``plt.savefig``.
        coords : Dict[Any, NodeCoordinate], optional
            Pre-computed coordinates.
        root_age : float, optional
            The root age.
        color_column : str, optional
            Column name for color mapping.
        size_column : str, optional
            Column name for size mapping.
        cmap : str
            Matplotlib colormap name.
        palette : List[str], optional
            List of colors for categorical variables.
        size_range : Tuple[float, float]
            Range for size scaling.
        figsize : Tuple[float, float]
            Figure size used when creating a new figure.

        Returns
        -------
        Dict[Any, NodeCoordinate]
            The computed coordinates.

        Raises
        ------
        ValueError
            ``coords`` 未提供时，根年龄无法确定或某节点到根的距离算不出来
            （见 :func:`compute_node_coordinates`）。宁可不画图，也不画出一张
            把所有节点压在 x=0 的假图（B-13 / C-13）。
        """
        # 自动适配 PhylogeneticTree 值对象：取其树后端（不硬编码 ete3）
        tree = as_ete_tree(tree)

        created_figure = False
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
            created_figure = True

        # 应用通过 map_color / map_size 设置的链式映射
        if color_column is None and self._color_mapping.get("column"):
            color_column = self._color_mapping["column"]
            cmap = self._color_mapping.get("cmap", cmap)
            palette = self._color_mapping.get("palette", palette)
        if size_column is None and self._size_mapping.get("column"):
            size_column = self._size_mapping["column"]
            size_range = self._size_mapping.get("range", size_range)

        if coords is None:
            coords = compute_node_coordinates(tree, root_age)

        segments, properties = get_tree_segments(tree, coords, root_age)

        colors = self._compute_colors(segments, properties, color_column, cmap, palette)

        linewidths = self._compute_linewidths(
            segments, properties, size_column, size_range
        )

        if segments:
            lc = LineCollection(
                segments,
                colors=colors,
                linewidths=linewidths,
                zorder=1,
            )
            ax.add_collection(lc)

        self._plot_nodes(ax, coords, tree)

        self._plot_labels(ax, coords, tree)

        self._set_axes_style(ax, coords)

        if output is not None:
            output_path = Path(output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_path, bbox_inches="tight")
            if created_figure:
                plt.close(fig)
        elif created_figure:
            # 未指定 output 且自动创建了 figure：关闭以避免句柄/内存泄漏。
            # 调用方如需保留 figure，应自行传入 ax。
            plt.close(fig)

        return coords

    def _compute_colors(
        self,
        segments: List[List[Tuple[float, float]]],
        properties: List[Dict[str, Any]],
        color_column: Optional[str],
        cmap: str,
        palette: Optional[List[str]],
    ) -> List[Union[str, Tuple[float, float, float, float]]]:
        """Compute colors for each segment."""
        import matplotlib.colors as mcolors

        default_color = self.branch_color

        if color_column is None:
            return [default_color] * len(segments)

        values: List[Any] = []
        for props in properties:
            value = props.get(color_column)
            values.append(value)

        unique_values: Set[Any] = set(v for v in values if v is not None)

        if all(isinstance(v, (int, float)) for v in unique_values if v is not None):
            numeric_values = [v if v is not None else 0 for v in values]
            min_val = min(numeric_values)
            max_val = max(numeric_values)

            if max_val == min_val:
                normalized = [0.5] * len(numeric_values)
            else:
                normalized = [
                    (v - min_val) / (max_val - min_val) for v in numeric_values
                ]

            colormap = plt.colormaps[cmap]
            return [colormap(n) for n in normalized]
        else:
            # matplotlib 的内置调色板是 RGBA 元组而不是颜色名，所以另用一个
            # 局部变量承接，不污染入参 ``palette`` 的 List[str] 契约。
            palette_or_default: Sequence[Any]
            if palette is None:
                palette_or_default = list(mcolors.TABLEAU_COLORS.values())
            else:
                palette_or_default = palette

            value_to_color: Dict[Any, Any] = {}
            for i, val in enumerate(sorted(unique_values)):
                value_to_color[val] = palette_or_default[i % len(palette_or_default)]

            return [value_to_color.get(v, default_color) for v in values]

    def _compute_linewidths(
        self,
        segments: List[List[Tuple[float, float]]],
        properties: List[Dict[str, Any]],
        size_column: Optional[str],
        size_range: Tuple[float, float],
    ) -> List[float]:
        """Compute linewidths for each segment."""
        if size_column is None:
            return [self.branch_linewidth] * len(segments)

        values = []
        for props in properties:
            value = props.get(size_column)
            if value is not None:
                try:
                    values.append(float(value))
                except (TypeError, ValueError):
                    values.append(0)
            else:
                values.append(0)

        min_val = min(values) if values else 0
        max_val = max(values) if values else 1

        if max_val == min_val:
            normalized = [0.5] * len(values)
        else:
            normalized = [(v - min_val) / (max_val - min_val) for v in values]

        min_size, max_size = size_range
        return [min_size + n * (max_size - min_size) for n in normalized]

    def _plot_nodes(
        self,
        ax: Axes,
        coords: Dict[Any, NodeCoordinate],
        tree: Any,
    ) -> None:
        """Plot node markers."""
        leaf_x = []
        leaf_y = []
        internal_x = []
        internal_y = []
        root_x = []
        root_y = []

        for _node, coord in coords.items():
            if coord.is_root:
                root_x.append(coord.x)
                root_y.append(coord.y)
            elif coord.is_leaf:
                leaf_x.append(coord.x)
                leaf_y.append(coord.y)
            else:
                internal_x.append(coord.x)
                internal_y.append(coord.y)

        if leaf_x and self.leaf_size > 0:
            ax.scatter(
                leaf_x,
                leaf_y,
                s=self.leaf_size**2,
                c=self.leaf_color,
                zorder=2,
            )

        if internal_x and self.internal_node_size > 0:
            ax.scatter(
                internal_x,
                internal_y,
                s=self.internal_node_size**2,
                c=self.branch_color,
                zorder=2,
            )

        if root_x:
            ax.scatter(
                root_x,
                root_y,
                s=(self.leaf_size * 1.5) ** 2,
                c=self.branch_color,
                zorder=2,
            )

    def _plot_labels(
        self,
        ax: Axes,
        coords: Dict[Any, NodeCoordinate],
        tree: Any,
    ) -> None:
        """Plot node labels."""
        for _node, coord in coords.items():
            if coord.is_leaf and self.show_leaf_labels:
                ax.text(
                    coord.x - self.label_offset,
                    coord.y,
                    coord.name,
                    fontsize=self.label_fontsize,
                    va="center",
                    ha="left",
                )
            elif not coord.is_leaf and self.show_internal_labels and coord.name:
                ax.text(
                    coord.x,
                    coord.y,
                    coord.name,
                    fontsize=self.label_fontsize - 1,
                    va="center",
                    ha="center",
                )

    def _set_axes_style(
        self,
        ax: Axes,
        coords: Dict[Any, NodeCoordinate],
    ) -> None:
        """Set axes style and limits."""
        if not coords:
            return

        x_values = [c.x for c in coords.values()]
        y_values = [c.y for c in coords.values()]

        x_min, x_max = min(x_values), max(x_values)
        y_min, y_max = min(y_values), max(y_values)

        x_padding = (x_max - x_min) * 0.05 + 1
        y_padding = (y_max - y_min) * 0.05

        ax.set_xlim(x_max + x_padding, max(0, x_min - x_padding))
        ax.set_ylim(y_min - y_padding, y_max + y_padding)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

        ax.tick_params(
            left=False, labelleft=False, bottom=False, top=False, labelbottom=False
        )
        ax.set_xticks([])

    def map_color(
        self,
        column_name: str,
        cmap: str = "viridis",
        palette: Optional[List[str]] = None,
        target: str = "branch",
    ) -> "TreePlotter":
        """
        Set up color mapping for a column.

        Parameters
        ----------
        column_name : str
            The column to map colors to.
        cmap : str
            Matplotlib colormap name.
        palette : List[str], optional
            List of colors for categorical variables.
        target : str
            Target to color: "branch", "label", or "both".

        Returns
        -------
        TreePlotter
            Self for method chaining.
        """
        self._color_mapping = {
            "column": column_name,
            "cmap": cmap,
            "palette": palette,
            "target": target,
        }
        return self

    def map_size(
        self,
        column_name: str,
        size_range: Tuple[float, float] = (2.0, 8.0),
    ) -> "TreePlotter":
        """
        Set up size mapping for a column.

        Parameters
        ----------
        column_name : str
            The column to map sizes to.
        size_range : Tuple[float, float]
            Range for size scaling.

        Returns
        -------
        TreePlotter
            Self for method chaining.
        """
        self._size_mapping = {
            "column": column_name,
            "range": size_range,
        }
        return self
