"""
约束 -> 目标软件格式渲染的单元测试（wave-1 #3 修复验证）

验证 AgeConstraint.to_software_format 在传入节点真实字段后：
1. 渲染结果不含任何占位符（{label}/{Name}/{t1}/{t2}/{taxa_str}）；
2. 包含真实 taxa/年龄字段。

这是将旧实现「返回未解析占位符模板、由适配器各自 .format(...)」改为
「模型内部一次性渲染真实字段」的核心正确性验证。
"""

import math
import warnings

import pytest

from phylodater.core.exceptions import SemanticDegradationWarning
from phylodater.models.constraints import (
    FixedAgeConstraint,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    SkewNormalConstraint,
    SkewTConstraint,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)

# 需要节点字段的 6 个目标软件（mcmctree 单独测试，无占位符）
TARGETS = ["r8s", "treepl", "pathd8", "lsd2", "wlogdate", "mdcat"]

# 与模板占位符对应的真实节点字段
CTX = dict(label="TaxA", name="Cal1", t1="Leaf1", t2="Leaf2", taxa_str="Leaf1+Leaf2")

CONSTRAINTS = [
    FixedAgeConstraint(fixed_age=50.0),
    UniformAgeConstraint(min_age=10.0, max_age=50.0),
    SoftLowerBoundConstraint(min_age=20.0),
    MaximumAgeConstraint(max_age=100.0),
    SoftBoundsConstraint(min_age=10.0, max_age=100.0),
    GammaPriorConstraint(alpha=2.0, beta=0.5),
    SkewNormalConstraint(location=50.0, scale=10.0, shape=2.0),
    SkewTConstraint(location=50.0, scale=10.0, shape=2.0, df=5.0),
]


def _render(constraint, target):
    """渲染并抑制降级警告（不影响占位符/字段断言）。"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SemanticDegradationWarning)
        return constraint.to_software_format(target, **CTX)


@pytest.mark.parametrize("constraint", CONSTRAINTS, ids=lambda c: type(c).__name__)
@pytest.mark.parametrize("target", TARGETS)
def test_rendered_output_has_no_placeholder(constraint, target):
    """渲染结果不得残留任何 { 占位符。"""
    out = _render(constraint, target)
    assert "{" not in out, f"{type(constraint).__name__}/{target} 仍含占位符: {out!r}"


@pytest.mark.parametrize("constraint", CONSTRAINTS, ids=lambda c: type(c).__name__)
@pytest.mark.parametrize("target", [t for t in TARGETS if t != "lsd2"])
def test_rendered_output_contains_real_node_fields(constraint, target):
    """渲染结果应包含真实节点字段（label/name/t1 之一）。

    lsd2 的日期文件格式仅含年龄边界、不含节点名，单独处理。
    """
    out = _render(constraint, target)
    assert ("TaxA" in out) or ("Cal1" in out) or ("Leaf1" in out)


def test_fixed_age_r8s_full_render():
    out = _render(FixedAgeConstraint(fixed_age=50.0), "r8s")
    assert out == "fixage taxon=TaxA age=50.0;"


def test_uniform_age_pathd8_full_render():
    out = _render(UniformAgeConstraint(min_age=10.0, max_age=50.0), "pathd8")
    assert "mrca: Leaf1, Leaf2, minage=10.0;" in out
    assert "mrca: Leaf1, Leaf2, maxage=50.0;" in out


def test_wlogdate_mdcat_taxa_and_age():
    out = _render(UniformAgeConstraint(min_age=10.0, max_age=50.0), "wlogdate")
    # 降级为区间中点 30.0，格式 Name=taxa_str\\tage
    assert out == "Cal1=Leaf1+Leaf2\t30.0"


def test_mcmctree_format_has_no_placeholder():
    """mcmctree 走 to_mcmctree_calib_string，无节点占位符。"""
    out = FixedAgeConstraint(fixed_age=1000.0).to_software_format("mcmctree")
    assert "{" not in out
    assert out.startswith("B(")


# --------------------------------------------------------------------------- #
# MCMCTree 先验串的"槽位数"快照（审阅项 A-3/A-4/B-1/B-2 的回归锁）
#
# 上游 treesub.c 的 npfossils[] 表（mcmctree.c:202-204）：
#     fossils[]   = {" ", "L", "U", "B", "G", "SN", "ST", "S2N"}
#     npfossils[] = { 0,    4,   2,   4,   2,   3,    4,     7 }
# 串里参数个数与槽位不符时，多写的值被静默忽略、少写的槽停在上游默认值。
# --------------------------------------------------------------------------- #
MCMCTREE_SLOT_COUNTS = {
    "B": 4,  # tL, tU, tailL, tailR
    "L": 4,  # tL, P, c, tailL
    "U": 2,  # tU, tailR
    "G": 2,  # alpha, beta（无 offset 槽）
    "SN": 3,  # loc, scale, shape
    "ST": 4,  # loc, scale, shape, df
}

MCMCTREE_SAMPLES = [
    ("B", FixedAgeConstraint(fixed_age=100.0)),
    ("B", UniformAgeConstraint(min_age=50.0, max_age=100.0)),
    ("B", SoftBoundsConstraint(min_age=50.0, max_age=100.0)),
    ("L", SoftLowerBoundConstraint(min_age=50.0)),
    ("U", MaximumAgeConstraint(max_age=100.0)),
    ("G", GammaPriorConstraint(alpha=2.0, beta=0.5)),
    ("SN", SkewNormalConstraint(location=50.0, scale=10.0, shape=2.0)),
    ("ST", SkewTConstraint(location=50.0, scale=10.0, shape=2.0, df=5.0)),
]


@pytest.mark.parametrize("keyword,constraint", MCMCTREE_SAMPLES, ids=lambda v: str(v))
def test_mcmctree_string_matches_upstream_slot_count(keyword, constraint):
    out = constraint.to_mcmctree_calib_string()
    assert out.startswith(keyword + "("), out
    params = out[len(keyword) + 1 : -1].split(",")
    assert len(params) == MCMCTREE_SLOT_COUNTS[keyword], f"{keyword} 槽位数不符: {out}"
    for token in params:
        value = float(token)
        assert math.isfinite(value), f"{keyword} 串含非有限值: {out}"
    # 区间型串必须严格递增（A-3：0.1/0.5/1.0 Ma 级别不得坍缩或颠倒）
    if keyword == "B":
        assert float(params[1]) > float(params[0]) >= 0.0, out


@pytest.mark.parametrize("keyword,constraint", MCMCTREE_SAMPLES, ids=lambda v: str(v))
def test_mcmctree_string_free_of_bound_hijack_characters(keyword, constraint):
    """注解串不得含 '<' / '>'：上游在匹配 L/U/B/G 之前先扫这两个字符（treesub.c:8690-8717）"""
    out = constraint.to_mcmctree_calib_string()
    assert "<" not in out and ">" not in out, out


def test_tiny_calibrations_do_not_collapse_when_rendered():
    """A-3 的实际危害面：第四纪/tip-dating 级别的校准（<5 Ma）渲染后仍可分辨"""
    for age in (0.05, 0.1, 0.5, 1.0, 2.0):
        out = FixedAgeConstraint(fixed_age=age).to_mcmctree_calib_string()
        parts = out[2:-1].split(",")
        assert float(parts[1]) > float(parts[0]), f"{age} Ma 坍缩: {out}"
    for lo, hi in ((0.01, 0.05), (0.1, 0.5), (1.0, 2.0)):
        out = UniformAgeConstraint(min_age=lo, max_age=hi).to_mcmctree_calib_string()
        parts = out[2:-1].split(",")
        assert float(parts[1]) > float(parts[0]), f"[{lo},{hi}] 坍缩: {out}"
