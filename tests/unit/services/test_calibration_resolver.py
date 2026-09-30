"""
CalibrationResolver 服务单元测试

测试校准点解析服务的各项功能。

重点覆盖全项目最核心的科学不变量——"祖先节点必须不早于后代节点"：

* 审阅项 B-7：冲突判定**不得依赖校准点在 YAML 里的书写顺序**（每一对都要查两个方向）。
* 审阅项 B-8：祖先关系必须由 **MRCA 的完整后代叶集合** 决定，而不是用户为求 MRCA
  指定的代表类群（``resolved_taxa``）。"鸟/恐龙混选"与"两种剑龙类"这类代表类群
  互不重叠的嵌套写法此前被系统性放过。
* 审阅项 B-9：本文件旧的 ``test_validate_ancestry`` 因为传参个数错误（被测函数只收
  一个位置参数）外加 ``except Exception: pass``，从来没有执行过被测代码的一行。
  下面的测试全部用**真实对象**（不 ``Mock(spec=...)``），并且真的调用
  ``validate_ancestry``。
"""

import itertools
from unittest.mock import Mock

import pytest

from phylodater.core.exceptions import CalibrationConflictError
from phylodater.models import (
    CalibrationPoint,
    FixedAgeConstraint,
    MaximumAgeConstraint,
    PhylogeneticTree,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)
from phylodater.models.constraints import SkewNormalConstraint, SkewTConstraint
from phylodater.services.calibration_resolver import CalibrationResolver

# --------------------------------------------------------------------------- #
# 测试拓扑
#
#   (((Triceratops, Tyrannosaurus), (Stegosaurus, Ankylosaurus)),
#    (Crocodylus, Archaeopteryx));
#
# 用**代表类群**指定两个节点：
#   祖先节点 Archosauria ← 代表类群 {Triceratops, Archaeopteryx}，其 MRCA 是全树根
#   后代节点 Stegosauria ← 代表类群 {Stegosaurus, Ankylosaurus}
# 两组代表类群**完全不相交**，但 Stegosauria 确实嵌在 Archosauria 之内——
# 这正是审阅项 B-8 指出的系统性假阴性场景（古生物学里最常见的写法）。
# --------------------------------------------------------------------------- #
TRICERATOPS = "Triceratops"
TYRANNOSAURUS = "Tyrannosaurus"
STEGOSAURUS = "Stegosaurus"
ANKYLOSAURUS = "Ankylosaurus"
CROCODYLUS = "Crocodylus"
ARCHAEOPTERYX = "Archaeopteryx"

NESTING_NEWICK = (
    "((({t1}:1.0,{t2}:1.0)Dinosauria:1.0,({t3}:1.0,{t4}:1.0)Stegosauria:1.0)"
    "Ornithischia:1.0,({t5}:1.0,{t6}:1.0)Archosauria_crown:1.0)Root:0.0;"
).format(
    t1=TRICERATOPS,
    t2=TYRANNOSAURUS,
    t3=STEGOSAURUS,
    t4=ANKYLOSAURUS,
    t5=CROCODYLUS,
    t6=ARCHAEOPTERYX,
)
ALL_TIPS = [
    TRICERATOPS,
    TYRANNOSAURUS,
    STEGOSAURUS,
    ANKYLOSAURUS,
    CROCODYLUS,
    ARCHAEOPTERYX,
]
# 由内到外的真实 clade（叶集合）
CLADES = [
    {TRICERATOPS, TYRANNOSAURUS},
    {STEGOSAURUS, ANKYLOSAURUS},
    {TRICERATOPS, TYRANNOSAURUS, STEGOSAURUS, ANKYLOSAURUS},
    {CROCODYLUS, ARCHAEOPTERYX},
    set(ALL_TIPS),
]


class _TopologyStub:
    """只提供 ``PhylogeneticTree`` 公开接口的一张固定拓扑表（不依赖任何树解析后端）。

    ``get_mrca`` 的契约与真实实现一致：返回该 MRCA 下的**全部**叶节点名，
    无法定位（含幽灵分类单元）时返回 ``None``。这样 B-7/B-8 的守卫逻辑可以在任何
    解释器上被确定性地测试；:class:`TestAncestryGuardRealTree` 再用真实
    ``PhylogeneticTree`` 把同一批场景跑一遍。
    """

    def __init__(self, newick, tips, clades):
        self.newick = newick
        self.tip_names = list(tips)
        self._tips = set(tips)
        self._clades = [set(c) for c in clades]
        self.mrca_calls = []

    def get_mrca(self, tip_names):
        taxa = set(tip_names)
        self.mrca_calls.append(sorted(taxa))
        if not taxa or not taxa.issubset(self._tips):
            return None  # 幽灵分类单元 → 无法定位（真实实现同样返回 None）
        for clade in sorted(self._clades, key=len):  # 最小的包含者 = MRCA
            if taxa.issubset(clade):
                return sorted(clade)
        return None

    def get_mrca_terminals(self, tip_names):
        terminals = self.get_mrca(tip_names)
        if not terminals:
            return None
        return (terminals[0], terminals[-1])

    def get_mrca_terminals_all(self, tip_names):
        return self.get_mrca(tip_names)

    def get_distance_matrix(self):
        return {}


@pytest.fixture
def stub_tree():
    return _TopologyStub(NESTING_NEWICK, ALL_TIPS, CLADES)


@pytest.fixture
def stub_resolver(stub_tree):
    return CalibrationResolver(stub_tree)


def _point(name, representatives, constraint, is_root=False):
    """构造真实 CalibrationPoint（不用 Mock，确保被测代码真的跑起来）。"""
    return CalibrationPoint(
        name=name,
        age_constraint=constraint,
        resolved_taxa=list(representatives),
        is_root_node=is_root,
    )


def _violating_pair():
    """祖先上界 150 Ma、后代下界 200 Ma —— 后代比祖先还老，物理上不可能。"""
    ancestor = _point(
        "Archosauria",
        [TRICERATOPS, ARCHAEOPTERYX],
        UniformAgeConstraint(100.0, 150.0),
    )
    descendant = _point(
        "Stegosauria",
        [STEGOSAURUS, ANKYLOSAURUS],
        UniformAgeConstraint(180.0, 200.0),
    )
    return ancestor, descendant


def _consistent_pair():
    ancestor = _point(
        "Archosauria",
        [TRICERATOPS, ARCHAEOPTERYX],
        UniformAgeConstraint(200.0, 250.0),
    )
    descendant = _point(
        "Stegosauria",
        [STEGOSAURUS, ANKYLOSAURUS],
        UniformAgeConstraint(150.0, 180.0),
    )
    return ancestor, descendant


class TestCalibrationResolver:
    """CalibrationResolver 测试类"""

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        return tree

    @pytest.fixture
    def resolver(self, mock_tree):
        """创建解析器实例"""
        return CalibrationResolver(mock_tree)

    def test_resolver_init(self, mock_tree):
        """测试解析器初始化"""
        resolver = CalibrationResolver(mock_tree)
        assert resolver.tree == mock_tree
        assert resolver.tip_names == ["human", "chimp", "mouse", "rat"]

    def test_resolver_init_none_tree(self):
        """测试None树初始化"""
        with pytest.raises(TypeError):
            CalibrationResolver(None)

    def test_resolve_mrca_pair(self, resolver):
        """测试解析 MRCA 对 - 使用存在的分类群"""
        # 使用单个存在的分类群进行测试
        # resolve方法需要找到至少2个匹配的叶节点
        with pytest.raises(Exception):
            resolver.resolve("human")

    def test_resolve_root_calibration(self, resolver):
        """测试解析根节点校准"""
        # "root"不是叶节点名称，应该失败
        with pytest.raises(Exception):
            resolver.resolve("root")

    def test_resolve_with_multiple_taxa(self, resolver):
        """测试解析包含多个分类群的目标"""
        # "primates"不是叶节点名称，应该失败
        with pytest.raises(Exception):
            resolver.resolve("primates")

    def test_resolve_missing_taxa(self, resolver):
        """测试解析缺失的分类群"""
        with pytest.raises(Exception):
            resolver.resolve("missing_taxon_xyz_12345")

    def test_resolve_empty_target(self, resolver):
        """测试解析空目标"""
        with pytest.raises(ValueError):
            resolver.resolve("")

    def test_resolve_none_target(self, resolver):
        """测试解析None目标"""
        with pytest.raises(TypeError):
            resolver.resolve(None)


class TestAncestryGuardOrdering:
    """B-7：冲突判定不得依赖校准点的书写顺序。"""

    def test_ancestor_older_than_descendant_passes(self, stub_resolver):
        """合法嵌套（祖先更老）→ 通过并返回 True。"""
        ancestor, descendant = _consistent_pair()

        assert stub_resolver.validate_ancestry([ancestor, descendant]) is True

    def test_violating_nesting_raises(self, stub_resolver):
        """违背不变量（后代 200 Ma > 祖先上界 150 Ma）→ 必须 CalibrationConflictError。"""
        ancestor, descendant = _violating_pair()

        with pytest.raises(CalibrationConflictError) as excinfo:
            stub_resolver.validate_ancestry([ancestor, descendant])

        message = str(excinfo.value)
        assert "Archosauria" in message and "Stegosauria" in message
        # 报错必须给出双方约束与具体边界，便于用户定位 YAML 的哪一行
        assert "150" in message and "180" in message

    def test_violating_nesting_raises_when_descendant_listed_first(self, stub_resolver):
        """B-7 核心：把后代写在祖先**之前**同样必须报错（旧实现此时完全放过）。"""
        ancestor, descendant = _violating_pair()

        with pytest.raises(CalibrationConflictError):
            stub_resolver.validate_ancestry([descendant, ancestor])

    @pytest.mark.parametrize("order", list(itertools.permutations(range(3))))
    def test_verdict_is_invariant_under_every_permutation(self, stub_resolver, order):
        """三种校准点的**全部 6 种书写顺序**都得出同一判定（报错）。

        这锁死的是"顺序敏感"这一整类缺陷，而不是某一种排列。
        """
        violator = _point(
            "Archosauria",
            [TRICERATOPS, ARCHAEOPTERYX],
            UniformAgeConstraint(100.0, 150.0),
        )
        violated = _point(
            "Stegosauria",
            [STEGOSAURUS, ANKYLOSAURUS],
            UniformAgeConstraint(180.0, 200.0),
        )
        unrelated = _point(
            "Dinosauria", [TRICERATOPS, TYRANNOSAURUS], MaximumAgeConstraint(4000.0)
        )
        points = [violator, violated, unrelated]

        with pytest.raises(CalibrationConflictError):
            stub_resolver.validate_ancestry([points[i] for i in order])

    @pytest.mark.parametrize("order", list(itertools.permutations(range(2))))
    def test_consistent_input_passes_in_every_order(self, stub_resolver, order):
        """反过来：合法的校准集合在任何顺序下都不该误报。"""
        ancestor, descendant = _consistent_pair()
        points = [ancestor, descendant]

        assert stub_resolver.validate_ancestry([points[i] for i in order]) is True

    def test_pair_comparison_asks_the_topology_backend(self, stub_resolver):
        """配对枚举确实调用了拓扑判定（而不是靠 resolved_taxa 猜）。"""
        ancestor, descendant = _violating_pair()

        with pytest.raises(CalibrationConflictError):
            stub_resolver.validate_ancestry([ancestor, descendant])

        queried = {frozenset(call) for call in stub_resolver.tree.mrca_calls}
        assert frozenset({TRICERATOPS, ARCHAEOPTERYX}) in queried
        assert frozenset({ANKYLOSAURUS, STEGOSAURUS}) in queried


class TestAncestryGuardTopology:
    """B-8：祖先判定必须基于 MRCA 的完整后代叶集合。"""

    def test_non_overlapping_representative_clades_are_detected_as_nested(
        self, stub_resolver
    ):
        """代表类群互不重叠的嵌套节点必须被识别并报错。"""
        ancestor, descendant = _violating_pair()
        assert set(ancestor.resolved_taxa).isdisjoint(descendant.resolved_taxa)

        assert stub_resolver._is_ancestor(ancestor, descendant) is True
        with pytest.raises(CalibrationConflictError):
            stub_resolver.validate_ancestry([ancestor, descendant])

    def test_old_superset_prefilter_would_have_missed_this(self, stub_resolver):
        """回归锁：旧的 ``resolved_taxa`` 超集预过滤在这组输入上直接 return False。

        一旦有人把那道预过滤加回来，本条就会退化（``_is_ancestor`` 变 False）。
        """
        ancestor, descendant = _violating_pair()
        set1, set2 = set(ancestor.resolved_taxa), set(descendant.resolved_taxa)
        assert not set1.issuperset(set2)  # 旧代码在这一步就判"无关"
        assert stub_resolver._is_ancestor(ancestor, descendant) is True

    def test_siblings_are_not_ancestors(self, stub_resolver):
        """姐妹群（互不包含）不得被当成祖先——修 B-8 不能顺手造成大量误报。"""
        left = _point(
            "Dinosauria",
            [TRICERATOPS, TYRANNOSAURUS],
            UniformAgeConstraint(100.0, 120.0),
        )
        right = _point(
            "Stegosauria",
            [STEGOSAURUS, ANKYLOSAURUS],
            UniformAgeConstraint(180.0, 200.0),
        )

        assert stub_resolver._is_ancestor(left, right) is False
        assert stub_resolver._is_ancestor(right, left) is False
        # 两个独立节点给出互不相容的年龄是允许的（它们没有时序关系）
        assert stub_resolver.validate_ancestry([left, right]) is True

    def test_same_node_conflicting_constraints_raise(self, stub_resolver):
        """B-8 关联问题：同一节点（后代叶集合相同）给了不相容的两条约束必须拦下。

        旧实现里 ``set1 != set2`` 让同节点对直接跳过比较，与 B-6 的"同节点后者
        优先、静默覆盖"合起来两道防线全放行。
        """
        first = _point(
            "Stegosauria",
            [STEGOSAURUS, ANKYLOSAURUS],
            UniformAgeConstraint(100.0, 120.0),
        )
        second = _point(
            # 另一组代表类群，但解析到同一个节点
            "Cycadsauria",
            [ANKYLOSAURUS, STEGOSAURUS],
            UniformAgeConstraint(150.0, 180.0),
        )

        with pytest.raises(CalibrationConflictError) as excinfo:
            stub_resolver.validate_ancestry([first, second])
        assert "SAME node" in str(excinfo.value)

    def test_same_node_overlapping_constraints_pass(self, stub_resolver):
        """同节点但区间相交（软下界 + 上限）→ 合法，不得误报。"""
        first = _point(
            "Stegosauria", [STEGOSAURUS, ANKYLOSAURUS], SoftLowerBoundConstraint(60.0)
        )
        second = _point(
            "Stegosauria2",
            [ANKYLOSAURUS, STEGOSAURUS],
            MaximumAgeConstraint(150.0),
        )

        assert stub_resolver.validate_ancestry([first, second]) is True

    def test_root_calibration_is_ancestor_of_every_node(self, stub_resolver):
        """根校准（is_root_node=True）是所有节点的祖先，冲突必须被拦。"""
        root = _point(
            "Root", ALL_TIPS, UniformAgeConstraint(100.0, 150.0), is_root=True
        )
        nested = _point(
            "Stegosauria",
            [STEGOSAURUS, ANKYLOSAURUS],
            UniformAgeConstraint(180.0, 200.0),
        )

        with pytest.raises(CalibrationConflictError):
            stub_resolver.validate_ancestry([root, nested])
        with pytest.raises(CalibrationConflictError):  # 顺序无关
            stub_resolver.validate_ancestry([nested, root])

    def test_undecidable_topology_is_reported_loudly(self, stub_resolver):
        """无法定位 MRCA（幽灵分类单元）时：不得声称是祖先，但必须 error 级留痕。

        "守卫什么都没检查"不能被伪装成通过（C-31 同族的 fail-open）。
        """
        logger = Mock()
        stub_resolver.logger = logger
        ghost = _point(
            "Ghost", ["GhostX", "GhostY"], UniformAgeConstraint(100.0, 150.0)
        )
        other = _point(
            "Stegosauria",
            [STEGOSAURUS, ANKYLOSAURUS],
            UniformAgeConstraint(180.0, 200.0),
        )

        assert stub_resolver.validate_ancestry([ghost, other]) is True
        assert logger.error.called, "无法完成守卫时必须 error 级披露，不能静默放行"


# --------------------------------------------------------------------------- #
# 同一批场景用**真实** PhylogeneticTree 再跑一遍。
# 当前工作区的解释器（Python 3.14）上 ete3 无法导入、``tree.get_mrca`` 返回 None
# （审阅项 A-5/B-10），因此在没有可用拓扑后端时这一组显式 skip ——
# 用 skip 而不是"改断言迁就实现"：models/tree.py 的后端回退修好后会自动生效。
# --------------------------------------------------------------------------- #
def _topology_backend_available() -> bool:
    probe = PhylogeneticTree("((A,B),(C,D));")
    try:
        return probe.get_mrca(["A", "B"]) is not None
    except Exception:
        return False


_REQUIRES_TOPOLOGY = pytest.mark.skipif(
    not _topology_backend_available(),
    reason=(
        "PhylogeneticTree.get_mrca() 无可用拓扑后端"
        "（ete3 在本解释器上不可导入，见审阅项 A-5/B-10）"
    ),
)


@_REQUIRES_TOPOLOGY
class TestAncestryGuardRealTree:
    """用真实解析出的树验证同一批不变量（B-7 + B-8）。"""

    @pytest.fixture
    def tree(self):
        return PhylogeneticTree(NESTING_NEWICK)

    @pytest.fixture
    def resolver(self, tree):
        return CalibrationResolver(tree)

    def test_stub_topology_matches_real_parser(self, resolver):
        """测试替身的拓扑表与真实解析器一致（防止 oracle 本身写错）。"""
        for representatives in (
            [TRICERATOPS, ARCHAEOPTERYX],
            [STEGOSAURUS, ANKYLOSAURUS],
            [TRICERATOPS, TYRANNOSAURUS],
        ):
            stub = _TopologyStub(NESTING_NEWICK, ALL_TIPS, CLADES)
            assert set(resolver.tree.get_mrca(representatives)) == set(
                stub.get_mrca(representatives)
            ), representatives

    def test_is_ancestor_uses_real_topology(self, resolver):
        ancestor, descendant = _violating_pair()

        assert resolver._is_ancestor(ancestor, descendant) is True
        assert resolver._is_ancestor(descendant, ancestor) is False

    def test_violation_raises_in_both_orders(self, resolver):
        ancestor, descendant = _violating_pair()

        with pytest.raises(CalibrationConflictError):
            resolver.validate_ancestry([ancestor, descendant])
        with pytest.raises(CalibrationConflictError):
            resolver.validate_ancestry([descendant, ancestor])

    def test_consistent_pair_passes(self, resolver):
        ancestor, descendant = _consistent_pair()

        assert resolver.validate_ancestry([ancestor, descendant]) is True

    def test_fixed_age_conflict_raises(self, resolver):
        """固定点校准同样受守卫保护（祖先 150 Ma 不可能比后代 200 Ma 年轻）。"""
        ancestor = _point(
            "Archosauria", [TRICERATOPS, ARCHAEOPTERYX], FixedAgeConstraint(150.0)
        )
        descendant = _point(
            "Stegosauria", [STEGOSAURUS, ANKYLOSAURUS], FixedAgeConstraint(200.0)
        )

        with pytest.raises(CalibrationConflictError):
            resolver.validate_ancestry([descendant, ancestor])


class TestCalibrationResolverEdgeCases:
    """CalibrationResolver 边界情况测试"""

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        return tree

    @pytest.fixture
    def resolver(self, mock_tree):
        return CalibrationResolver(mock_tree)

    def test_exact_match(self, resolver):
        """测试精确匹配"""
        result = resolver._exact_match("human")
        assert result is not None

    def test_word_boundary_match(self, resolver):
        """测试词边界匹配"""
        result = resolver._word_boundary_match("human")
        assert result is not None

    def test_validate_ancestry_empty_or_single_is_noop(self, resolver):
        """空列表 / 单点列表没有可比较的对，直接通过。"""
        assert resolver.validate_ancestry([]) is True
        assert resolver.validate_ancestry(None) is True

    def test_undecidable_mock_tree_never_claims_ancestry(self, resolver):
        """树后端给不出可用拓扑时，``_is_ancestor`` 必须 False（不猜）。"""
        cal1 = _point("a", ["human", "chimp"], UniformAgeConstraint(100.0, 120.0))
        cal2 = _point("b", ["mouse", "rat"], UniformAgeConstraint(150.0, 180.0))

        assert resolver._is_ancestor(cal1, cal2) is False
        assert resolver._is_ancestor(cal2, cal1) is False

    def test_validate_ancestry_receives_one_positional_argument(self):
        """签名契约：``validate_ancestry`` 只收**一个**位置参数（校准点列表）。

        审阅项 B-9 记录旧测试写成 ``resolver.validate_ancestry(cal1, cal2)``
        → 必然 TypeError → 被测代码一行都没跑过。
        """
        import inspect

        signature = inspect.signature(CalibrationResolver.validate_ancestry)
        positional = [
            p
            for p in signature.parameters.values()
            if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
        ]
        assert [p.name for p in positional] == ["self", "calibration_points"]

    def test_skew_t_bound_uses_documented_normal_approximation(self, stub_resolver):
        """C-3：skew-t 没有精确分位数可用（scipy 从不提供 skewt 分布），

        旧代码 ``scipy.stats.skewt.ppf`` 恒抛异常、被 ``except`` 吞掉后静默返回
        ``location*0.5`` 这个凭空中点。现在必须走 constraints.py 的
        "Azzalini 精确矩 + 正态近似"，并且发 warning 说明这是近似。
        """
        constraint = SkewTConstraint(location=100.0, scale=10.0, shape=0.5, df=10.0)
        logger = Mock()
        stub_resolver.logger = logger

        lower = stub_resolver._get_min_age(constraint)
        upper = stub_resolver._get_max_age(constraint)
        mean, sd = constraint.moment_approximation()
        z = CalibrationResolver._NORMAL_975_Z

        assert lower == pytest.approx(mean - z * sd)
        assert upper == pytest.approx(mean + z * sd)
        assert lower != 100.0 * 0.5 and upper != 100.0 * 2.0
        assert logger.warning.called

    def test_skew_normal_bound_still_uses_exact_scipy_quantile(self, stub_resolver):
        """skewnorm 是 scipy 真实存在的分布，不应被一并降级成近似。"""
        constraint = SkewNormalConstraint(location=100.0, scale=10.0, shape=0.5)

        assert stub_resolver._get_max_age(constraint) > 100.0
        assert stub_resolver._get_min_age(constraint) < 100.0

    def test_skew_t_conflict_is_caught_through_validate_ancestry(self, stub_resolver):
        """先验型约束也必须参与不变量检查（此前 skewt 分支返回假中点）。"""
        ancestor = _point(
            "Archosauria",
            [TRICERATOPS, ARCHAEOPTERYX],
            FixedAgeConstraint(150.0),
        )
        descendant = _point(
            "Stegosauria",
            [STEGOSAURUS, ANKYLOSAURUS],
            SkewTConstraint(location=400.0, scale=5.0, shape=1.0, df=5.0),
        )

        with pytest.raises(CalibrationConflictError):
            stub_resolver.validate_ancestry([descendant, ancestor])


# --------------------------------------------------------------------------- #
# B-20：``mrca_leaf_pair is None`` 是**合法**输入（根校准），
# 但 auto_calibrator 的三条汇报路径都无条件写 ``mrca_leaf_pair[0]``。
# --------------------------------------------------------------------------- #
class TestAutoCalibratorRootCalibration:
    """根校准（无叶节点对）不得让汇报路径崩成 TypeError。"""

    @pytest.fixture
    def tree(self):
        return PhylogeneticTree(NESTING_NEWICK)

    @pytest.fixture
    def calibrator(self, tree):
        from phylodater.services.auto_calibrator import AutoCalibrator

        return AutoCalibrator(tree)

    @pytest.fixture
    def root_point(self):
        """报告里实测的崩溃输入：``mrca_leaf_pair=None`` + ``is_root_node=True``。"""
        return CalibrationPoint(
            name="LUCA",
            age_constraint=FixedAgeConstraint(fixed_age=3800.0),
            resolved_taxa=[],
            mrca_leaf_pair=None,
            is_root_node=True,
        )

    def test_describe_mrca_location_is_none_safe(self, calibrator):
        from phylodater.services.auto_calibrator import AutoCalibrator

        assert AutoCalibrator.describe_mrca_location(("A", "B")) == "(A, B)"
        assert AutoCalibrator.describe_mrca_location(None, True) == "(root)"
        assert AutoCalibrator.describe_mrca_location(None, False) == "(unresolved)"

    def test_summary_of_root_calibration_does_not_crash(self, calibrator, root_point):
        summary = calibrator.get_calibration_summary([root_point])

        assert "(root)" in summary
        assert "3800" in summary

    def test_write_calibrated_tree_with_root_calibration(
        self, calibrator, root_point, tree, tmp_path
    ):
        from phylodater.services.auto_calibrator import CalibratedTree

        out = tmp_path / "calibrated.nwk"
        calibrator.write_calibrated_tree(
            CalibratedTree(tree=tree, calibrations=[root_point]), out
        )

        text = out.read_text(encoding="utf-8")
        assert "MRCA (root)" in text
        assert text.strip().endswith(";")  # 树本体照旧写出

    def test_end_to_end_root_auto_calibration(self, calibrator):
        """``--auto-calibrate "ROOT:3800"`` 这条合法通路必须走得通。"""
        calibrated = calibrator.auto_calibrate(["ROOT:3800"])

        assert len(calibrated.calibrations) == 1
        point = calibrated.calibrations[0]
        assert point.is_root_node is True
        # 汇报路径（旧实现在这里 TypeError）
        assert calibrator.get_calibration_summary(calibrated.calibrations)

    def test_equal_range_becomes_fixed_point(self, calibrator):
        """C-28：``TAXON:2500-2500`` 用户想要的就是固定点，不再报内部变量名错误。"""
        parsed = calibrator.parse_calibration_string("LUCA:2500-2500")

        assert parsed.is_range is False
        assert parsed.age == 2500.0

    def test_reversed_range_still_raises(self, calibrator):
        """上下界写反不会被悄悄调转（与 §六.33 的自律一致）。"""
        parsed = calibrator.parse_calibration_string("LUCA:3000-2500")

        assert parsed.is_range is True
        with pytest.raises(ValueError):
            calibrator.resolve_calibration(parsed)

    def test_garbage_fragment_raises_instead_of_being_dropped(self):
        """C-27：不含 ``:`` 的片段以前被静默丢弃（用户以为传了两个约束）。"""
        from phylodater.services.auto_calibrator import parse_auto_calibrate_args

        assert parse_auto_calibrate_args(["Bacteria:3500", "Cyanobacteriota:2500"]) == [
            "Bacteria:3500",
            "Cyanobacteriota:2500",
        ]
        with pytest.raises(ValueError) as excinfo:
            parse_auto_calibrate_args(["Bacteria:3500 typo"])
        assert "typo" in str(excinfo.value)
