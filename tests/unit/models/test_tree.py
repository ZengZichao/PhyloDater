"""
PhylogeneticTree 模型单元测试

主要覆盖与定年方法相关的树转换方法。

这些测试**不得**要求 ete3 可导入（审阅项 A-5/B-10）：ete3 3.1.x 导入被 Python 3.13
移除的标准库 ``cgi``，因此在本工作区的解释器上根本无法导入，而校准标注写树正是
MCMCTree 通路的核心步骤。
"""

import sys

import pytest

from phylodater.core.exceptions import (
    CalibrationConflictError,
    CalibrationResolutionError,
)
from phylodater.models import CalibrationPoint, PhylogeneticTree
from phylodater.models.constraints import (
    FixedAgeConstraint,
    MaximumAgeConstraint,
    SoftBoundsConstraint,
    UniformAgeConstraint,
)

# constraints.py 的实际出口（审阅项 B-1：MCMCTree 的 B() 有 4 个槽，
# 均匀界用极小非零尾部 NEAR_HARD_TAIL_PROB=0.0001 逼近"硬界"）
UNIFORM_6_8_IN_GA = "B(0.0060, 0.0080, 0.0001, 0.0001)"


@pytest.fixture
def ete3_blocked():
    """把 ete3 伪装成"本解释器上不可导入"，并复位模型的惰性探测缓存。

    这是关键的一条环境隔离：本工作区（Python 3.14）上 ete3 恰好装不上，
    只靠"当前环境碰巧没有"来通过测试是不够的——必须在能导入 ete3 的
    解释器（3.12）上同样证明核心功能不依赖它。fixture 放在模块层，
    让所有需要屏蔽 ete3 的测试类都能复用。
    """
    import phylodater.models.tree as tree_module

    saved = (
        tree_module._ete_import_state,
        tree_module._ete_tree_class,
        tree_module._ete_import_error,
    )
    saved_module = sys.modules.get("ete3")
    sys.modules["ete3"] = None  # 使 `from ete3 import Tree` 抛 ImportError
    tree_module._ete_import_state = tree_module._ETE_UNTESTED
    tree_module._ete_tree_class = None
    tree_module._ete_import_error = None
    try:
        yield tree_module
    finally:
        (
            tree_module._ete_import_state,
            tree_module._ete_tree_class,
            tree_module._ete_import_error,
        ) = saved
        if saved_module is None:
            sys.modules.pop("ete3", None)
        else:
            sys.modules["ete3"] = saved_module


class TestNoEte3Requirement:
    """模型层不得硬依赖 ete3（A-5/B-10）"""

    def test_importing_tree_module_does_not_import_ete3(self, module_imports_cleanly):
        """``phylodater.models.tree`` 被导入时不应把 ete3 拉进进程。

        探针跑在一个新进程里：本批测试中先跑的可视化测试会真正导入 ete3，
        直接断言当前进程的 ``sys.modules`` 只会测出“测试顺序”，而不是产品行为。
        """
        assert "phylodater.models.tree" in sys.modules
        assert module_imports_cleanly(
            "phylodater.models.tree", "ete3"
        ), "模型层在模块导入阶段就拉起了 ete3；ete3 必须是真正可选的后端"

    def test_every_topology_op_works_without_ete3(self, ete3_blocked):
        """MRCA / 单系性 / 定根 / 校准标注在 ete3 彻底不可导入时全部可用。"""
        assert ete3_blocked.ete3_available() is False
        tree = PhylogeneticTree("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")
        cal = CalibrationPoint(
            name="Primates",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        assert tree.get_mrca(["human", "chimp"]) == ["chimp", "human"]
        assert tree.is_monophyletic(["human", "chimp"]) is True
        assert tree.is_monophyletic(["human", "mouse"]) is False
        assert tree.get_mrca_terminals_all(["human", "chimp"]) == ["chimp", "human"]
        assert "B(" in tree.with_calibration_annotations([cal]).newick
        assert tree.with_paml48_calibration_annotations([cal]).newick.startswith(
            "4 1\n"
        )
        assert tree.rerooted_at_midpoint().newick.count(":") == 6
        assert tree.rerooted_with_outgroup(["mouse"])[0].tip_names.count("mouse") == 1
        assert tree.rerooted_at_min_variance()[0].tip_names == tree.tip_names

    def test_ete3_backend_raises_import_error_when_unavailable(self, ete3_blocked):
        """ete3 不可用时，可视化层依赖的后端入口必须抛 ImportError（可被捕获降级）。"""
        tree = PhylogeneticTree("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")

        with pytest.raises(ImportError, match="ete3"):
            tree._ete3_tree
        with pytest.raises(ImportError, match="ete3"):
            tree._ete4_tree


class TestCalibrationAnnotations:
    """测试树校准标注方法"""

    @pytest.fixture
    def simple_tree(self):
        return PhylogeneticTree("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")

    def test_with_calibration_annotations_removes_branch_lengths(self, simple_tree):
        """with_calibration_annotations 应输出无分支长度的树"""
        cal = CalibrationPoint(
            name="Primates",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        annotated = simple_tree.with_calibration_annotations([cal])

        assert ":" not in annotated.newick
        assert "human" in annotated.newick
        assert "chimp" in annotated.newick

    def test_with_calibration_annotations_includes_constraint(self, simple_tree):
        """约束应以内联节点标注形式出现"""
        cal = CalibrationPoint(
            name="Primates",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        annotated = simple_tree.with_calibration_annotations([cal])

        # FixedAgeConstraint(70.0 Ma) -> constraints.py 的 4 槽窄软界（Ga 单位）
        assert "B(" in annotated.newick
        assert "B(0.0693, 0.0707, 0.005, 0.005)" in annotated.newick

    def test_with_calibration_annotations_skips_root(self, simple_tree):
        """根节点约束不应写入树文件"""
        root_cal = CalibrationPoint(
            name="Root",
            is_root_node=True,
            age_constraint=MaximumAgeConstraint(max_age=100.0),
        )
        internal_cal = CalibrationPoint(
            name="Primates",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=UniformAgeConstraint(min_age=6.0, max_age=8.0),
        )

        annotated = simple_tree.with_calibration_annotations([root_cal, internal_cal])

        # 根约束不应出现在树标注中
        assert "U(" not in annotated.newick
        # 内部节点约束应出现（串由 constraints.py 的公开渲染器给出）
        assert UNIFORM_6_8_IN_GA in annotated.newick

    def test_constraint_string_is_what_constraints_py_emits(self, simple_tree):
        """写进树的串必须逐字等于 constraints.py 的出口，模型层不得自己拼 B()。

        钉住"集中式渲染"这条链路：一旦 models/tree.py 里重新出现手写的
        ``.4f`` Ma→Ga 换算（审阅项 A-3 的病因），本测试即失败。
        """
        for constraint in (
            FixedAgeConstraint(fixed_age=70.0),
            UniformAgeConstraint(min_age=6.0, max_age=8.0),
            SoftBoundsConstraint(min_age=60.0, max_age=90.0),
        ):
            cal = CalibrationPoint(
                name="n", mrca_leaf_pair=("human", "chimp"), age_constraint=constraint
            )
            expected = constraint.to_mcmctree_calib_string()
            assert expected in simple_tree.with_calibration_annotations([cal]).newick
            assert expected, "渲染器不应产出空串"

    def test_fixed_age_below_5_ma_does_not_collapse(self, simple_tree):
        """0.1 Ma 级别的固定校准不得被格式化压成零宽度（审阅项 A-3）。"""
        cal = CalibrationPoint(
            name="tiny",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=FixedAgeConstraint(fixed_age=0.1),
        )

        annotated = simple_tree.with_calibration_annotations([cal])

        assert "B(0.000099, 0.000101, 0.005, 0.005)" in annotated.newick
        assert "B(0.0000, 0.0000" not in annotated.newick

    def test_with_paml48_calibration_annotations_includes_header(self, simple_tree):
        """PAML 4.8 格式应自带 PHYLIP 头"""
        cal = CalibrationPoint(
            name="Primates",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        annotated = simple_tree.with_paml48_calibration_annotations([cal])

        lines = annotated.newick.split("\n")
        assert lines[0] == "4 1"
        assert "B(" in annotated.newick

    def test_paml48_root_constraint_uses_central_ga_renderer(self, simple_tree):
        """4.8 的根上界后缀 '<…' 必须走 constraints.py 的 Ga 渲染，不得 .4f 坍缩。"""
        root_cal = CalibrationPoint(
            name="Root",
            is_root_node=True,
            age_constraint=MaximumAgeConstraint(max_age=100.0),
        )
        tiny_root = CalibrationPoint(
            name="Root",
            is_root_node=True,
            age_constraint=MaximumAgeConstraint(max_age=3.0),
        )
        internal = CalibrationPoint(
            name="Primates",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        assert (
            simple_tree.with_paml48_calibration_annotations([internal, root_cal])
            .newick.rstrip()
            .endswith("'<0.1000'")
        )
        # 3 Ma → 0.003 Ga：4 位小数仍够，但绝不能是 0.0000
        assert (
            "<0.0000'"
            not in simple_tree.with_paml48_calibration_annotations(
                [internal, tiny_root]
            ).newick
        )

    def test_both_annotation_methods_are_equivalent_without_header(self, simple_tree):
        """除 PHYLIP 头与空格外，两种标注方法生成的 Newick 应等价"""
        cal = CalibrationPoint(
            name="Primates",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        t_4108 = simple_tree.with_calibration_annotations([cal])
        t_48_full = simple_tree.with_paml48_calibration_annotations([cal])
        t_48_body = "\n".join(t_48_full.newick.split("\n")[1:])

        # 4.10.8+ 格式在逗号后保留空格，4.8 格式不带空格；语义等价即可
        assert t_4108.newick.replace(" ", "") == t_48_body.replace(" ", "")

    def test_multiple_calibrations_land_on_their_own_nodes(self, simple_tree):
        """两条落在不同节点的校准都必须出现（一条都不能丢）。"""
        primates = CalibrationPoint(
            name="Primates",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )
        rodentia = CalibrationPoint(
            name="Rodentia",
            mrca_leaf_pair=("mouse", "rat"),
            age_constraint=UniformAgeConstraint(min_age=6.0, max_age=8.0),
        )

        annotated = simple_tree.with_calibration_annotations(
            [primates, rodentia]
        ).newick

        assert annotated.count("B(") == 2
        assert "B(0.0693, 0.0707, 0.005, 0.005)" in annotated
        assert UNIFORM_6_8_IN_GA in annotated

    # ---------------- B-6：静默丢点 / 同节点覆盖 ---------------- #

    def test_missing_mrca_pair_raises_instead_of_silent_drop(self, simple_tree):
        """非根校准缺 mrca_leaf_pair 时必须报错，而不是静默 continue（B-6）。"""
        cal = CalibrationPoint(
            name="NoPair",
            resolved_taxa=["human", "chimp"],
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        with pytest.raises(CalibrationResolutionError, match="mrca_leaf_pair"):
            simple_tree.with_calibration_annotations([cal])

    def test_unlocatable_mrca_pair_raises_instead_of_warning_only(self, simple_tree):
        """叶名不在树中时同样必须报错（此前只留一条 warning 然后丢点）。"""
        cal = CalibrationPoint(
            name="Ghost",
            mrca_leaf_pair=("human", "not_in_tree"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        with pytest.raises(CalibrationResolutionError, match="无法在树中定位"):
            simple_tree.with_calibration_annotations([cal])

    def test_leaf_level_mrca_pair_raises(self, simple_tree):
        """两个代表类群定位到叶节点（而非内部节点）时不能写出无效标注。"""
        cal = CalibrationPoint(
            name="SameLeaf",
            mrca_leaf_pair=("human", "human"),
            age_constraint=FixedAgeConstraint(fixed_age=70.0),
        )

        with pytest.raises(CalibrationResolutionError, match="而不是内部节点"):
            simple_tree.with_calibration_annotations([cal])

    def test_same_node_calibrations_are_merged_by_intersection(self, simple_tree):
        """同一节点上的多条校准：区间求交（下界取最大、上界取最小），而非后者覆盖。"""
        wide = CalibrationPoint(
            name="Wide",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=UniformAgeConstraint(min_age=60.0, max_age=100.0),
        )
        narrow = CalibrationPoint(
            name="Narrow",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=UniformAgeConstraint(min_age=70.0, max_age=80.0),
        )

        annotated = simple_tree.with_calibration_annotations([wide, narrow]).newick

        assert "B(0.0700, 0.0800, 0.0001, 0.0001)" in annotated
        # 任何一条原始区间都不允许"赢了"（那正是静默覆盖的形态）
        assert "B(0.0600, 0.1000" not in annotated
        assert "B(0.0700, 0.0800, 0.0001, 0.0001)" in annotated
        assert annotated.count("B(") == 1

    def test_same_node_detection_uses_topology_not_leaf_pair(self, simple_tree):
        """叶对不同但节点相同（嵌套分类群）也要被识别为同一条枝。

        ``(human,chimp)`` 与 ``(chimp,human)`` 是同一节点的两种锚定写法。
        """
        first = CalibrationPoint(
            name="A",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=UniformAgeConstraint(min_age=60.0, max_age=100.0),
        )
        second = CalibrationPoint(
            name="B",
            mrca_leaf_pair=("chimp", "human"),
            age_constraint=UniformAgeConstraint(min_age=50.0, max_age=70.0),
        )

        annotated = simple_tree.with_calibration_annotations([first, second]).newick

        assert annotated.count("B(") == 1
        assert "B(0.0600, 0.0700, 0.0001, 0.0001)" in annotated

    def test_empty_intersection_raises_conflict(self, simple_tree):
        """交集为空（互相矛盾的先验）必须报错——MCMCTree 也会以 'fossil bounds
        in tree incorrect' 中止。"""
        old = CalibrationPoint(
            name="Old",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=UniformAgeConstraint(min_age=200.0, max_age=300.0),
        )
        young = CalibrationPoint(
            name="Young",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=UniformAgeConstraint(min_age=10.0, max_age=20.0),
        )

        with pytest.raises(CalibrationConflictError, match="交集为空"):
            simple_tree.with_calibration_annotations([old, young])

    def test_same_node_merge_is_recorded_as_warning(self, simple_tree, capsys):
        """合并必须留痕：至少一条 warning 说明"曾经会是静默覆盖"。"""
        first = CalibrationPoint(
            name="A",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=UniformAgeConstraint(min_age=60.0, max_age=100.0),
        )
        second = CalibrationPoint(
            name="B",
            mrca_leaf_pair=("human", "chimp"),
            age_constraint=UniformAgeConstraint(min_age=70.0, max_age=80.0),
        )

        simple_tree.with_calibration_annotations([first, second])

        captured = capsys.readouterr()
        messages = captured.out + captured.err
        assert "WARNING" in messages, "同节点合并只留 info 级是不够的"
        assert "SAME node" in messages


class TestDistanceMatrixCache:
    """测试距离矩阵缓存（#15 性能修复）"""

    @pytest.fixture
    def tree(self):
        return PhylogeneticTree("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")

    def test_get_distance_matrix_is_deterministic(self, tree):
        """多次调用应返回一致的矩阵（缓存命中后结果不变）"""
        m1 = tree.get_distance_matrix()
        m2 = tree.get_distance_matrix()

        assert m1 == m2

    def test_get_distance_matrix_uses_cache(self, tree):
        """重复调用应命中实例级缓存，而非每次重建矩阵。

        通过对比缓存字典在第二次调用后大小不变来验证（无新增计算产物）。
        """
        first = tree.get_distance_matrix()
        cache_after_first = dict(tree._distance_matrix_cache)

        second = tree.get_distance_matrix()
        cache_after_second = dict(tree._distance_matrix_cache)

        # 缓存条目数量不应因重复调用而增长
        assert len(cache_after_first) == len(cache_after_second)
        assert first == second

    def test_get_distance_matrix_values_symmetric(self, tree):
        """距离矩阵应对称：dist(a,b) == dist(b,a) 且无负距离。

        矩阵为 {(tip1, tip2): dist} 扁平字典，包含全部对称条目。
        """
        matrix = tree.get_distance_matrix()
        tips = set()
        for a, b in matrix:
            tips.add(a)
            tips.add(b)

        assert len(matrix) > 0  # 非空（4 个 tip -> 12 条对称距离）
        for a in tips:
            for b in tips:
                if a == b:
                    continue
                assert matrix[(a, b)] == matrix[(b, a)]
                assert matrix[(a, b)] >= 0.0


class TestLeadingCommentTolerance:
    """树字符串开头的 ``[&R]`` / ``[&U]`` / NHX 注释不得导致解析失败。

    BEAST、IQ-TREE、DendroPy 等主流工具写出的树常常以有根标记开头（``[&R] (…)``）。
    旧的最小游走器在这里只读标签，遇到 ``[`` 就停，于是整棵树被解成一个没有叶子的
    根：MCMCTree 与 treePL 通路直接报 "Cannot parse Newick for tree operations"，
    祖先关系校验也退化为"无法解析 MRCA 叶集合"。ML 分支长度变体的基准测试暴露了这个缺陷。
    """

    CASES = [
        (
            "[&R] (Dog:0.09,((Rabbit:0.10,(Mouse:0.02,Rat:0.01):0.09):0.03,Human:0.01):0.09);",
            ["Dog", "Human", "Mouse", "Rabbit", "Rat"],
        ),
        ("[&&R] (A:1,(B:1,C:1)BC:1);", ["A", "B", "C"]),
        ("[&U] (A:1,B:1);", ["A", "B"]),
        ("[&nhx key=1] ((A:1,B:1):2,C:3);", ["A", "B", "C"]),
    ]

    @pytest.mark.parametrize("newick,expected_tips", CASES)
    def test_tips_survive_a_leading_comment(self, newick, expected_tips):
        tree = PhylogeneticTree(newick)

        assert sorted(tree.tip_names) == sorted(expected_tips)
        # 拓扑操作（而不是只取叶名）也必须可用
        assert tree.iter_clades()
        assert tree.get_descendant_leaves(expected_tips[:2]) is not None

    def test_rerooting_works_on_a_commented_tree(self):
        tree = PhylogeneticTree("[&R] ((A:1,B:1):1,(C:1,D:1):1);")

        rerooted = tree.rerooted_with_outgroup(["C"])[0]
        assert sorted(rerooted.tip_names) == ["A", "B", "C", "D"]


class TestCladeApi:
    """公开的"分支/单系群"枚举 API（供比较层退役其本地 Newick 解析器）

    跨文件需求：多方法比较要按**拓扑**对齐节点，而 ``PhylogeneticTree`` 的旧公开
    API 里没有"枚举所有分支节点及其后代叶集合"这一件，于是
    ``core/comparison_reporter.py`` 自带了一个最小 Newick 读取器。这里补上，并
    顺便修正 ``get_node_topology_id`` 的语义（此前它只是把入参排序拼接，因此
    锚定同一节点的两种叶对写法拿不到同一个 ID）。
    """

    @pytest.fixture
    def tree(self):
        return PhylogeneticTree(
            "((human:0.1,chimp:0.1)Primates:0.2,(mouse:0.3,rat:0.3)Rodentia:0.1);"
        )

    @pytest.fixture
    def ladder_tree(self):
        return PhylogeneticTree("(((A:1,B:1):1,C:1):1,(D:1,E:1):1);")

    def test_iter_clades_enumerates_every_node_with_its_full_leaf_set(self, tree):
        clades = tree.iter_clades()

        # 1 根 + 2 内部 + 4 叶
        assert len(clades) == 7
        assert all(
            set(c) == {"label", "leaves", "is_leaf", "depth", "topology_id"}
            for c in clades
        )
        by_label = {c["label"]: c for c in clades}
        assert by_label["Primates"]["leaves"] == ("human", "chimp")
        assert by_label["Primates"]["depth"] == 1
        assert by_label["Primates"]["is_leaf"] is False
        assert by_label["human"]["depth"] == 2
        assert by_label["human"]["is_leaf"] is True
        assert by_label[""]["leaves"] == ("human", "chimp", "mouse", "rat")
        assert by_label[""]["depth"] == 0

    def test_iter_clades_is_preorder_with_monotone_depth(self, tree):
        """前序：父先于子，深度沿路径单调不减（报告层按此自顶向下切分子树）。"""
        clades = tree.iter_clades()
        positions = {c["label"]: i for i, c in enumerate(clades)}

        assert positions[""] < positions["Primates"]
        assert positions["Primates"] < positions["human"]
        assert [c["depth"] for c in clades[:3]] == [0, 1, 2]

    def test_topology_id_is_node_identity_not_a_reformatted_query(self, ladder_tree):
        """锚定同一节点的不同叶对 -> 同一个 ID（这是"修正语义"的核心断言）。

        在 ``(((A,B),C),(D,E))`` 上，叶对 (A,C) 与 (B,C) 分属节点 ``((A,B),C)`` 的
        两个不同子树，因此它们的 MRCA 都是那个节点——旧实现按入参排序拼接，
        会得到 "A,C" 与 "B,C" 两个不同 ID，同一节点因此无法跨校准/跨方法对齐。
        """
        ladder = ladder_tree

        first = ladder.get_node_topology_id(["A", "C"])
        second = ladder.get_node_topology_id(["B", "C"])
        assert first == second == ladder.get_node_topology_id(["A", "B", "C"])
        # (A,B) 锚定的是更深的那个节点 ((A,B))，必须是另一个 ID
        assert first != ladder.get_node_topology_id(["A", "B"])
        # 另一侧节点也必须不同
        assert first != ladder.get_node_topology_id(["D", "E"])
        # 旧实现的形态（入参排序拼接）会让前两者互不相同：
        assert first != "A,C"

    def test_topology_id_marks_unlocatable_queries(self, tree):
        """无法按拓扑定位时 ID 带 '?' 前缀，不伪装成正常 ID。"""
        assert tree.get_node_topology_id(["ghost"]).startswith("?")
        assert tree.get_node_topology_id(["human", "ghost"]).startswith("?")

    def test_get_descendant_leaves_matches_get_mrca(self, tree):
        assert tree.get_descendant_leaves(["human", "chimp"]) == ["chimp", "human"]
        assert tree.get_descendant_leaves(["ghost"]) is None
        assert tree.get_descendant_leaves([]) == []

    def test_get_node_depth_counts_edges_from_root(self, tree):
        assert tree.get_node_depth(["human", "chimp"]) == 1
        assert tree.get_node_depth(["mouse", "rat"]) == 1
        # 跨根两侧的叶对，其 MRCA 就是根 -> 深度 0
        assert tree.get_node_depth(["human", "mouse"]) == 0
        assert tree.get_node_depth(["ghost"]) == 0
        assert tree.get_node_depth([]) == 0

    def test_clade_api_survives_unparseable_input(self):
        broken = PhylogeneticTree("((A,B")

        assert broken.iter_clades() == []
        assert broken.get_descendant_leaves(["A"]) is None
        assert broken.get_node_depth(["A"]) == 0

    def test_clade_api_needs_no_ete3(self, ete3_blocked):
        """该 API 只用零依赖游走器（在 ete3 不可导入时仍然可用）。

        不断言当前进程的 ``sys.modules``：同一批运行里其他测试可以合法地导入
        ete3；真正要证明的是“ete3 导不进来时该 API 仍然跑”，所以直接屏蔽它。
        """
        assert ete3_blocked.ete3_available() is False
        tree = PhylogeneticTree("((A:1,B:1):2,(C:1,D:1):3);")

        assert len(tree.iter_clades()) == 7
