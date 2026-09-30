"""
TreeRootingService 单元测试

测试树定根服务的各项功能。

这些测试**不 mock ete3**（审阅项 A-5/B-10）：ete3 在 Python >= 3.13 上无法导入，
而旧测试靠 ``patch("ete3.Tree")`` 才跑得起来，等于"7 个测试证明不了任何真实行为"。
现在直接跑真实拓扑，断言的是定根操作的**不变量**：叶集合不变、根确实是二叉、
根的分裂确实是所请求的那个外群、tip-to-tip 距离不被重根改变。

MAD 一组是 C-32 的回归：输出文件名显式交给外部程序，拿不到时靠**本次专用临时目录**
确定性定位，实在没有才回退 stdout 且先过 Newick 形状检查；"没有输出文件"、"文件里
不是树"、"是树但不是这棵树"三条失败路径各给各的错误信息。
"""

import shlex
import stat
from pathlib import Path

import pytest

from phylodater.core.exceptions import TreeValidationError
from phylodater.models import PhylogeneticTree
from phylodater.services.tree_rooting import (
    RootingMethod,
    RootingResult,
    TreeRootingService,
    _looks_like_newick,
)

ROOTED_NWK = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
UNROOTED_NWK = "(human:0.1,chimp:0.1,mouse:0.3,rat:0.3);"
STAR_NWK = "(A:1,B:1,C:1,D:1);"


def _tip_to_tip(tree: PhylogeneticTree, tip_a: str, tip_b: str) -> float:
    """两叶之间的路径长度（DendroPy 的超度量无关距离，作为独立参照实现）。"""
    matrix = tree.get_distance_matrix()
    return matrix[(tip_a, tip_b)]


def _root_heights(tree: PhylogeneticTree) -> dict:
    """根到各叶的距离——用 DendroPy 独立算一遍，交叉校验本模块的 Newick 游走器。"""
    from dendropy import Tree

    parsed = Tree.get(
        data=tree.newick,
        schema="newick",
        preserve_underscores=True,
        rooting="force-rooted",
    )
    return {
        leaf.taxon.label: float(leaf.distance_from_root())
        for leaf in parsed.leaf_node_iter()
    }


def _variance(values) -> float:
    values = list(values)
    mean = sum(values) / len(values)
    return sum((value - mean) ** 2 for value in values) / len(values)


@pytest.fixture
def rerooter():
    """创建定根服务实例"""
    return TreeRootingService()


@pytest.fixture
def tree():
    return PhylogeneticTree(ROOTED_NWK)


@pytest.fixture
def temp_dir(tmp_path):
    """临时输出目录"""
    return tmp_path


class TestNoEte3InRootingService:
    """定根服务不得依赖 ete3（A-5/B-10）"""

    def test_module_never_imports_ete3(self, module_imports_cleanly):
        from pathlib import Path

        import phylodater.services.tree_rooting as module

        assert module_imports_cleanly(
            "phylodater.services.tree_rooting", "ete3"
        ), "定根服务在模块导入阶段就拉起了 ete3"
        source = Path(module.__file__).read_text(encoding="utf-8")
        assert "from ete3" not in source
        assert "import ete3" not in source

    def test_is_rooted_fails_closed_when_the_tree_is_unreadable(self, rerooter):
        """C-31：旧实现在探测失败（ete3 ``ImportError``）时 ``return True``。

        那种 fail-open 让任何用 ``is_rooted()`` 做决策的调用方都得到"已经定过根了"，
        与 B-22 的"看似已定根就跳过"叠加后重根会被静默跳过。现在看不懂一律答
        "未定根"，由上层按 ``require_rooted`` 决定报错还是告警。
        """

        class _Junk:
            newick = "this is not a tree at all!!"

        assert rerooter.is_rooted(PhylogeneticTree(ROOTED_NWK)) is True
        assert rerooter.is_rooted(_Junk()) is False

    def test_all_methods_run_on_this_interpreter(self):
        """本解释器上（ete3 装不上）三种内置定根方法都必须真能跑。"""
        service = TreeRootingService()
        tree = PhylogeneticTree(ROOTED_NWK)

        for method in (
            RootingMethod.MIDPOINT,
            RootingMethod.MINVAR,
            RootingMethod.OUTGROUP,
        ):
            kwargs = {}
            if method is RootingMethod.OUTGROUP:
                kwargs["outgroup_taxa"] = ["mouse"]
            result = service.root_tree(tree, method=method, **kwargs)
            assert isinstance(result, RootingResult)
            assert sorted(result.rooted_tree.tip_names) == sorted(tree.tip_names)


class TestTreeRootingService:
    """TreeRootingService 测试类（fixtures 在模块级：定根服务与树对各类共用）"""

    # ---------------- 中点定根 ---------------- #

    def test_root_tree_midpoint(self, rerooter, tree):
        """测试中点定根"""
        result = rerooter.root_tree(tree, method=RootingMethod.MIDPOINT)

        assert isinstance(result, RootingResult)
        assert result.method == RootingMethod.MIDPOINT
        rooted = result.rooted_tree
        # 叶不变
        assert sorted(rooted.tip_names) == ["chimp", "human", "mouse", "rat"]
        # 结果真的是有根二叉
        assert rooted.is_rooted is True
        # 最长路径 human↔rat = 0.7，中点 0.35 落在根上
        assert rooted.is_monophyletic(["mouse", "rat"]) is True
        assert rooted.is_monophyletic(["human", "chimp"]) is True

    def test_midpoint_does_not_skip_binary_rooted_input(self, rerooter, tree):
        """B-22 回归：输入以二叉形式书写（"看起来已有根"）也必须真的被重根。

        旧实现第一件事就是 ``if len(t.get_children()) == 2: return tree``，
        于是用户显式请求中点定根却什么也没发生，只留一条 info 日志。
        """
        assert rerooter.is_rooted(tree) is True  # 前提：它"看起来"已经有根

        result = rerooter.root_tree(tree, method=RootingMethod.MIDPOINT)

        assert result.rooted_tree is not tree
        assert result.was_rerooted is True
        assert result.skipped_reason is None
        # 旧根把树分成 {human,chimp} | {mouse,rat}；中点把根挪进了鼠枝内部，
        # 于是 {mouse,rat} 这一侧的枝被劈成两段（0.05 + 0.25）。
        assert rerooter.get_root_children(tree) == [
            ["human", "chimp"],
            ["mouse", "rat"],
        ]
        assert sorted(
            sorted(child) for child in rerooter.get_root_children(result.rooted_tree)
        ) == [sorted(["human", "chimp"]), sorted(["mouse", "rat"])]
        assert result.rooted_tree.newick != tree.newick

    def test_midpoint_preserves_all_pairwise_distances(self, rerooter, tree):
        """重根是纯拓扑操作：任意两叶距离必须逐字不变。"""
        result = rerooter.root_tree(tree, method=RootingMethod.MIDPOINT)
        rooted = result.rooted_tree

        for tip_a, tip_b in [
            ("human", "chimp"),
            ("mouse", "rat"),
            ("human", "rat"),
            ("chimp", "mouse"),
        ]:
            assert _tip_to_tip(rooted, tip_a, tip_b) == pytest.approx(
                _tip_to_tip(tree, tip_a, tip_b)
            )

    def test_midpoint_root_equalizes_the_diameter_ends(self, rerooter):
        """中点定根的定义性判据（用 DendroPy 独立量一遍根高）。

        ``((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1)`` 的最长 tip-to-tip
        路径是 0.7（如 human↔rat），因此新根到两端都得是 0.35；旧实现根本不看枝长
        就直接返回输入（输入根高是 0.3 / 0.4，并不等距）。
        """
        result = rerooter.root_tree(
            PhylogeneticTree(ROOTED_NWK), method=RootingMethod.MIDPOINT
        )
        heights = _root_heights(result.rooted_tree)

        assert max(heights.values()) == pytest.approx(0.35)
        assert sorted(heights.values()) == pytest.approx(
            sorted([0.35, 0.35, 0.35, 0.35])
        )
        # 原来的根并不等距（0.3 / 0.3 / 0.4 / 0.4），这就是"被跳过"的证据
        original_heights = _root_heights(PhylogeneticTree(ROOTED_NWK))
        assert max(original_heights.values()) == pytest.approx(0.4)
        assert _variance(original_heights.values()) > _variance(heights.values())

    def test_midpoint_splits_edge_so_both_halves_are_equal(self, rerooter):
        """直径中点落在枝内部时要把那一条枝劈成两段。"""
        tree = PhylogeneticTree("((A:1,B:6):2,(C:1,D:1):2);")

        rooted = rerooter.root_tree(tree, method=RootingMethod.MIDPOINT).rooted_tree

        heights = _root_heights(rooted)
        # 最长路径 B↔C(或 B↔D) = 6+2+2+1 = 11 ⇒ 根高处处为 5.5
        assert heights["B"] == pytest.approx(5.5)
        assert heights["C"] == pytest.approx(5.5)
        assert heights["D"] == pytest.approx(5.5)
        assert heights["A"] == pytest.approx(1.5)
        # 距离不变
        assert _tip_to_tip(rooted, "B", "D") == pytest.approx(11.0)

    def test_midpoint_zero_branch_lengths_is_not_silently_claimed(self, rerooter):
        """枝长全缺失时中点不唯一：必须显式披露"未重根"，不能假装重过。"""
        tree = PhylogeneticTree("((A,B),(C,D));")

        result = rerooter.root_tree(tree, method=RootingMethod.MIDPOINT)

        assert result.was_rerooted is False
        assert result.skipped_reason
        assert result.rooted_tree is tree

    def test_midpoint_writes_output_file(self, rerooter, tree, temp_dir):
        out = temp_dir / "mid.nwk"

        result = rerooter.root_tree(
            tree, method=RootingMethod.MIDPOINT, output_path=out
        )

        assert out.exists()
        assert out.read_text(encoding="utf-8").strip() == result.rooted_tree.newick

    # ---------------- 外群定根 ---------------- #

    def test_root_tree_outgroup(self, rerooter, tree):
        """测试外群定根"""
        result = rerooter.root_tree(
            tree, method=RootingMethod.OUTGROUP, outgroup_taxa=["mouse"]
        )

        assert isinstance(result, RootingResult)
        assert result.method == RootingMethod.OUTGROUP
        assert result.outgroup_taxa == ["mouse"]
        rooted = result.rooted_tree
        assert sorted(rooted.tip_names) == sorted(tree.tip_names)
        assert rooted.is_rooted is True
        # 外群必须是根的两个子节点之一
        assert ["mouse"] in rerooter.get_root_children(rooted)

    def test_outgroup_forces_reroot_on_binary_tree(self, rerooter, tree):
        """显式指定外群时即使"看起来已定根"也必须重挂到指定外群。"""
        result = rerooter.root_tree(
            tree, method=RootingMethod.OUTGROUP, outgroup_taxa=["rat"]
        )

        assert result.was_rerooted is True
        assert ["rat"] in rerooter.get_root_children(result.rooted_tree)
        for pair in (("human", "chimp"), ("mouse", "rat")):
            assert _tip_to_tip(result.rooted_tree, *pair) == pytest.approx(
                _tip_to_tip(tree, *pair)
            )

    def test_outgroup_multiple_taxa_use_their_mrca(self, rerooter, tree):
        """多个外群：按它们的 MRCA 定根。"""
        result = rerooter.root_tree(
            tree, method=RootingMethod.OUTGROUP, outgroup_taxa=["mouse", "rat"]
        )

        assert result.rooted_tree.is_monophyletic(["mouse", "rat"]) is True
        assert ["mouse", "rat"] in [
            sorted(child) for child in rerooter.get_root_children(result.rooted_tree)
        ]

    def test_outgroup_on_already_correct_root_reports_no_change(self, rerooter, tree):
        """树本来就在指定外群所在的枝上：如实报告未改变，而不是静默。"""
        result = rerooter.root_tree(
            tree, method=RootingMethod.OUTGROUP, outgroup_taxa=["mouse", "rat"]
        )

        assert result.rooted_tree.is_rooted is True
        # 根分裂与输入一致 -> was_rerooted 必须为 False（但仍按请求做了计算）
        assert result.was_rerooted is False

    def test_root_tree_outgroup_not_found(self, rerooter, tree):
        """测试外群不存在时抛出异常"""
        with pytest.raises(TreeValidationError, match="None of the outgroup taxa"):
            rerooter.root_tree(
                tree,
                method=RootingMethod.OUTGROUP,
                outgroup_taxa=["nonexistent"],
            )

    def test_outgroup_covering_whole_tree_reports_no_change(self, rerooter):
        """外群覆盖整棵树时无法定根：显式披露，不返回一棵伪重根的树。"""
        tree = PhylogeneticTree(UNROOTED_NWK)

        result = rerooter.root_tree(
            tree,
            method=RootingMethod.OUTGROUP,
            outgroup_taxa=["human", "chimp", "mouse", "rat"],
        )

        assert result.was_rerooted is False
        assert result.rooted_tree is tree
        assert "覆盖整棵树" in result.skipped_reason

    def test_root_tree_outgroup_missing_taxa(self, rerooter, tree):
        """测试外群定根缺少外群分类单元时抛出异常"""
        with pytest.raises(ValueError, match="Outgroup taxa required"):
            rerooter.root_tree(tree, method=RootingMethod.OUTGROUP)

    def test_root_tree_unknown_method(self, rerooter, tree):
        """测试未知定根方法抛出异常"""
        with pytest.raises(ValueError, match="Unknown rooting method"):
            rerooter.root_tree(tree, method="invalid")

    def test_outgroup_writes_output_file(self, rerooter, tree, temp_dir):
        out = temp_dir / "og.nwk"

        rerooter.root_tree(
            tree,
            method=RootingMethod.OUTGROUP,
            outgroup_taxa=["mouse"],
            output_path=out,
        )

        assert out.exists()

    # ---------------- 最小方差定根 ---------------- #

    def test_minvar_reroots_binary_tree(self, rerooter):
        """B-22 回归：minvar 同样不得因"看起来已有根"而跳过。

        旧实现在 ``len(t.get_children()) == 2`` 时直接原样返回输入树。
        ``rooted`` 在这里是**计算之后**才发现现有根已经是最优位置，
        与"根本没算"是两件事——所以用一棵根放歪了的二叉树来证明它真的算了。
        """
        tree = PhylogeneticTree("(A:1,(B:1,(C:1,D:1):1):1);")
        assert rerooter.is_rooted(tree) is True  # 前提：二叉，旧实现会直接跳过

        result = rerooter.root_tree(tree, method=RootingMethod.MINVAR)
        rooted = result.rooted_tree

        assert result.was_rerooted is True
        assert rooted.newick != tree.newick
        assert sorted(rooted.tip_names) == ["A", "B", "C", "D"]
        # 根到叶距离方差必须由 0.6875 降到 0.1875
        assert _variance(_root_heights(tree).values()) == pytest.approx(0.6875)
        assert _variance(_root_heights(rooted).values()) == pytest.approx(0.1875)
        # 距离不变
        for pair in (("A", "C"), ("B", "D"), ("A", "B")):
            assert _tip_to_tip(rooted, *pair) == pytest.approx(_tip_to_tip(tree, *pair))

    def test_minvar_is_never_worse_than_the_input_root(self, rerooter):
        """最小方差判据本身：结果方差不高于输入（含"输入已最优"的恒等情形）。"""
        cases = [
            ROOTED_NWK,  # 现有根已是最优 -> 方差相同
            "((A:1,B:6):2,(C:1,D:1):2);",  # 现有根明显偏
            "(((A:1,B:1):1,C:1):1,D:1);",
        ]
        for newick in cases:
            tree = PhylogeneticTree(newick)
            rooted = rerooter.root_tree(tree, method=RootingMethod.MINVAR).rooted_tree
            assert (
                _variance(_root_heights(rooted).values())
                <= _variance(_root_heights(tree).values()) + 1e-12
            ), newick

    def test_minvar_reports_already_optimal_root(self, rerooter):
        """现有根已是最优时，如实报告未改变（而不是静默、也不是假装重了根）。"""
        result = rerooter.root_tree(
            PhylogeneticTree(ROOTED_NWK), method=RootingMethod.MINVAR
        )

        assert result.was_rerooted is False
        assert result.skipped_reason, "根未移动时必须给出理由，不能只留布尔位"
        assert sorted(result.rooted_tree.tip_names) == [
            "chimp",
            "human",
            "mouse",
            "rat",
        ]

    def test_minvar_needs_at_least_two_leaves(self, rerooter):
        """叶少于 2 时无从计算方差：显式披露"未重根"。"""
        result = rerooter.root_tree(
            PhylogeneticTree("A:1;"), method=RootingMethod.MINVAR
        )

        assert result.was_rerooted is False
        assert "叶节点少于 2" in result.skipped_reason

    def test_minvar_writes_output_file(self, rerooter, temp_dir):
        out = temp_dir / "minvar.nwk"

        result = rerooter.root_tree(
            PhylogeneticTree("(A:1,(B:1,(C:1,D:1):1):1);"),
            method=RootingMethod.MINVAR,
            output_path=out,
        )

        assert out.exists()
        assert out.read_text(encoding="utf-8").strip() == result.rooted_tree.newick

    # ---------------- 只读查询 ---------------- #

    def test_is_rooted_true(self, rerooter):
        """测试已定时检测"""
        assert rerooter.is_rooted(PhylogeneticTree(ROOTED_NWK)) is True

    def test_is_rooted_false(self, rerooter):
        """测试未定时检测（三叉根 = 未定根的通行表示）"""
        assert rerooter.is_rooted(PhylogeneticTree(UNROOTED_NWK)) is False
        assert rerooter.is_rooted(PhylogeneticTree(STAR_NWK)) is False

    def test_get_root_children(self, rerooter):
        """测试获取根节点子分支"""
        result = rerooter.get_root_children(PhylogeneticTree(ROOTED_NWK))

        assert sorted(sorted(child) for child in result) == [
            ["chimp", "human"],
            ["mouse", "rat"],
        ]

    def test_get_root_children_not_binary(self, rerooter):
        """测试非二叉树获取根节点子分支"""
        assert (
            rerooter.get_root_children(PhylogeneticTree("(A:0.1,B:0.1,C:0.1);")) == []
        )

    def test_get_root_children_preserves_original_shape(self, rerooter):
        """返回形态保持 List[List[str]]（外部调用方依赖）。"""
        result = rerooter.get_root_children(PhylogeneticTree(ROOTED_NWK))

        assert isinstance(result, list) and len(result) == 2
        assert all(isinstance(child, list) for child in result)
        assert all(isinstance(name, str) for child in result for name in child)

    def test_accepts_duck_typed_tree(self, rerooter):
        """只带 ``.newick`` 的鸭子类型输入照样能定根（不给 AttributeError）。"""

        class _OnlyNewick:
            newick = ROOTED_NWK

        result = rerooter.root_tree(_OnlyNewick(), method=RootingMethod.MIDPOINT)

        assert result.was_rerooted is True
        assert sorted(result.rooted_tree.tip_names) == [
            "chimp",
            "human",
            "mouse",
            "rat",
        ]

    def test_rejects_tree_without_newick(self, rerooter):
        class _Nothing:
            pass

        with pytest.raises(TreeValidationError, match="PhylogeneticTree"):
            rerooter.root_tree(_Nothing(), method=RootingMethod.MIDPOINT)


class TestMadRooting:
    """MAD 通路仍走外部程序（不依赖 ete3），并统一经过 RootingResult 收尾。

    C-32 的关键断言：输出**按名字**交给 MAD，找不到时先做 Newick 形状检查，
    而"没有产出文件"必须当场报出来，不能推迟成毫不相干的解析错误。
    """

    #: 一个"真的换了根位置"、且叶集合与 :data:`ROOTED_NWK` 完全一致的输出
    ROOTED_OUTPUT = "((rat:0.3,mouse:0.3):0.1,(human:0.1,chimp:0.1):0.2);"
    #: shell 片段里的输出目标：我们显式给的名字 / MAD 自己的 ``{input}_rooted.nwk``
    NAMED = '"$2"'
    GUESSED = '"${1%.nwk}_rooted.nwk"'

    @staticmethod
    def _fake_mad(tmp_path, name, body):
        """写一个假 MAD：``$1`` 为输入树，``$2`` 为我们显式指定的输出文件名。"""
        wrapper = tmp_path / name
        wrapper.write_text("#!/bin/sh\n" + body + "\n", encoding="utf-8")
        wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC)
        return wrapper

    @classmethod
    def _write_tree(cls, tree_text: str, target: str) -> str:
        """shell 片段：把 ``tree_text`` 写进 ``target``（target 原样拼进命令）。"""
        return f"printf '%s\\n' {shlex.quote(tree_text)} > {target}"

    @pytest.fixture
    def fake_mad(self, tmp_path):
        """按约定：认第二个位置参数作为输出文件名。"""
        return self._fake_mad(
            tmp_path, "mad", self._write_tree(self.ROOTED_OUTPUT, self.NAMED)
        )

    def test_mad_root_via_external_program(self, tmp_path, fake_mad):
        service = TreeRootingService()
        tree = PhylogeneticTree(ROOTED_NWK)

        result = service.root_tree(
            tree,
            method=RootingMethod.MAD,
            mad_bin=str(fake_mad),
            output_path=tmp_path / "mad.nwk",
        )

        assert result.method == RootingMethod.MAD
        assert sorted(result.rooted_tree.tip_names) == sorted(tree.tip_names)
        assert (tmp_path / "mad.nwk").exists()

    def test_mad_is_told_the_output_filename_explicitly(self, tmp_path):
        """C-32：不再猜 ``{stem}_rooted.nwk``，而是把名字传给它。"""
        recorder = tmp_path / "argv.txt"
        body = (
            f"printf '%s\\n' \"$@\" > {shlex.quote(str(recorder))}\n"
            + self._write_tree(self.ROOTED_OUTPUT, self.NAMED)
        )
        fake = self._fake_mad(tmp_path, "mad", body)

        result = TreeRootingService().root_tree(
            PhylogeneticTree(ROOTED_NWK), method=RootingMethod.MAD, mad_bin=str(fake)
        )

        argv = recorder.read_text(encoding="utf-8").split()
        assert len(argv) == 2, "第二个参数必须是输出文件名"
        assert argv[0].endswith("input.nwk")
        assert Path(argv[1]).name == "rooted.nwk"
        assert sorted(result.rooted_tree.tip_names) == [
            "chimp",
            "human",
            "mouse",
            "rat",
        ]

    def test_mad_that_ignores_output_argument_is_still_found(self, tmp_path):
        """只会按自己命名规则出文件的 MAD 版本：靠私有临时目录确定性定位。"""
        fake = self._fake_mad(
            tmp_path,
            "mad",
            self._write_tree(self.ROOTED_OUTPUT, self.GUESSED),
        )

        result = TreeRootingService().root_tree(
            PhylogeneticTree(ROOTED_NWK), method=RootingMethod.MAD, mad_bin=str(fake)
        )

        assert sorted(result.rooted_tree.tip_names) == [
            "chimp",
            "human",
            "mouse",
            "rat",
        ]

    def test_mad_stdout_fallback_requires_newick_shape(self, tmp_path):
        """只往 stdout 打树的版本：形状检查通过才交给解析器。"""
        fake = self._fake_mad(tmp_path, "mad", f"printf '{self.ROOTED_OUTPUT}\\n'")

        result = TreeRootingService().root_tree(
            PhylogeneticTree(ROOTED_NWK), method=RootingMethod.MAD, mad_bin=str(fake)
        )

        assert sorted(result.rooted_tree.tip_names) == [
            "chimp",
            "human",
            "mouse",
            "rat",
        ]

    def test_mad_with_no_output_fails_loudly_not_as_parser_error(self, tmp_path):
        """真实成因必须当场报：既没出文件、stdout 又只是日志。"""
        fake = self._fake_mad(
            tmp_path,
            "mad",
            "echo 'MAD 1.0.3 reading input.nwk'; echo 'done in 0.01 s'",
        )

        with pytest.raises(TreeValidationError, match="produced no rooted tree") as exc:
            TreeRootingService().root_tree(
                PhylogeneticTree(ROOTED_NWK),
                method=RootingMethod.MAD,
                mad_bin=str(fake),
            )
        # 错误信息自己就把"三条定位尝试"讲清楚，不再把锅甩给 Newick 解析器
        assert "not a Newick tree" in str(exc.value)

    def test_mad_named_output_with_non_newick_content_is_reported_as_such(
        self, tmp_path
    ):
        """失败路径二：文件按名字写了，但内容不是树——与"没有文件"区分开。"""
        fake = self._fake_mad(
            tmp_path,
            "mad",
            f"printf '%s\\n' 'MAD log line one' 'second line' > {self.NAMED}",
        )

        with pytest.raises(TreeValidationError, match="content is not a Newick tree"):
            TreeRootingService().root_tree(
                PhylogeneticTree(ROOTED_NWK),
                method=RootingMethod.MAD,
                mad_bin=str(fake),
            )

    def test_mad_returning_a_different_tree_is_rejected(self, tmp_path):
        """失败路径三：形状像 Newick 但不是这棵树（定根不改叶集合）。"""
        fake = self._fake_mad(
            tmp_path,
            "mad",
            self._write_tree(
                "((dog:0.3,cat:0.3):0.1,(human:0.1,chimp:0.1):0.2);", self.NAMED
            ),
        )

        with pytest.raises(TreeValidationError, match="not the same tree as the input"):
            TreeRootingService().root_tree(
                PhylogeneticTree(ROOTED_NWK),
                method=RootingMethod.MAD,
                mad_bin=str(fake),
            )

    def test_mad_missing_binary_raises(self, rerooter, tree):
        """外部程序不存在时必须包装成 TreeValidationError（§六.39 的既有纪律）。"""
        with pytest.raises(TreeValidationError):
            rerooter.root_tree(
                tree, method=RootingMethod.MAD, mad_bin="definitely-not-installed-mad"
            )

    def test_mad_nonzero_exit_reports_stderr(self, tmp_path):
        fake = self._fake_mad(tmp_path, "mad", "echo 'bad newick in input' >&2; exit 3")

        with pytest.raises(TreeValidationError, match="MAD failed") as exc:
            TreeRootingService().root_tree(
                PhylogeneticTree(ROOTED_NWK),
                method=RootingMethod.MAD,
                mad_bin=str(fake),
            )
        assert "bad newick in input" in str(exc.value)


class TestNewickShapeCheck:
    """C-32 的形状闸门本身：既不能放过日志，也不能放过合法的树"""

    @pytest.mark.parametrize(
        "text",
        [
            "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);",
            "(A:1,B:1);",
            "A:1;",
            "(((A:1,B:1):1,C:1):1,D:1);",
            "((A:1,B:1)AB:0.5,(C:1,D:1)CD:0.5);",
            "((A:1,\n B:1):0.5,C:1.5);",  # 折行书写的树
            "('Homo sapiens':0.1,Pan:0.2);",  # 带引号的空格标签
            "((A:1,B:1)root:0.5,C:1.5);",
            "((A:1,B:1)[&NHX foo=bar]:0.5,C:1.5);",  # 带注释的树
            "(菌种A:1,菌种B:1);",  # 非 ASCII 标签不能被误判成"不是树"
        ],
    )
    def test_accepts_real_newick(self, text):
        assert _looks_like_newick(text) is True

    @pytest.mark.parametrize(
        "text",
        [
            "",
            "   ",
            "MAD 1.0.3 reading input.nwk",
            "done in 0.01 s;",  # 句子冒充语句
            "MAD done (0.1s);",  # 括号配平但带空格
            "((A:1,B:1);",  # 括号不配平
            "(A:1,B:1);(C:1,D:1);",  # 两条语句
            "tree 'x' = (A:1,B:1);",  # Nexus 前缀
            "Error: cannot parse the input tree;",
            "MAD finished; wrote 4 nodes",  # 日志里带分号
            "reading tree (A,B); done",  # 分号不在结尾
            "100.0",  # 只有一个数字
        ],
    )
    def test_rejects_non_newick(self, text):
        assert _looks_like_newick(text) is False


class TestNoTempFileLeak:
    """MAD 的临时文件清理（保持既有正确行为，§六.39）"""

    @staticmethod
    def _leaked():
        from tempfile import gettempdir

        return set(Path(gettempdir()).glob("*_rooted.nwk")) | set(
            Path(gettempdir()).glob("phylodater_mad_*")
        )

    def test_temp_files_removed(self, tmp_path):
        service = TreeRootingService()
        tree = PhylogeneticTree(ROOTED_NWK)
        before = self._leaked()

        with pytest.raises(TreeValidationError):
            service.root_tree(
                tree, method=RootingMethod.MAD, mad_bin="definitely-not-installed-mad"
            )

        assert self._leaked() == before

    def test_temp_dir_removed_on_success(self, tmp_path):
        """成功路径（含"按自己命名规则出文件"的版本）也不能留下中间文件。"""
        fake = TestMadRooting._fake_mad(
            tmp_path,
            "mad",
            TestMadRooting._write_tree(
                TestMadRooting.ROOTED_OUTPUT, TestMadRooting.GUESSED
            ),
        )
        before = self._leaked()

        service = TreeRootingService()
        service.root_tree(
            PhylogeneticTree(ROOTED_NWK), method=RootingMethod.MAD, mad_bin=str(fake)
        )

        assert self._leaked() == before
