"""可视化层回归测试：地质年代缩写（#17）+ 根年龄推断（B-13）+ 尖端年龄（B-14）。

关于测试替身
------------
真正的绘图路径依赖 ete3/ete4 树对象，而 ete3 3.1.x 在本工作区的解释器上根本无法
导入（Python >=3.13 移除了 stdlib ``cgi``，见审阅报告 A-5），ete4 也未安装。
因此这里用一个**最小 ete3 形状的树替身**（``.name/.dist/.up/.children/.is_leaf()/
.is_root()/.traverse()/.get_distance()/.add_feature()/.ladderize()``）驱动 ``viz``
层——可视化层本来就按鸭子类型访问节点，不需要任何 ete 依赖。替身里 root-to-tip
距离按"沿 ``up`` 链累加枝长"计算，与 ete3 语义一致。
"""

import re

import matplotlib

matplotlib.use("Agg")  # headless

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from phylodater.viz.dating_figure import DatingFigure  # noqa: E402
from phylodater.viz.geo_plot import (  # noqa: E402
    GeoPlotter,
    get_geo_abbreviations,
    load_ics_data,
)
from phylodater.viz.tree_plot import (  # noqa: E402
    _WARNED_KEYS,
    NodeCoordinate,
    as_ete_tree,
    compute_node_coordinates,
    distance_from_root,
    get_recorded_node_age,
    has_plot_tree_shape,
    infer_root_age_from_branch_lengths,
    resolve_root_age,
)

# --------------------------------------------------------------------------- #
# 最小 ete3 形状树替身
# --------------------------------------------------------------------------- #


class EteLikeNode:
    """ete3 TreeNode 的最小等价物（无枝长的 Newick 记为 0.0，与 ete3 一致）。"""

    def __init__(self, name="", dist=0.0):
        self.name = name
        self.dist = dist
        self.up = None
        self.children = []
        self.features = set()

    # --- ete3 接口 ---
    def is_leaf(self):
        return not self.children

    def is_root(self):
        return self.up is None

    def add_feature(self, key, value):
        setattr(self, key, value)
        self.features.add(key)

    def traverse(self, strategy="postorder"):
        for child in self.children:
            yield from child.traverse(strategy=strategy)
        yield self

    def get_distance(self, target):
        """根到后代节点的枝长累加（替身里 self 始终用作根）。"""
        total = 0.0
        current = target
        while current is not None:
            if current is self:
                return total
            total += float(current.dist)
            current = current.up
        raise ValueError(f"节点 {getattr(target, 'name', target)!r} 不在该树中")

    def ladderize(self, reverse=False):
        def size(node):
            return 1 if node.is_leaf() else sum(size(c) for c in node.children)

        self.children.sort(key=size, reverse=reverse)
        for child in self.children:
            child.ladderize(reverse=reverse)

    # --- 断言辅助 ---
    def by_name(self, name):
        for node in self.traverse():
            if node.name == name:
                return node
        raise KeyError(name)


class NoGetDistanceNode(EteLikeNode):
    """没有 ``get_distance()`` 的树后端：只能走 ``up`` 链累加。"""

    get_distance = None


_LABEL_RE = re.compile(r"[^,()();:]*")
_LENGTH_RE = re.compile(r":([-+0-9.eE]+)")


def parse_newick(newick, missing_branch_length=0.0, node_class=EteLikeNode):
    """极简 Newick 解析器（仅供测试；标签内不含 ``:``，如 ``[&&NHX:…]`` 写法）。"""
    text = newick.strip().rstrip(";")
    pos = 0

    def parse_subtree():
        nonlocal pos
        node = node_class()
        if pos < len(text) and text[pos] == "(":
            pos += 1
            while True:
                child = parse_subtree()
                child.up = node
                node.children.append(child)
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
                raise AssertionError(f"Newick 解析失败，位置 {pos}: {text!r}")
        label = _LABEL_RE.match(text, pos)
        node.name = label.group(0)
        pos = label.end()
        length = _LENGTH_RE.match(text, pos)
        if length:
            node.dist = float(length.group(1))
            pos = length.end()
        else:
            node.dist = missing_branch_length
        return node

    root = parse_subtree()
    assert pos == len(text), f"Newick 末尾有未消费内容: {text[pos:]!r}"
    return root


#: 严格超度量、根高 100 Ma 的四尖端时间树。
ULTRAMETRIC_100 = "((A:40,B:40)AB:60,(C:70,D:70)CD:30)root:0;"
#: 同一拓扑，但 (C,D) 这对尖端比 (A,B) 新 60 Ma（含灭绝类群 / 尖端定年的形态）。
NON_ULTRAMETRIC = "((A:40,B:40)AB:60,(C:10,D:10)CD:30)root:0;"


class SpyLogger:
    """记录 PhyloDaterLogger 调用的最小替身。"""

    def __init__(self):
        self.messages = {"debug": [], "info": [], "warning": [], "error": []}

    def __getattr__(self, level):
        def record(message, **_kwargs):
            self.messages[level].append(str(message))

        return record

    def any(self, needle, level="warning"):
        return any(needle in message for message in self.messages[level])


@pytest.fixture(autouse=True)
def _clean_state():
    """每个用例：清空逐节点告警去重表 + 关闭 matplotlib figure。"""
    _WARNED_KEYS.clear()
    yield
    _WARNED_KEYS.clear()
    plt.close("all")


@pytest.fixture
def spy_tree_plot(monkeypatch):
    """截获 ``viz.tree_plot`` 的日志输出。"""
    import phylodater.viz.tree_plot as tree_plot_module

    spy = SpyLogger()
    monkeypatch.setattr(tree_plot_module, "logger", spy)
    return spy


# --------------------------------------------------------------------------- #
# #17 地质年代命名碰撞（原有回归测试，保持不变）
# --------------------------------------------------------------------------- #


class TestGeoAbbreviations:
    """确保地质年代缩写无命名碰撞"""

    def test_no_duplicate_abbreviations(self):
        """任意两个地质年代不得映射到相同缩写（如 Miaolingian/Miocene 曾碰撞）"""
        abbr = get_geo_abbreviations()
        values = list(abbr.values())
        assert len(values) == len(set(values)), "存在重复缩写: " + str(
            [v for v in set(values) if values.count(v) > 1]
        )

    def test_miaolingian_distinct_from_miocene(self):
        """Miaolingian 与 Miocene 必须使用不同缩写"""
        abbr = get_geo_abbreviations()
        assert abbr["Miaolingian"] != abbr["Miocene"]
        assert abbr["Miaolingian"] == "Mia"
        assert abbr["Miocene"] == "Mi"


# --------------------------------------------------------------------------- #
# B-13：root_age 缺省时必须从分支长度真实推断
# --------------------------------------------------------------------------- #


class TestRootAgeInference:
    def test_documented_example_without_root_age(self, tmp_path):
        """模块 docstring 里那段不传 root_age 的用法示例不得再塌成 x=0"""
        tree = parse_newick(ULTRAMETRIC_100)

        figure = DatingFigure(width=6, height=4, dpi=60)
        figure.add_tree(tree).add_geo_scale().render()

        assert figure._root_age == pytest.approx(100.0)

        coords = figure.get_node_coordinates()
        ages = {coord.name: coord.age for coord in coords.values()}
        assert ages["root"] == pytest.approx(100.0)
        assert ages["AB"] == pytest.approx(40.0)  # 100 - 60
        assert ages["CD"] == pytest.approx(70.0)  # 100 - 30
        for tip in ("A", "B", "C", "D"):
            assert ages[tip] == pytest.approx(0.0)

        # 时间轴必须张得开，且方向仍是"老在左、现在在右"（§六.17 不得回归）
        tree_xlim = figure.ax_tree.get_xlim()
        geo_xlim = figure.ax_geo.get_xlim()
        assert tree_xlim[0] > tree_xlim[1]
        assert tree_xlim[0] > 100.0 > tree_xlim[1] >= 0.0
        assert geo_xlim == pytest.approx(tree_xlim)

        output = tmp_path / "dating_figure.png"
        assert figure.save(output).exists()
        assert output.stat().st_size > 0
        figure.close()

    def test_explicit_root_age_still_overrides(self):
        """显式 root_age 优先级最高（既有契约不得破坏）"""
        tree = parse_newick(ULTRAMETRIC_100)
        figure = DatingFigure()
        figure.add_tree(tree, root_age=500.0)
        assert figure._root_age == pytest.approx(500.0)
        assert resolve_root_age(tree, 500.0) == pytest.approx(500.0)

    def test_root_age_conflicting_with_branch_lengths_warns(self, monkeypatch):
        """枝长与 root_age 不一致时尖端不再落在 0：必须说出来（B-14 的副作用）"""
        import phylodater.viz.tree_plot as tree_plot_module

        spy = SpyLogger()
        monkeypatch.setattr(tree_plot_module, "logger", spy)

        coords = compute_node_coordinates(parse_newick(ULTRAMETRIC_100), root_age=500.0)
        by_name = {coord.name: coord for coord in coords.values()}
        assert by_name["A"].age == pytest.approx(400.0)  # 500 - (60 + 40)
        assert spy.any("与树自身分支长度隐含的")

    def test_inference_is_max_root_to_tip_path(self):
        tree = parse_newick(NON_ULTRAMETRIC)
        assert infer_root_age_from_branch_lengths(tree) == pytest.approx(100.0)

    def test_non_ultrametric_inference_warns_about_time_tree(self, spy_tree_plot):
        """非超度量输入：按最长 root-to-tip 取值并告警（B-13 + B-18）"""
        root_age = resolve_root_age(parse_newick(NON_ULTRAMETRIC))
        assert root_age == pytest.approx(100.0)
        assert spy_tree_plot.any("不是时间树")

    def test_cladogram_without_branch_lengths_raises(self):
        """纯拓扑树：拒绝出图，而不是把根年龄当成 0（B-13 的核心症状）"""
        tree = parse_newick("((A,B)AB,(C,D)CD)root;")
        with pytest.raises(ValueError, match="无法确定根年龄"):
            DatingFigure().add_tree(tree)
        with pytest.raises(ValueError, match="无法确定根年龄"):
            compute_node_coordinates(tree)

    def test_resolve_root_age_reads_root_annotation_first(self):
        tree = parse_newick(ULTRAMETRIC_100)
        tree.add_feature("age", 250.0)
        assert resolve_root_age(tree) == pytest.approx(250.0)

    def test_resolve_root_age_falls_back_to_recorded_node_ages(self, spy_tree_plot):
        """无枝长也无根注解时，用树上已记录节点年龄的最大值定标并告警"""
        tree = parse_newick("((A,B)AB,(C,D)CD)root;", missing_branch_length=None)
        tree.by_name("AB").add_feature("age", 80.0)
        assert resolve_root_age(tree) == pytest.approx(80.0)
        assert spy_tree_plot.any("按树上已有节点年龄的最大值推断根年龄")

    def test_negative_and_nan_root_age_rejected(self):
        tree = parse_newick(ULTRAMETRIC_100)
        with pytest.raises(ValueError, match="为负"):
            resolve_root_age(tree, -5.0)
        with pytest.raises(ValueError, match="有限"):
            resolve_root_age(tree, float("nan"))

    def test_distance_falls_back_to_up_chain_without_get_distance(self):
        """后端没有 get_distance() 时按 up 链累加，结果必须一致"""
        tree = parse_newick(ULTRAMETRIC_100, node_class=NoGetDistanceNode)
        assert distance_from_root(tree, tree.by_name("A")) == pytest.approx(100.0)
        assert resolve_root_age(tree) == pytest.approx(100.0)

    def test_recorded_root_age_from_result_wins_after_add_tree(self):
        """add_tree 之后的 add_result/注解仍优先于推断出的根年龄（原契约）"""
        tree = parse_newick(ULTRAMETRIC_100)
        figure = DatingFigure()
        figure.add_tree(tree)  # 推断 100 并写回根节点
        tree.add_feature("age", 150.0)  # 相当于 add_result 把根年龄改成 150
        coords = compute_node_coordinates(tree, root_age=figure._root_age)
        assert coords[tree].age == pytest.approx(150.0)

    def test_documented_example_with_a_real_ete_backend(self):
        """同一用例跑在真实 ete 后端上（A-5 环境下自动跳过）。

        替身与真实 ete 对象的鸭子类型可能漂移，因此只要后端可导入就再跑一遍。
        """
        tree_class = None
        for module_name in ("ete4", "ete3"):
            try:
                tree_class = __import__(module_name, fromlist=["Tree"]).Tree
                break
            except ImportError:
                continue
        if tree_class is None:
            pytest.skip("ete4/ete3 均不可导入（见审阅报告 A-5）")

        tree = tree_class(ULTRAMETRIC_100, format=1)
        figure = DatingFigure(width=6, height=4, dpi=60)
        figure.add_tree(tree).add_geo_scale().render()
        ages = {c.name: c.age for c in figure.get_node_coordinates().values()}
        assert ages["root"] == pytest.approx(100.0)
        assert ages["AB"] == pytest.approx(40.0)
        assert ages["CD"] == pytest.approx(70.0)
        assert ages["A"] == pytest.approx(0.0)
        assert figure.ax_tree.get_xlim()[0] > figure.ax_tree.get_xlim()[1]
        figure.close()

    def test_phylodater_tree_object_uses_public_backend(self):
        """PhylogeneticTree 走模型层入口取树后端，不在 viz 里 import ete3"""
        PhylogeneticTree = pytest.importorskip("phylodater.models").PhylogeneticTree
        tree = PhylogeneticTree(ULTRAMETRIC_100)
        figure = DatingFigure(width=6, height=4, dpi=60)
        try:
            figure.add_tree(tree)
        except ImportError as exc:  # A-5：ete3 在 3.13+ 无法导入
            pytest.skip(f"本机没有可用的 ete 树后端: {exc}")
        assert figure._root_age == pytest.approx(100.0)


# --------------------------------------------------------------------------- #
# B-14：叶/尖端年龄不再硬编码为 0
# --------------------------------------------------------------------------- #


class TestTipAges:
    def test_ultrametric_time_tree_tips_still_land_on_zero(self):
        """回归保护：超度量时间树上尖端自然回到 0，与修复前一致"""
        coords = compute_node_coordinates(parse_newick(ULTRAMETRIC_100))
        tips = [c for c in coords.values() if c.is_leaf]
        assert tips
        assert all(c.age == pytest.approx(0.0) for c in tips)

    def test_non_ultrametric_tips_use_same_formula_as_internal_nodes(self):
        """含灭绝类群 / 尖端定年：尖端落在枝长隐含的时间点上，而不是"现在"。"""
        coords = compute_node_coordinates(parse_newick(NON_ULTRAMETRIC))
        by_name = {coord.name: coord for coord in coords.values()}
        assert by_name["A"].age == pytest.approx(0.0)  # 100 - (60 + 40)
        assert by_name["C"].age == pytest.approx(60.0)  # 100 - (30 + 10)
        assert by_name["CD"].age == pytest.approx(70.0)  # 100 - 30

    def test_extinct_tip_uses_recorded_age_feature(self):
        """树里记录了化石尖端年龄时必须用它（报告 B-14 的"优先使用该特征"）"""
        tree = parse_newick(ULTRAMETRIC_100)
        tree.by_name("A").add_feature("age", 66.0)  # 地层年龄 66 Ma
        coords = compute_node_coordinates(tree)
        by_name = {coord.name: coord for coord in coords.values()}
        assert by_name["A"].age == pytest.approx(66.0)
        assert by_name["B"].age == pytest.approx(0.0)

    def test_tip_label_date_suffix_is_used(self):
        """``A|-0.045`` —— 沿用本项目的采样日期轴约定（负值 = 过去，单位 Ma）"""
        tree = parse_newick("((A|-0.045:40,B:40)AB:60,(C:70,D:70)CD:30)root:0;")
        coords = compute_node_coordinates(tree)
        by_name = {coord.name: coord for coord in coords.values()}
        assert by_name["A|-0.045"].age == pytest.approx(0.045)

    def test_nhx_tip_date_annotation_is_used(self):
        tree = parse_newick("((A[&&DATE=-2.5]:40,B:40)AB:60,(C:70,D:70)CD:30)root:0;")
        tip = tree.by_name("A[&&DATE=-2.5]")
        assert get_recorded_node_age(tip, node_kind="尖端") == pytest.approx(2.5)

    def test_recorded_tip_age_is_reported_in_log(self, spy_tree_plot):
        """报告要求"报告有多少叶被赋值" """
        tree = parse_newick(ULTRAMETRIC_100)
        tree.by_name("A").add_feature("age", 100.0)
        compute_node_coordinates(tree)
        assert spy_tree_plot.any("1/4 个尖端按树内记录的采样日期", level="info")

    def test_calendar_tip_date_is_not_guessed(self, spy_tree_plot):
        """日历式日期要先假设"现在是哪一年"：不猜，只告警并回退到枝长口径"""
        tree = parse_newick("((A|2019-12-31:40,B:40)AB:60,(C:70,D:70)CD:30)root:0;")
        coords = compute_node_coordinates(tree)
        by_name = {coord.name: coord for coord in coords.values()}
        assert by_name["A|2019-12-31"].age == pytest.approx(0.0)
        assert spy_tree_plot.any("无法解释为距今 Ma")

    def test_bare_four_digit_tip_date_is_not_guessed(self, spy_tree_plot):
        tree = parse_newick("((A|2019:40,B:40)AB:60,(C:70,D:70)CD:30)root:0;")
        compute_node_coordinates(tree)
        assert spy_tree_plot.any("裸四位数")

    def test_recorded_negative_age_is_not_folded_to_positive(self, spy_tree_plot):
        """负的 age 注解是方向/单位写错的信号：不折成正数（同 C-55 的取向）"""
        tree = parse_newick(ULTRAMETRIC_100)
        tree.by_name("A").add_feature("age", -66.0)
        coords = compute_node_coordinates(tree)
        by_name = {coord.name: coord for coord in coords.values()}
        assert by_name["A"].age == pytest.approx(0.0)  # 回退到枝长公式
        assert spy_tree_plot.any("年龄轴上不可能")

    def test_missing_branch_length_falls_back_to_zero_with_warning(self, spy_tree_plot):
        """既无年龄注解、整条路径也无枝长：按 0 绘制但必须告警（不再静默）"""
        tree = parse_newick(
            "((A:40,B:40)AB:60,(C,D)CD:30)root:0;", missing_branch_length=None
        )
        tree.by_name("CD").dist = None
        coords = compute_node_coordinates(tree)
        by_name = {coord.name: coord for coord in coords.values()}
        assert by_name["A"].age == pytest.approx(0.0)
        assert by_name["C"].age == 0.0
        assert spy_tree_plot.any("只能按'现在'(0 Ma) 绘制")


# --------------------------------------------------------------------------- #
# C-13 / C-14：算不出的距离与退化时间轴不再伪装成数据
# --------------------------------------------------------------------------- #


class TestFailureModesAreVisible:
    def test_unreachable_node_raises_instead_of_plotting_zero(self):
        root = parse_newick(ULTRAMETRIC_100)
        orphan = EteLikeNode(name="孤儿", dist=1.0)  # up 链不指向任何根
        with pytest.raises(ValueError, match="无法计算节点 '孤儿'"):
            distance_from_root(root, orphan)

    def test_broken_parent_chain_aborts_coordinate_computation(self):
        root = parse_newick(ULTRAMETRIC_100)
        broken = EteLikeNode(name="断链", dist=5.0)
        root.children[0].children.append(broken)  # broken.up 仍是 None
        with pytest.raises(ValueError, match="无法计算节点 '断链'"):
            compute_node_coordinates(root, root_age=100.0)

    def test_nan_branch_length_is_not_silently_zero(self):
        root = parse_newick(ULTRAMETRIC_100)
        root.by_name("A").dist = float("nan")
        with pytest.raises(ValueError, match="非有限值"):
            compute_node_coordinates(root, root_age=100.0)

    def test_age_range_without_coordinates_raises(self):
        with pytest.raises(ValueError, match="尚未计算出任何节点坐标"):
            DatingFigure()._get_age_range()

    def test_degenerate_age_range_raises(self):
        """C-14：零宽度区间不得被兜底值伪装成"有数据"。"""
        figure = DatingFigure()
        figure._node_coords = {
            "n1": NodeCoordinate(0.0, 1.0, "A", True, False, 0.0, {}),
            "n2": NodeCoordinate(0.0, 2.0, "B", True, False, 0.0, {}),
        }
        with pytest.raises(ValueError, match="时间轴退化"):
            figure._get_age_range()

    def test_non_finite_ages_rejected_by_age_range(self):
        figure = DatingFigure()
        figure._node_coords = {
            "n1": NodeCoordinate(0.0, 1.0, "A", True, False, float("nan"), {}),
        }
        with pytest.raises(ValueError, match="没有任何有限年龄值"):
            figure._get_age_range()

    def test_render_rejects_degenerate_axis_even_without_geo_scale(self):
        """根年龄被写坏时（枝长全 0 的相对时间树），出图前必须中止"""
        tree = parse_newick("((A,B)AB,(C,D)CD)root;")
        for node in tree.traverse():
            node.add_feature("age", 1.0)
        figure = DatingFigure()
        figure.add_tree(tree)
        with pytest.raises(ValueError, match="时间轴退化"):
            figure.render()


# --------------------------------------------------------------------------- #
# 地质时标：DatingFigure 传的是降序 (老界, 新界)，两个朝向都得画对
# --------------------------------------------------------------------------- #


class TestGeoTimescaleRange:
    @staticmethod
    def _render(age_range, **kwargs):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        try:
            GeoPlotter(split_rows=False, **kwargs).plot(
                ax=ax, age_range=age_range, geo_data=load_ics_data()
            )
            spans = sorted(
                (round(p.get_x(), 6), round(p.get_x() + p.get_width(), 6))
                for p in ax.patches
            )
            return ax, spans
        finally:
            plt.close(fig)

    def test_ascending_and_descending_ranges_render_identically(self):
        _, ascending = self._render((0.0, 120.0))
        axis, descending = self._render((120.0, 0.0))

        assert ascending, "时标一个色块都没画出来"
        assert ascending == descending
        assert axis.get_xlim()[0] > axis.get_xlim()[1]  # 老在左

    def test_intervals_are_clamped_into_viewport(self):
        axis, spans = self._render((100.0, 0.0))
        assert spans
        for left, right in spans:
            assert left >= -1e-9
            assert right <= 100.0 + 1e-9
            assert right > left
        assert axis.get_xlim()[0] > axis.get_xlim()[1]

    def test_age_ticks_cover_viewport(self):
        axis, _ = self._render((100.0, 0.0), age_tick_interval=20.0)
        ticks = [float(t) for t in axis.get_xticks()]
        in_view = [t for t in ticks if 0.0 <= t <= 100.0]
        assert len(in_view) >= 5

    def test_split_rows_also_renders_both_bands(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        try:
            GeoPlotter(split_rows=True).plot(
                ax=ax, age_range=(100.0, 0.0), geo_data=load_ics_data()
            )
            heights = {round(p.get_y(), 6) for p in ax.patches}
            assert ax.patches
            assert len(heights) == 2, "split_rows 应该画出上下两条带"
        finally:
            plt.close(fig)

    def test_get_intervals_for_age_finds_containing_intervals(self):
        intervals = GeoPlotter().get_intervals_for_age(100.0)
        names = {row["name"] for row in intervals}
        assert "Phanerozoic" in names
        assert "Mesozoic" in names
        assert "Cenozoic" not in names
        assert "Proterozoic" not in names


# --------------------------------------------------------------------------- #
# B-14 补：Newick 标签里的年龄/日期注解（ETE NHX 块与 phylocom 风格）
# --------------------------------------------------------------------------- #


class TestNewickAgeAnnotations:
    @pytest.mark.parametrize(
        "label,expected",
        [
            ("A{[&age=55]}", 55.0),  # 花括号包注解（某些上游写法）
            ("A[&age=55]", 55.0),
            ("A[&&NHX:age=55]", 55.0),
            ("A[&&NHX:S=homo:age=55:t=1]", 55.0),
            ("A[&&NHX:tip_age=55]", 55.0),
            ("A[&date=-0.045]", 0.045),  # 采样日期轴：负值 = 过去
            ("A[&&NHX:date=-2.5]", 2.5),
            # 与年龄无关 / 不可解释的注解：一律回退到分支长度口径
            ("A[&color=blue]", None),
            ("A[&&NHX:S=homo]", None),
            ("A[&age=notanumber]", None),
            ("A[&age=-55]", None),
            ("A", None),
        ],
    )
    def test_annotation_parsing(self, label, expected):
        node = EteLikeNode(name=label, dist=40.0)
        age = get_recorded_node_age(node, node_kind="尖端")
        if expected is None:
            assert age is None
        else:
            assert age == pytest.approx(expected)

    def test_annotated_tip_is_drawn_at_its_recorded_age(self):
        """端到端：标签里的 ``[&age=…]`` 直接决定尖端画在哪儿（B-14）。"""
        tree = parse_newick("((A[&age=55]:40,B:40)AB:60,(C:70,D:70)CD:30)root:0;")
        ages = {c.name: c.age for c in compute_node_coordinates(tree).values()}
        assert ages["A[&age=55]"] == pytest.approx(55.0)
        assert ages["B"] == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# B-14 §3："枝干在时间轴上倒挂" / 两套时间尺度混用必须可见
# --------------------------------------------------------------------------- #


class TestAgeAxisConsistencyWarnings:
    def test_ultrametric_time_tree_produces_no_warning(self, spy_tree_plot):
        compute_node_coordinates(parse_newick(ULTRAMETRIC_100))
        assert not spy_tree_plot.messages["warning"]

    def test_child_older_than_parent_warns(self, spy_tree_plot):
        tree = parse_newick(ULTRAMETRIC_100)
        tree.by_name("CD").add_feature("age", 20.0)  # 父节点 20 Ma
        tree.by_name("C").add_feature("age", 55.0)  # 子孙 55 Ma
        compute_node_coordinates(tree)
        assert spy_tree_plot.any("时间轴倒挂")
        assert spy_tree_plot.any("比其父节点")

    def test_recorded_root_age_different_from_scale_warns(self, spy_tree_plot):
        """add_result 只改写部分节点时，图上是两套时间尺度——必须说出来。"""
        tree = parse_newick(ULTRAMETRIC_100)
        coords = compute_node_coordinates(tree)  # 推断 100，并写回根节点 age
        tree.by_name("root").add_feature("age", 300.0)
        compute_node_coordinates(tree, root_age=100.0)
        assert coords[tree].age == pytest.approx(100.0)
        assert spy_tree_plot.any("图上会同时出现两套时间尺度")


# --------------------------------------------------------------------------- #
# A-5/B-10：可视化层取树后端的形状契约（不自行 import ete）
# --------------------------------------------------------------------------- #


class DendropyLikeTree:
    """没有 ``traverse()``/``root`` 的树对象（例如 DendroPy 的 Tree）。"""

    def __init__(self):
        self.node_nodes = []


class ModelLikeWithWrongBackend:
    def to_ete_tree(self):
        return DendropyLikeTree()


class ModelLikeWithBrokenBackend:
    def to_ete_tree(self):
        raise ImportError("No module named 'cgi'")


class ModelLikeWithNoBackend:
    pass


class TestTreeBackendPlumbing:
    def test_ete_like_tree_is_passed_through(self):
        tree = parse_newick(ULTRAMETRIC_100)
        assert as_ete_tree(tree) is tree

    def test_phylodater_shape_check_accepts_ete_like_only(self):
        assert has_plot_tree_shape(parse_newick(ULTRAMETRIC_100))
        assert not has_plot_tree_shape(DendropyLikeTree())
        assert not has_plot_tree_shape(object())

    def test_backend_without_plot_shape_fails_loudly(self):
        """模型层若返回非 ete 形状的树，要立刻报错而不是画出半张图。"""
        with pytest.raises(TypeError, match="不具备绘图所需的树接口"):
            as_ete_tree(ModelLikeWithWrongBackend())

    def test_unimportable_backend_reports_the_model_layer(self):
        with pytest.raises(ImportError, match="可视化用的树后端") as exc:
            as_ete_tree(ModelLikeWithBrokenBackend())
        assert "可视化层不自行导入 ete" in str(exc.value)

    def test_object_without_any_backend_raises_type_error(self):
        with pytest.raises(TypeError, match="无法为绘图取得"):
            as_ete_tree(ModelLikeWithNoBackend())
