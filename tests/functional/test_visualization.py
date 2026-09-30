"""
模块十一：可视化模块功能测试

验证树绘图、地质年代图、报告产物等，并覆盖审阅项 B-13（``root_age`` 缺省时
必须按分支长度推断，而不是当成 0）与 B-14（叶/尖端年龄不再硬编码为 0）的
**端到端**行为：断言取自真正画到画布上的图元（scatter 偏移、axes 限值、存盘
文件内容），而不只是内部字典。

关于测试替身：本工作区的解释器上 ete3 根本无法导入（3.1.x 依赖 Python 3.13 起
被移除的 stdlib ``cgi``，见审阅报告 A-5），ete4 未安装，因此这里复用
``tests/unit/test_viz_geo.py`` 里的"最小 ete3 形状"树替身来驱动真实绘图路径
——可视化层按鸭子类型访问节点，替身足以跑完整坐标计算与出图。
涉及 ``PhylogeneticTree`` 的用例仍走模型层真实入口，只在树后端缺失时 skip。
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402

from phylodater import PhylogeneticTree, TreePlotter  # noqa: E402
from phylodater.models import DatingResult, NodeAgeEstimate  # noqa: E402
from phylodater.models.results import CIType  # noqa: E402
from phylodater.viz.dating_figure import DatingFigure  # noqa: E402
from phylodater.viz.tree_plot import _WARNED_KEYS  # noqa: E402

# 复用单元测试里的最小 ete3 形状树替身（理由见模块 docstring 的 A-5 说明）
from tests.unit.test_viz_geo import (  # noqa: E402
    NON_ULTRAMETRIC,
    ULTRAMETRIC_100,
    SpyLogger,
    parse_newick,
)

#: 严格超度量树，但尖端 C 是 55 Ma 的化石类群（Newick 注解形式）
FOSSIL_ANNOTATED_NEWICK = "((A:40,B:40)AB:60,(C{[&age=55]}:45,D:45)CD:25)root:0;"


@pytest.fixture(autouse=True)
def _clean_state():
    """每个用例：清空逐节点告警去重表 + 关闭 matplotlib figure。"""
    _WARNED_KEYS.clear()
    yield
    _WARNED_KEYS.clear()
    plt.close("all")


def _ages_by_name(coords):
    return {coord.name: coord.age for coord in coords.values()}


def _plotted_x(ax):
    """画布上真正被画出来的节点标记横坐标（= 年龄 Ma）。"""
    points = []
    for collection in ax.findobj(PathCollection):
        points.extend(float(x) for x, _y in collection.get_offsets())
    return points


class TestVisualization:
    """可视化测试"""

    def test_viz_001_tree_plotter_import(self):
        """VIZ-001: 可视化模块导入"""
        from phylodater import GeoPlotter, TreePlotter

        assert TreePlotter is not None
        assert GeoPlotter is not None

    def test_viz_002_dating_tree_plot(self, tmp_path):
        """VIZ-002: 定年树绘图（真实模型对象 + 真实出图产物）"""
        # root-to-tip 分别为 0.4 / 0.5 / 0.4 —— 最长路径 0.5（非超度量）
        tree = PhylogeneticTree.from_newick("((A:0.1,B:0.2):0.3,C:0.4);")
        plotter = TreePlotter()
        output = tmp_path / "tree.png"

        try:
            coords = plotter.plot(tree, output=output)
        except ImportError as exc:
            # 只允许"模型层没有可用树后端"这一种跳过理由（A-5，属 models/tree.py
            # 的范围）。可视化层自身的任何异常都必须让测试失败，不能被
            # ``except Exception`` 吞掉——那正是审阅项 B-9 批评的"永不可能失败的空调用"。
            pytest.skip(f"模型层树后端不可用（审阅项 A-5，非 viz 层问题）: {exc}")

        assert output.exists()
        assert output.stat().st_size > 0

        # B-13：未提供 root_age 时按最长 root-to-tip 路径推断，而不是当成 0
        ages = _ages_by_name(coords)
        assert max(ages.values()) == pytest.approx(0.5)
        # B-14：尖端按自身枝长落到隐含采样点，而不是全部钉在 0
        assert ages["B"] == pytest.approx(0.0)
        assert ages["A"] == pytest.approx(0.1)
        assert ages["C"] == pytest.approx(0.1)

    def test_viz_015_backend_switch(self):
        """VIZ-015: 后端切换检查"""
        backend = matplotlib.get_backend()
        assert backend is not None

    def test_viz_016_empty_result_visualization(self, monkeypatch):
        """VIZ-016: 空/错位的定年结果不得被画成一张"看起来正常"的图。

        修复前本用例是 ``reporter = type(...); assert True``——一句永不可能失败的
        断言（同审阅项 B-9 的形态），现在改为真的调用 ``add_result``。
        """
        import phylodater.viz.dating_figure as dating_figure_module

        spy = SpyLogger()
        monkeypatch.setattr(dating_figure_module, "logger", spy)

        figure = DatingFigure(width=6, height=4, dpi=60)
        figure.add_tree(parse_newick(ULTRAMETRIC_100))

        empty = DatingResult(
            method_name="none",
            run_id="empty_run",
            dated_tree_newick=ULTRAMETRIC_100,
            node_ages={},
        )
        with pytest.raises(ValueError, match="no node ages"):
            figure.add_result(empty)

        unmatched = DatingResult(
            method_name="wrong",
            run_id="unmatched_run",
            dated_tree_newick=ULTRAMETRIC_100,
            node_ages={
                "node_does_not_exist": NodeAgeEstimate(
                    mean_age=99.0, median_age=99.0, ci_type=CIType.RANGE
                )
            },
        )
        figure.add_result(unmatched)  # 不得崩溃
        assert spy.any("not found in tree"), "结果节点与树对不上时必须逐条披露"
        # 一个点都没落上 => 坐标仍是推断出的 100 Ma，没有被 99 静默污染
        ages = _ages_by_name(figure.render().get_node_coordinates())
        assert ages["root"] == pytest.approx(100.0)
        figure.close()


class TestRootAgeInferenceEndToEnd:
    """B-13：``add_tree()`` 不带 root_age 的整条链路必须出真图。"""

    def test_module_docstring_example_renders_a_real_pdf(self, tmp_path):
        """照模块 docstring 第一段代码敲（不传 root_age）：不得塌成 x=0。"""
        tree = parse_newick(ULTRAMETRIC_100)
        figure = (
            DatingFigure(width=8, height=6, dpi=72)
            .add_tree(tree)
            .add_geo_scale()
            .render()
        )

        assert figure._root_age == pytest.approx(100.0)
        ages = _ages_by_name(figure.get_node_coordinates())
        assert ages["root"] == pytest.approx(100.0)
        assert ages["AB"] == pytest.approx(40.0)
        assert ages["CD"] == pytest.approx(70.0)
        assert all(ages[t] == pytest.approx(0.0) for t in ("A", "B", "C", "D"))

        # 地质时标真的画到了画布上（C-14 修复前会兜底成一条没有数据的 0-100 Ma 轴）
        assert len(figure.ax_geo.patches) > 0
        xlim = figure.ax_tree.get_xlim()
        assert xlim[0] > xlim[1]  # 老在左、现在在右（§六.17，不得回归）
        assert xlim[0] > 100.0 and xlim[1] >= 0.0
        assert figure.ax_geo.get_xlim() == pytest.approx(xlim)

        output = tmp_path / "figure.pdf"
        assert figure.save(output).exists()
        assert output.stat().st_size > 1024
        figure.close()

    def test_explicit_root_age_still_wins(self):
        tree = parse_newick(ULTRAMETRIC_100)
        figure = DatingFigure(width=6, height=4, dpi=60)
        figure.add_tree(tree, root_age=500.0).render()
        assert figure._root_age == pytest.approx(500.0)
        ages = _ages_by_name(figure.get_node_coordinates())
        assert ages["root"] == pytest.approx(500.0)
        figure.close()

    def test_cladogram_is_refused_instead_of_collapsed(self):
        """没有任何枝长/年龄信息时拒绝出图（修复前：所有节点画在原点且无警告）。"""
        tree = parse_newick("((A,B)AB,(C,D)CD)root;")
        with pytest.raises(ValueError, match="无法确定根年龄"):
            DatingFigure().add_tree(tree)

    def test_degenerate_axis_never_reaches_the_file(self, tmp_path):
        """全树同一年龄（零宽度时标）时不得留下任何图形产物。"""
        tree = parse_newick(ULTRAMETRIC_100)
        for node in tree.traverse():
            node.add_feature("age", 42.0)
        figure = DatingFigure(width=6, height=4, dpi=60)
        figure.add_tree(tree)
        output = tmp_path / "never.pdf"
        with pytest.raises(ValueError, match="时间轴退化"):
            figure.render()
        assert not output.exists()
        figure.close()

    def test_dating_result_ages_take_precedence(self, monkeypatch):
        """``add_result`` 写进树的年龄优先于枝长推断（既有契约）。

        只覆盖了部分节点时，其余节点仍按旧尺度定位 —— 这种"两套时间尺度"必须
        告警（B-14 §3 的"枝干在时间轴上倒挂"观感就是这么来的）。
        """
        import phylodater.viz.tree_plot as tree_plot_module

        spy = SpyLogger()
        monkeypatch.setattr(tree_plot_module, "logger", spy)

        tree = parse_newick(ULTRAMETRIC_100)
        figure = DatingFigure(width=6, height=4, dpi=60)
        figure.add_tree(tree)
        result = DatingResult(
            method_name="mcmctree",
            run_id="run_1",
            dated_tree_newick=ULTRAMETRIC_100,
            node_ages={
                "root": NodeAgeEstimate(
                    mean_age=300.0,
                    median_age=300.0,
                    ci_lower=280.0,
                    ci_upper=320.0,
                    ci_type=CIType.HPD95,
                ),
                "CD": NodeAgeEstimate(
                    mean_age=210.0, median_age=210.0, ci_type=CIType.NONE
                ),
            },
        )
        figure.add_result(result).render()
        ages = _ages_by_name(figure.get_node_coordinates())
        assert ages["root"] == pytest.approx(300.0)
        assert ages["CD"] == pytest.approx(210.0)
        # AB 没有被结果覆盖 => 仍按 add_tree 推断出的 root_age=100 定位
        assert ages["AB"] == pytest.approx(40.0)
        assert spy.any("图上会同时出现两套时间尺度")
        figure.close()

    def test_age_inversion_is_reported(self, monkeypatch):
        """尖端比其父节点还老（注解写错单位/方向）时必须告警，而不是静默出图。"""
        import phylodater.viz.tree_plot as tree_plot_module

        spy = SpyLogger()
        monkeypatch.setattr(tree_plot_module, "logger", spy)

        tree = parse_newick(ULTRAMETRIC_100)
        tree.by_name("CD").add_feature("age", 20.0)  # 父节点只有 20 Ma
        tree.by_name("C").add_feature("age", 55.0)  # 子孙却 55 Ma
        coords = TreePlotter().plot(tree)
        ages = _ages_by_name(coords)
        assert ages["CD"] == pytest.approx(20.0)
        assert ages["C"] == pytest.approx(55.0)
        assert ages["A"] == pytest.approx(0.0)  # 未带注解的节点仍按枝长推算
        assert spy.any("时间轴倒挂")


class TestTipAgesEndToEnd:
    """B-14：灭绝类群 / 带采样日期的尖端不得被画到"现在"。"""

    def test_extinct_tip_is_drawn_at_its_fossil_age_on_canvas(self):
        tree = parse_newick(ULTRAMETRIC_100)
        tree.by_name("C").add_feature("age", 55.0)  # 地层年龄 55 Ma
        figure = DatingFigure(width=8, height=6, dpi=72)
        figure.add_tree(tree).render()

        ages = _ages_by_name(figure.get_node_coordinates())
        assert ages["C"] == pytest.approx(55.0)
        assert ages["D"] == pytest.approx(0.0)
        # 真正画到画布上的尖端标记里必须出现 55 Ma 这个位置
        plotted = _plotted_x(figure.ax_tree)
        assert any(abs(x - 55.0) < 1e-6 for x in plotted)
        figure.close()

    def test_newick_age_annotation_is_used(self):
        """``[&age=55]`` / ``[&&NHX:age=55]`` 形式的注解也走真实年龄。"""
        tree = parse_newick(FOSSIL_ANNOTATED_NEWICK)
        coords = TreePlotter().plot(tree)
        ages = _ages_by_name(coords)
        assert ages["C{[&age=55]}"] == pytest.approx(55.0)
        assert ages["A"] == pytest.approx(0.0)
        assert ages["D"] == pytest.approx(30.0)  # 100 - (25 + 45)

    def test_non_ultrametric_tips_are_not_glued_to_present(self):
        coords = TreePlotter().plot(parse_newick(NON_ULTRAMETRIC))
        ages = _ages_by_name(coords)
        assert ages["A"] == pytest.approx(0.0)
        assert ages["C"] == pytest.approx(60.0)  # 而非修复前硬编码的 0

    def test_tip_date_suffix_survives_the_render(self):
        """尾型定年 / 古 DNA：``A|-0.045`` 的尖端画在 45 ka，而不是 0。"""
        tree = parse_newick("((A|-0.045:40,B:40)AB:60,(C:70,D:70)CD:30)root:0;")
        figure = DatingFigure(width=8, height=6, dpi=72)
        figure.add_tree(tree).render()
        ages = _ages_by_name(figure.get_node_coordinates())
        assert ages["A|-0.045"] == pytest.approx(0.045)
        figure.close()

    def test_tip_without_any_age_information_warns(self, monkeypatch):
        """默认按 0 绘制只允许发生在"确实没有任何信息"时，且必须告警。"""
        import phylodater.viz.tree_plot as tree_plot_module

        spy = SpyLogger()
        monkeypatch.setattr(tree_plot_module, "logger", spy)

        tree = parse_newick(
            "((A:40,B:40)AB:60,(C,D)CD:30)root:0;", missing_branch_length=None
        )
        tree.by_name("CD").dist = None
        ages = _ages_by_name(TreePlotter().plot(tree))
        assert ages["C"] == 0.0
        assert ages["A"] == pytest.approx(0.0)
        assert spy.any("只能按'现在'(0 Ma) 绘制")


class TestVizHasNoDirectEteImport:
    """A-5/B-10：可视化层不得自行 import ete3/ete4，一律经模型层取树后端。"""

    def test_viz_modules_do_not_import_ete(self):
        import re
        from pathlib import Path

        import phylodater.viz as viz_package

        pattern = re.compile(r"^\s*(?:from|import)\s+ete[34]\b", re.MULTILINE)
        source_dir = Path(viz_package.__file__).parent
        offenders = [
            path.name
            for path in sorted(source_dir.glob("*.py"))
            if pattern.search(path.read_text(encoding="utf-8"))
        ]
        assert offenders == [], f"viz 层仍在模块内直接导入 ete: {offenders}"

    def test_newick_string_goes_through_model_layer(self):
        """传 Newick 字符串时也必须走 PhylogeneticTree 的树后端入口。"""
        figure = DatingFigure(width=6, height=4, dpi=60)
        try:
            figure.add_tree(ULTRAMETRIC_100)
        except ImportError as exc:
            # 本工作区无可用 ete 后端（A-5）：错误必须来自模型层入口并给出
            # 可执行的下一步，而不是裸的 ``ModuleNotFoundError: No module named 'cgi'``
            assert "可视化用的树后端" in str(exc)
            assert "可视化层不自行导入 ete" in str(exc)
        else:  # pragma: no cover - 取决于本机是否装了 ete4/ete3
            assert figure._root_age == pytest.approx(100.0)
            figure.close()

    def test_non_tree_input_fails_loudly(self):
        with pytest.raises(TypeError, match="无法为绘图取得"):
            TreePlotter().plot(object())
