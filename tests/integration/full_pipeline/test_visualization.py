#!/usr/bin/env python3
"""
PhyloDater 时间树绘制功能测试

集成口径：走 ``PhylogeneticTree`` -> 模型层树后端 -> ``DatingFigure`` 的完整链路，
并覆盖审阅项 B-13（``root_age`` 缺省时按分支长度推断，而不是当成 0）与 B-14
（叶/尖端年龄不再硬编码为 0）。

模型层目前只提供 ete3 后端，而 ete3 3.1.x 在本工作区的解释器（Python 3.14）上
根本无法导入（依赖已被移除的 stdlib ``cgi``，见审阅报告 A-5）。因此凡是需要真实
``PhylogeneticTree`` 树后端的用例，只在后端缺失时 ``skip`` 并写明原因（属
``models/tree.py`` 的范围），任何来自可视化层自身的异常都仍会让测试失败；同一批
断言另外用"最小 ete 形状"的树替身再跑一遍，保证 B-13/B-14 在任意环境下都有
真实的出图证据。
"""

import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 非交互式后端

import pytest  # noqa: E402

from phylodater.core.pipeline import DatingPipeline, PipelineConfig  # noqa: E402
from phylodater.models import (  # noqa: E402
    DatingResult,
    NodeAgeEstimate,
    PhylogeneticTree,
)
from phylodater.models.results import CIType  # noqa: E402
from phylodater.viz import DatingFigure, TreePlotter  # noqa: E402
from phylodater.viz.tree_plot import as_ete_tree, compute_node_coordinates  # noqa: E402

# 最小 ete 形状树替身（可视化层按鸭子类型访问节点，理由见模块 docstring）
from tests.unit.test_viz_geo import parse_newick  # noqa: E402

NEWICK = (
    "((human:6,(chimp:2,bonobo:2):4):9,"
    "(gorilla:10,(orangutan:7,sumatran:7):3):5)root:0;"
)
TIP_NAMES = ("human", "chimp", "bonobo", "gorilla", "orangutan", "sumatran")

#: 树自身分支长度隐含的根高度（最长 root-to-tip 路径）：
#: human 6+9、chimp/bonobo 2+4+9、gorilla 10+5、orangutan/sumatran 7+3+5 都等于 15
ROOT_AGE_OF_NEWICK = 15.0
N_NODES = 11  # 6 尖端 + 4 内部 + 根


def _ages_by_name(coords):
    return {coord.name: coord.age for coord in coords.values()}


def _plot_coords(tree):
    """``TreePlotter().plot(tree)`` 的简写（自建 figure，画完即关）。"""
    return TreePlotter().plot(tree)


def _model_backend_available() -> bool:
    """模型层能否给出可视化可用的（ete 形状）树后端。"""
    try:
        as_ete_tree(PhylogeneticTree.from_newick(NEWICK))
    except (ImportError, TypeError):
        return False
    return True


requires_model_backend = pytest.mark.skipif(
    not _model_backend_available(),
    reason=(
        "models/tree.py 只提供 ete3 后端，而 ete3 在本机解释器上无法导入"
        "（审阅项 A-5，属模型层范围，非 viz 层问题）"
    ),
)


def test_tree_plotter():
    """测试 TreePlotter：坐标必须由分支长度真实推出来。"""
    tree = PhylogeneticTree(_canonical_newick=NEWICK)
    assert tree is not None
    assert TreePlotter() is not None

    # 真实模型对象这一段只在树后端可用时才跑（A-5）
    try:
        backend = as_ete_tree(tree)
    except ImportError as exc:
        assert "可视化用的树后端" in str(exc), "错误应来自模型层入口"
    else:
        coords = compute_node_coordinates(backend, root_age=ROOT_AGE_OF_NEWICK)
        assert len(coords) == N_NODES
        ages = _ages_by_name(coords)
        assert ages["root"] == pytest.approx(ROOT_AGE_OF_NEWICK)
        assert all(ages[name] == pytest.approx(0.0) for name in TIP_NAMES)

    # 同一套断言跑在鸭子类型的树替身上：B-13 的推断与 ete 无关
    coords = compute_node_coordinates(parse_newick(NEWICK))
    assert len(coords) == N_NODES
    ages = _ages_by_name(coords)
    assert ages["root"] == pytest.approx(ROOT_AGE_OF_NEWICK)
    assert ages["root"] != 0.0, "B-13：未给 root_age 时不得当成 0"
    assert all(ages[name] == pytest.approx(0.0) for name in TIP_NAMES)


def test_dating_figure(tmp_path):
    """测试 DatingFigure 完整流程（含 B-13 / B-14 与真实存盘产物）。"""
    result = DatingResult(
        method_name="test",
        run_id="test_run_1",
        dated_tree_newick=NEWICK,
        node_ages={
            "root": NodeAgeEstimate(
                mean_age=ROOT_AGE_OF_NEWICK,
                median_age=ROOT_AGE_OF_NEWICK,
                ci_lower=14.0,
                ci_upper=16.0,
                ci_type=CIType.HPD95,
            ),
            "human_chimp": NodeAgeEstimate(
                mean_age=6.0,
                median_age=6.0,
                ci_lower=5.0,
                ci_upper=7.0,
                ci_type=CIType.HPD95,
            ),
        },
    )

    output_path = tmp_path / "dating_figure.png"

    # (1) 真实模型对象 + 显式 root_age + 定年结果：既有契约不得破坏
    figure = DatingFigure(width=12, height=8, dpi=150)
    try:
        figure.add_tree(PhylogeneticTree(_canonical_newick=NEWICK), root_age=15)
    except ImportError as exc:
        assert "可视化用的树后端" in str(exc), "错误应来自模型层入口"
    else:
        figure.add_result(result).add_geo_scale().render()
        assert figure.save(output_path).exists()
        assert output_path.stat().st_size > 0
        ages = _ages_by_name(figure.get_node_coordinates())
        assert ages["root"] == pytest.approx(ROOT_AGE_OF_NEWICK)
        assert figure.ax_tree.get_xlim()[0] > figure.ax_tree.get_xlim()[1]
        assert len(figure.ax_geo.patches) > 0
        figure.close()
    output_path.unlink(missing_ok=True)

    # (2) 模块文档里那段"不传 root_age"的用法：根年龄按分支长度推断（B-13）
    stub_figure = (
        DatingFigure(width=12, height=8, dpi=150)
        .add_tree(parse_newick(NEWICK))
        .add_geo_scale()
        .render()
    )
    assert stub_figure._root_age == pytest.approx(ROOT_AGE_OF_NEWICK)
    stub_ages = _ages_by_name(stub_figure.get_node_coordinates())
    assert stub_ages["root"] == pytest.approx(ROOT_AGE_OF_NEWICK)
    assert stub_ages["chimp"] == pytest.approx(0.0)
    assert stub_ages["sumatran"] == pytest.approx(0.0)
    assert stub_figure.save(output_path).exists()
    assert output_path.stat().st_size > 0
    stub_figure.close()

    # (3) B-14：灭绝类群 / 带采样日期的尖端用树里记录的年龄，不再钉在"现在"
    fossil = parse_newick(NEWICK)
    fossil.by_name("sumatran").add_feature("age", 0.5)
    fossil_ages = _ages_by_name(_plot_coords(fossil))
    assert fossil_ages["sumatran"] == pytest.approx(0.5)
    assert fossil_ages["gorilla"] == pytest.approx(0.0)

    dated = parse_newick(NEWICK)
    dated.by_name("human").add_feature("date", -0.045)  # LSD2 / 尾型定年口径
    dated_ages = _ages_by_name(_plot_coords(dated))
    assert dated_ages["human"] == pytest.approx(0.045)


@requires_model_backend
def test_dating_figure_with_model_tree_object(tmp_path):
    """真实 ``PhylogeneticTree`` -> 树后端 -> 出图（后端缺失时 skip，见 A-5）。"""
    tree = PhylogeneticTree.from_newick(NEWICK)
    figure = DatingFigure(width=10, height=7, dpi=120)
    figure.add_tree(tree).add_geo_scale().render()  # 不传 root_age：按枝长推断

    assert figure._root_age == pytest.approx(ROOT_AGE_OF_NEWICK)
    ages = _ages_by_name(figure.get_node_coordinates())
    assert ages["root"] == pytest.approx(ROOT_AGE_OF_NEWICK)
    assert ages["bonobo"] == pytest.approx(0.0)

    output = tmp_path / "from_model.pdf"
    assert figure.save(output).exists()
    assert output.stat().st_size > 1024
    figure.close()


def test_geo_plotter():
    """测试 GeoPlotter"""
    from phylodater.viz.geo_plot import GeoPlotter, load_ics_data

    geo_data = load_ics_data()
    assert len(geo_data) > 0

    plotter = GeoPlotter()
    assert plotter is not None


def test_pipeline_flow():
    """测试 Pipeline 完整流程"""
    config = PipelineConfig(
        methods=["pathd8"],
        output_dir=Path(tempfile.mkdtemp()),
        threads=1,
        dry_run=True,
    )
    assert config is not None

    pipeline = DatingPipeline(config)
    assert pipeline is not None
