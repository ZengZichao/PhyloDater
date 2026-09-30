"""
模块十三：结果比较与报告功能测试

验证 DatingResult 构建、比较表（按拓扑对齐）、TSV 区间语义、柱状图误差线、
校准对账与摘要。

覆盖的审阅项：
  * B-15 —— 多方法比较必须按拓扑（后代叶集合）对齐，而不是按各方法自己的键字符串；
    解析不出拓扑的键必须显式标为"不可比较"，绝不被静默并排或静默当作"无分歧"。
  * C-16 —— TSV 逐列保留区间语义（HPD95 / CI95 / RANGE / NONE / NA）。
  * C-15 —— 比较柱状图按真实区间画误差线，并披露被省略/无区间的格子。
  * B-26 —— 摘要里出现"请求 N / 生效 M"的校准对账。
"""

from pathlib import Path

import pytest

from phylodater import ComparisonReporter, DatingResult
from phylodater.models import (
    CalibrationPoint,
    CIType,
    NodeAgeEstimate,
    PhylogeneticTree,
)


class TestDatingResult:
    """DatingResult 测试"""

    def test_rpt_001_single_method_result(self):
        """RPT-001: 单方法结果字段完整"""
        result = DatingResult(
            method_name="pathd8",
            run_id="test_run_001",
            dated_tree_newick="((A:1,B:1):0);",
            node_ages={
                "root": NodeAgeEstimate(mean_age=100.0),
            },
        )
        assert result.method_name == "pathd8"
        assert "root" in result.node_ages

    def test_rpt_002_empty_result(self, tmp_path):
        """RPT-002: 空结果处理"""
        reporter = ComparisonReporter(tmp_path / "report")
        assert reporter is not None


class TestComparisonReport:
    """比较报告测试"""

    @pytest.fixture
    def sample_results(self):
        return {
            "pathd8": DatingResult(
                method_name="pathd8",
                run_id="run_pathd8",
                dated_tree_newick="((A:1,B:1):0);",
                node_ages={
                    "root": NodeAgeEstimate(mean_age=100.0),
                },
            ),
            "mcmctree": DatingResult(
                method_name="mcmctree",
                run_id="run_mcmctree",
                dated_tree_newick="((A:1,B:1):0);",
                node_ages={
                    "root": NodeAgeEstimate(mean_age=105.0),
                },
            ),
        }

    def test_rpt_003_txt_table(self, sample_results, tmp_path):
        """RPT-003: 比较表 TXT 格式"""
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(sample_results)
        out_file = tmp_path / "comparison_table.txt"
        assert out_file.exists()
        content = out_file.read_text(encoding="utf-8")
        assert "pathd8" in content
        assert "mcmctree" in content

    def test_rpt_004_tsv_table(self, sample_results, tmp_path):
        """RPT-004: 比较表 TSV 格式"""
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(sample_results)
        out_file = tmp_path / "comparison_table.tsv"
        assert out_file.exists()
        content = out_file.read_text(encoding="utf-8")
        assert "\t" in content
        assert "pathd8" in content

    def test_rpt_006_summary(self, sample_results, tmp_path):
        """RPT-006: 摘要信息"""
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(sample_results)
        out_file = tmp_path / "summary.txt"
        assert out_file.exists()
        content = out_file.read_text(encoding="utf-8")
        assert "pathd8" in content or "mcmctree" in content


# ==================== B-15：按拓扑对齐（本项目核心功能的有效性） ====================


class TestTopologyAlignment:
    """多方法比较必须按节点拓扑对齐，而不是按键名字符串"""

    # 输入树：((A,B),(C,D)) 两个冠群
    TREE = PhylogeneticTree("(((A,B),(C,D)),E);")

    @staticmethod
    def _two_methods_differing_key_strings():
        """同一分支：mcmctree 用校准点名，r8s 用 node_<id>（审阅报告的实测形态）。"""
        return {
            "mcmctree": DatingResult(
                method_name="mcmctree",
                run_id="r_mcm",
                # mcmctree 的定年树把该分支标为 Mammalia
                dated_tree_newick="(((A,B)Mammalia,(C,D)Reptilia),E);",
                node_ages={
                    "Mammalia": NodeAgeEstimate(
                        mean_age=80.0,
                        ci_lower=70.0,
                        ci_upper=90.0,
                        ci_type=CIType.HPD95,
                    )
                },
            ),
            "r8s": DatingResult(
                method_name="r8s",
                run_id="r_r8s",
                # r8s 的定年树把**同一个分支**标为 node_37
                dated_tree_newick="(((A,B)node_37,(C,D)node_38),E);",
                node_ages={"node_37": NodeAgeEstimate(mean_age=95.0)},
            ),
        }

    def test_rpt_005_same_branch_two_key_strings_share_one_row(self, tmp_path):
        """RPT-005（B-15 主证）：同一节点的两种键名必须落在**同一行**并给出分歧值"""
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(self._two_methods_differing_key_strings(), tree=self.TREE)

        tsv = (tmp_path / "comparison_table.tsv").read_text(encoding="utf-8")
        lines = [line for line in tsv.splitlines() if line.strip()]
        data_lines = [
            line
            for line in lines[1:]
            if line.split("\t")[3] == "topology"  # Alignment 列
        ]
        # 两个方法、一个共同分支 → 恰好一行拓扑对齐记录
        assert len(data_lines) == 1, data_lines
        cells = data_lines[0].split("\t")
        header = lines[0].split("\t")
        as_dict = dict(zip(header, cells))
        assert as_dict["mcmctree_mean"] == "80.0"
        assert as_dict["r8s_mean"] == "95.0"
        # 旧行为是两行、每行只有一列有值；这里必须能看到 15 Ma 的分歧
        txt = (tmp_path / "comparison_table.txt").read_text(encoding="utf-8")
        assert "15.0 Ma" in txt
        assert "mcmctree vs r8s" in txt

    def test_rpt_007_positional_keys_are_never_merged(self, tmp_path):
        """RPT-007（B-15/C-41）：位置序号键不得凭字符串相等被并到同一行"""
        results = {
            "r8s": DatingResult(
                method_name="r8s",
                run_id="r1",
                # r8s 树里没有 internal_node_1 这个标签 → 无法解析成 clade
                dated_tree_newick="(((A,B),(C,D)),E);",
                node_ages={"internal_node_1": NodeAgeEstimate(mean_age=400.0)},
            ),
            "mdcat": DatingResult(
                method_name="mdcat",
                run_id="r2",
                # mdcat 确实把 internal_node_1 写在 (A,B) 上 → 可解析
                dated_tree_newick="(((A,B)internal_node_1,(C,D)),E);",
                node_ages={"internal_node_1": NodeAgeEstimate(mean_age=80.0)},
            ),
        }
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(results, tree=self.TREE)
        tsv = (tmp_path / "comparison_table.tsv").read_text(encoding="utf-8")
        rows = [line.split("\t") for line in tsv.splitlines()[1:]]
        header = tsv.splitlines()[0].split("\t")
        aligned = [r for r in rows if r[header.index("Alignment")] == "topology"]
        uncomparable = [r for r in rows if r[header.index("Comparable")] == "FALSE"]
        # 只有 mdcat 那一支被解析成拓扑行；r8s 的同名键单独成"不可比较"行
        assert len(aligned) == 1
        assert aligned[0][header.index("mdcat_mean")] == "80.0"
        assert aligned[0][header.index("r8s_mean")] == "NA"
        assert len(uncomparable) == 1
        assert uncomparable[0][header.index("r8s_mean")] == "400.0"
        summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
        assert "Uncomparable node keys" in summary

    def test_rpt_008_leaf_entries_filtered_for_named_and_positional_keys(
        self, tmp_path
    ):
        """RPT-008（B-15(2)）：叶节点过滤对 node_<id> 那批方法同样生效（fail-closed）"""
        results = {
            "mcmctree": DatingResult(
                method_name="mcmctree",
                run_id="r1",
                dated_tree_newick="(((A,B)Mammalia,(C,D)),E);",
                node_ages={
                    "Mammalia": NodeAgeEstimate(mean_age=80.0),
                    "A": NodeAgeEstimate(mean_age=0.0),  # 叶：根到叶距离
                    "GB_GCA_000252485_1": NodeAgeEstimate(mean_age=1.0),  # 短名叶
                },
            )
        }
        tree = PhylogeneticTree(
            "(((GB_GCA_000252485.1_d_Bacteria,B)Mammalia,(C,D)),E);"
        )
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(results, tree=tree)
        # 全名叶 + 登录号短名叶都必须被当成叶节点丢掉，而不是当成 0 Ma / 1 Ma 的节点
        assert reporter._is_leaf_node("A", tree) is True
        assert reporter._is_leaf_node("GB_GCA_000252485_1", tree) is True
        assert reporter._is_leaf_node("Mammalia", tree) is False
        txt = (tmp_path / "comparison_table.txt").read_text(encoding="utf-8")
        assert "GB_GCA_000252485_1" not in txt.split("Per-row source keys")[0]
        assert "leaf-node entries dropped as root-to-tip distances: 2" in txt

    def test_rpt_009_topology_id_and_depth_are_meaningful(self, tmp_path):
        """RPT-009（B-15(3)）：Topology_ID / Depth 两列在被显示的行上真的有值"""
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(self._two_methods_differing_key_strings(), tree=self.TREE)
        txt = (tmp_path / "comparison_table.txt").read_text(encoding="utf-8")
        body = txt.split("Per-row source keys")[0]
        assert "N/A" not in body  # 旧实现里 Depth 恒为 N/A
        assert "2L:" in body  # 2 个后代叶 + 叶集合哈希
        # 不再是"节点名截断 25 字符"这种假拓扑标识
        assert "Mammalia" in body
        row_line = [ln for ln in body.splitlines() if ln.startswith("Mammalia")][0]
        assert "2L:" in row_line and " topology" in row_line

    def test_rpt_010_calibration_names_aligned_via_mrca(self, tmp_path):
        """RPT-010：给出校准点后，纯校准点名也能被解析成真实 clade（不靠标签）

        mcmctree 的 FigTree 回退通路把节点年龄以**校准点名**入表，而它的定年树里
        并没有这个内部标签；r8s 则把同一分支标成 ``node_37``。不传校准点时前者只能
        是 name-only、后者是 topology，两行永远并不到一起；传了校准点，前者也能按
        MRCA 的后代叶集合定位，于是两方法真正并排。
        """
        results = {
            "mcmctree": DatingResult(
                method_name="mcmctree",
                run_id="r1",
                dated_tree_newick="(((A,B),(C,D)),E);",  # 无任何内部标签
                node_ages={"Mammalia": NodeAgeEstimate(mean_age=80.0)},
            ),
            "r8s": DatingResult(
                method_name="r8s",
                run_id="r2",
                dated_tree_newick="(((A,B)node_37,(C,D)node_38),E);",
                node_ages={"node_37": NodeAgeEstimate(mean_age=99.0)},
            ),
        }
        calibrations = [CalibrationPoint(name="Mammalia", mrca_leaf_pair=("A", "B"))]

        without = ComparisonReporter(tmp_path / "without")
        without.generate(results, tree=self.TREE)
        header, rows = self._tsv(tmp_path / "without")
        split_rows = [r for r in rows if r[header.index("Alignment")] == "topology"]
        assert len(split_rows) == 1  # 只有 r8s 那一行被按拓扑解析
        assert split_rows[0][header.index("mcmctree_mean")] == "NA"

        with_cals = ComparisonReporter(tmp_path / "with")
        with_cals.generate(results, tree=self.TREE, calibrations=calibrations)
        header, rows = self._tsv(tmp_path / "with")
        aligned = [r for r in rows if r[header.index("Alignment")] == "topology"]
        assert len(aligned) == 1
        row = dict(zip(header, aligned[0]))
        assert row["mcmctree_mean"] == "80.0"
        assert row["r8s_mean"] == "99.0"
        assert row["Depth"] == "2"  # 深度现在来自输入树的真实拓扑

    @staticmethod
    def _tsv(directory: Path):
        lines = (directory / "comparison_table.tsv").read_text().splitlines()
        return lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]

    def test_rpt_011_missing_tree_is_disclosed(self, tmp_path):
        """RPT-011：没有输入树时，报告自己说明"对齐口径退化"，不假装按拓扑比对"""
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(self._two_methods_differing_key_strings(), tree=None)
        txt = (tmp_path / "comparison_table.txt").read_text(encoding="utf-8")
        assert "no input tree was supplied" in txt
        # node_37 / Mammalia 无法再被判为同一节点：只能各占一行，且被标为不可核验
        assert "node-only" in txt or "name-only" in txt


# ==================== C-16：TSV 必须保留逐格区间语义 ====================


class TestIntervalSemanticsInTsv:
    """C-16：HPD95 / CI95 / RANGE / NONE / NA 不能压成同名不同义的一列"""

    @pytest.fixture
    def mixed_results(self):
        return {
            "mcmctree": DatingResult(
                method_name="mcmctree",
                run_id="r1",
                dated_tree_newick="(((A,B)Mammalia,(C,D)Reptilia),E);",
                node_ages={
                    "Mammalia": NodeAgeEstimate(
                        mean_age=80.0,
                        ci_lower=70.0,
                        ci_upper=90.0,
                        ci_type=CIType.HPD95,
                    ),
                    "Reptilia": NodeAgeEstimate(
                        mean_age=250.0,
                        ci_lower=240.0,
                        ci_upper=260.0,
                        ci_type=CIType.RANGE,
                    ),
                },
            ),
            "r8s": DatingResult(
                method_name="r8s",
                run_id="r2",
                dated_tree_newick="(((A,B)Mammalia,(C,D)Reptilia),E);",
                node_ages={"Mammalia": NodeAgeEstimate(mean_age=95.0)},
            ),
        }

    def test_c16_tsv_has_ci_type_columns(self, mixed_results, tmp_path):
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(mixed_results, tree=PhylogeneticTree("(((A,B),(C,D)),E);"))
        lines = (tmp_path / "comparison_table.tsv").read_text().splitlines()
        header = lines[0].split("\t")
        assert "mcmctree_ci_type" in header
        assert "r8s_ci_type" in header
        rows = [ln.split("\t") for ln in lines[1:]]
        by_node = {r[header.index("Node")]: r for r in rows}
        mammalia = dict(zip(header, by_node["Mammalia"]))
        reptilia = dict(zip(header, by_node["Reptilia"]))
        # 同一列位置上，"贝叶斯 HPD" / "均值±SD 的 RANGE" / "无区间" / "无该估计" 互不混淆
        assert mammalia["mcmctree_ci_type"] == "HPD95"
        assert mammalia["r8s_ci_type"] == "NONE"
        assert reptilia["mcmctree_ci_type"] == "RANGE"
        assert reptilia["r8s_ci_type"] == "NA"
        assert mammalia["r8s_ci_lower"] == "NA"

    def test_c16_summary_flags_mixed_interval_semantics(self, mixed_results, tmp_path):
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(mixed_results, tree=PhylogeneticTree("(((A,B),(C,D)),E);"))
        summary = (tmp_path / "summary.txt").read_text()
        assert "RANGE=1" in summary
        assert "MIXED interval semantics" in summary


# ==================== C-15：比较柱状图必须画误差线并披露缺口 ====================


class TestComparisonChart:
    """C-15：误差线来自真实区间；缺失/无区间必须显式披露"""

    @pytest.fixture
    def results(self):
        return {
            "mcmctree": DatingResult(
                method_name="mcmctree",
                run_id="r1",
                dated_tree_newick="(((A,B)Mammalia,(C,D)Reptilia),E);",
                node_ages={
                    "Mammalia": NodeAgeEstimate(
                        mean_age=80.0,
                        ci_lower=70.0,
                        ci_upper=90.0,
                        ci_type=CIType.HPD95,
                    ),
                    "Reptilia": NodeAgeEstimate(
                        mean_age=250.0,
                        ci_lower=240.0,
                        ci_upper=260.0,
                        ci_type=CIType.HPD95,
                    ),
                },
            ),
            "treepl": DatingResult(
                method_name="treepl",
                run_id="r2",
                dated_tree_newick="(((A,B)Mammalia,(C,D)Reptilia),E);",
                node_ages={"Mammalia": NodeAgeEstimate(mean_age=95.0)},
            ),
        }

    def test_c15_error_bars_drawn_and_gaps_disclosed(
        self, results, tmp_path, monkeypatch
    ):
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.axes
        import numpy as np

        captured = []
        real_bar = matplotlib.axes.Axes.bar

        def spy(self, x, height, *args, **kwargs):
            captured.append(kwargs)
            return real_bar(self, x, height, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "bar", spy)

        reporter = ComparisonReporter(tmp_path)
        info = reporter.generate(results, tree=PhylogeneticTree("(((A,B),(C,D)),E);"))

        assert captured, "bar chart was not drawn"
        yerr = [np.asarray(call["yerr"]) for call in captured if "yerr" in call]
        assert yerr, "no yerr passed to ax.bar (C-15 regression)"
        assert any((arr > 0).any() for arr in yerr), "error bars are all zero"
        stats = info["plot"]
        assert stats["error_bars_per_method"]["mcmctree"] == 2
        assert stats["error_bars_per_method"]["treepl"] == 0
        # treepl 没有 Reptilia 这一格：柱子成 NaN，必须被数出来而不是静默消失
        assert stats["missing_cells"] == 1
        assert stats["cells_without_interval"] == 1
        assert (tmp_path / "comparison_plot.png").exists()

    def test_c17_visualization_failure_is_a_warning_not_debug(
        self, results, tmp_path, monkeypatch
    ):
        """C-17：图生成失败必须至少 warning，且成功语列出实际产物清单"""
        logged = []
        reporter = ComparisonReporter(tmp_path)
        # 单例 logger：只在本用例内替换 warning，用例结束自动还原
        monkeypatch.setattr(
            reporter.logger, "warning", lambda msg, **kw: logged.append(str(msg))
        )
        monkeypatch.setattr(
            ComparisonReporter,
            "generate_visualization",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")),
        )

        info = reporter.generate(results)

        assert any("boom" in message for message in logged), logged
        assert info["failures"]
        assert not (tmp_path / "comparison_plot.png").exists()
        # 仍然要产出表格，并且成功语列出实际产出的文件
        assert (tmp_path / "comparison_table.tsv").exists()
        assert "comparison_table.tsv" in " ".join(info["files"])


# ==================== B-26：校准"请求 N / 生效 M"必须在报告里可见 ====================


class TestCalibrationReconciliation:
    """B-26：适配器已经暴露在 result.metadata 里的对账信息，报告必须复述"""

    def test_b26_requested_versus_effective_in_summary(self, tmp_path):
        results = {
            "pathd8": DatingResult(
                method_name="pathd8",
                run_id="r1",
                dated_tree_newick="(((A,B)Mammalia,(C,D)),E);",
                node_ages={"Mammalia": NodeAgeEstimate(mean_age=80.0)},
                metadata={
                    "constraint_reconciliation": {
                        "calibrations_supplied": 5,
                        "constrained_nodes_written": 3,
                        "calibrations_not_written": ["OldNode", "RootNode"],
                        "directives_written": {"fixage": 2, "minage": 1},
                    },
                    "clock_tests": {"accepted": 3, "total": 5},
                },
            ),
            "mcmctree": DatingResult(
                method_name="mcmctree",
                run_id="r2",
                dated_tree_newick="(((A,B)Mammalia,(C,D)),E);",
                node_ages={
                    "Mammalia": NodeAgeEstimate(
                        mean_age=80.0,
                        ci_lower=70.0,
                        ci_upper=90.0,
                        ci_type=CIType.RANGE,
                    )
                },
                metadata={
                    "node_age_interval_kinds": {
                        "Mammalia": "多链外包络（RANGE），不是单链 HPD95"
                    },
                    "chain_interval_conflicts": {"Mammalia": [[1, 2]]},
                },
            ),
        }
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(results, tree=PhylogeneticTree("(((A,B),(C,D)),E);"))
        summary = (tmp_path / "summary.txt").read_text()
        assert "Calibration reconciliation" in summary
        assert "requested=5 / effective=3" in summary
        assert "['OldNode', 'RootNode']" in summary
        assert "clock tests accepted=3/5" in summary
        assert "chain interval conflicts on 1 node(s)" in summary
        # 1 个节点被算出，但请求了 5 个校准 —— 缩水必须看得见
        assert "only 1 node age(s) parsed for 5 requested" in summary

    def test_b26_methods_without_reconciliation_say_so(self, tmp_path):
        results = {
            "lsd2": DatingResult(
                method_name="lsd2",
                run_id="r1",
                dated_tree_newick="(((A,B)Mammalia,(C,D)),E);",
                node_ages={"Mammalia": NodeAgeEstimate(mean_age=80.0)},
            )
        }
        reporter = ComparisonReporter(tmp_path)
        reporter.generate(results, tree=PhylogeneticTree("(((A,B),(C,D)),E);"))
        summary = (tmp_path / "summary.txt").read_text()
        assert "no reconciliation metadata published by this adapter" in summary
