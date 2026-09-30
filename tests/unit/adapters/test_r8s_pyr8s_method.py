"""
R8sPyr8sMethod 适配器单元测试

测试 r8s/pyr8s 定年方法的各项功能

除了配置/语法/解析断言外，还覆盖 C-41 的纪律：任何"节点从结果里消失"的路径
都必须点名报出来并随 ``DatingResult.warnings`` 落盘，不允许静默丢数据。
"""

import json
import shutil
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.adapters.r8s_pyr8s_method import R8sPyr8sMethod
from phylodater.core.exceptions import SemanticDegradationWarning
from phylodater.infrastructure.configuration import CommonConfig, R8sConfig
from phylodater.models import (
    CalibrationPoint,
    FixedAgeConstraint,
    PhylogeneticTree,
    UniformAgeConstraint,
)


class TestR8sPyr8sMethod:
    """R8sPyr8sMethod 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="r8s_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def r8s_config(self):
        """创建 r8s 配置"""
        return R8sConfig(method="PL", algorithm="TN")

    @pytest.fixture
    def common_config(self):
        """创建通用配置"""
        return CommonConfig(nthreads=4)

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        tree.num_tips = 4
        tree.without_internal_labels = Mock(return_value=tree)
        tree.write = Mock()
        return tree

    @pytest.fixture
    def mock_calibrations(self):
        """创建模拟校准点"""
        calib1 = Mock(spec=CalibrationPoint)
        calib1.name = "Primates"
        calib1.resolved_taxa = ["human", "chimp"]
        calib1.mrca_leaf_pair = ("human", "chimp")
        calib1.is_root_node = False
        calib1.age_constraint = UniformAgeConstraint(min_age=6.0, max_age=8.0)

        return [calib1]

    def test_init_creates_work_directory(self, temp_dir, r8s_config, common_config):
        """测试初始化创建工作目录"""
        method = R8sPyr8sMethod(r8s_config, temp_dir, common_config=common_config)

        assert method.work_dir.exists()
        assert method.method_name == "r8s"

    def test_method_name_property(self, temp_dir, r8s_config, common_config):
        """测试方法名称属性"""
        method = R8sPyr8sMethod(r8s_config, temp_dir, common_config=common_config)
        assert method.method_name == "r8s"

    def test_common_config_usage(self, temp_dir, r8s_config, common_config):
        """测试 common_config 使用"""
        method = R8sPyr8sMethod(r8s_config, temp_dir, common_config=common_config)
        assert method.common_config.nthreads == 4

    def test_config_parameters(self, temp_dir, r8s_config, common_config):
        """测试配置参数"""
        method = R8sPyr8sMethod(r8s_config, temp_dir, common_config=common_config)

        assert method.config.method == "PL"
        assert method.config.algorithm == "TN"

    @patch("phylodater.adapters.r8s_pyr8s_method.ProcessRunner")
    def test_validate_environment_success(
        self, mock_runner_class, temp_dir, r8s_config, common_config
    ):
        """测试环境验证成功"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = "1.0"
        mock_runner_class.return_value = mock_runner

        method = R8sPyr8sMethod(r8s_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is True

    @patch("phylodater.adapters.r8s_pyr8s_method.ProcessRunner")
    @patch.dict("sys.modules", {"pyr8s": None, "pyr8s.core": None})
    def test_validate_environment_failure(
        self, mock_runner_class, temp_dir, r8s_config, common_config
    ):
        """测试环境验证失败"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = False
        mock_runner_class.return_value = mock_runner

        method = R8sPyr8sMethod(r8s_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is False

    def test_method_pl(self, temp_dir, common_config):
        """测试 PL 方法"""
        config = R8sConfig(method="PL")
        method = R8sPyr8sMethod(config, temp_dir, common_config=common_config)
        assert method.config.method == "PL"

    def test_method_nprs(self, temp_dir, common_config):
        """测试 NPRS 方法"""
        config = R8sConfig(method="NPRS")
        method = R8sPyr8sMethod(config, temp_dir, common_config=common_config)
        assert method.config.method == "NPRS"

    def test_work_dir_isolation(self, temp_dir, r8s_config, common_config):
        """测试工作目录隔离"""
        method1 = R8sPyr8sMethod(r8s_config, temp_dir, common_config=common_config)
        method2 = R8sPyr8sMethod(
            r8s_config, temp_dir / "other", common_config=common_config
        )

        assert method1.work_dir != method2.work_dir


class TestR8sConfig:
    """r8s 配置测试"""

    def test_default_config(self):
        """测试默认配置"""
        config = R8sConfig()
        # 默认 NPRS：pyr8s 后端仅支持 NPRS，默认 PL 会被静默替换(F14)
        assert config.method == "NPRS"
        assert config.algorithm == "TN"

    def test_custom_config(self):
        """测试自定义配置"""
        config = R8sConfig(method="NPRS", smoothing=10.0)
        assert config.method == "NPRS"
        assert config.smoothing == 10.0

    def test_method_validation(self):
        """测试方法验证"""
        with pytest.raises(ValueError):
            R8sConfig(method="INVALID")

    def test_smoothing_validation(self):
        """测试 smoothing 验证"""
        with pytest.raises(ValueError):
            R8sConfig(smoothing=-1.0)


class TestR8sOriginalSubprocessGrammar:
    """A-6：原版 r8s 子进程通路的调用形式与输入语法必须匹配 Sanderson r8s 1.7 手册。

    这些断言针对上游手册里真实存在的命令（blformat / set smoothing / divtime
    method=…algorithm=… / showage shownamed=yes / describe plot=chrono_description），
    而不是适配器自造的 blength / method= / algorithm= / smoothing=。
    """

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="r8s_grammar_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def _make(self, temp_dir):
        method = R8sPyr8sMethod(
            R8sConfig(method="PL", smoothing=50.0),
            temp_dir,
            common_config=CommonConfig(nthreads=2),
        )
        method._is_pyr8s = False
        return method

    def _mock_tree(self):
        tree = Mock(spec=PhylogeneticTree)
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        tree.newick = "((human:0.1,chimp:0.1)A:0.2,(mouse:0.3,rat:0.3)B:0.1);"
        tree.with_renamed_leaves = Mock(return_value=tree)
        tree.without_internal_labels = Mock(return_value=tree)
        return tree

    def test_prepare_inputs_writes_upstream_grammar(self, temp_dir):
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        method = self._make(temp_dir)
        root = CalibrationPoint(
            name="Root",
            age_constraint=FixedAgeConstraint(fixed_age=500.0),
            resolved_taxa=["human", "chimp", "mouse", "rat"],
            mrca_leaf_pair=None,
            is_root_node=True,
        )
        prim = CalibrationPoint(
            name="Primates",
            age_constraint=FixedAgeConstraint(fixed_age=80.0),
            resolved_taxa=["human", "chimp"],
            mrca_leaf_pair=("human", "chimp"),
            is_root_node=False,
        )
        method.prepare_inputs(self._mock_tree(), [root, prim], alignment_path=None)
        nex = (method.work_dir / "r8s_input.nex").read_text()

        assert "blformat lengths=persite nsites=1000 ultrametric=no;" in nex
        assert "set smoothing=50.0;" in nex
        assert "divtime method=PL algorithm=TN;" in nex
        assert "showage shownamed=yes;" in nex
        assert "describe plot=chrono_description;" in nex
        # 上游根本不存在的自造指令必须消失
        assert "blength nsites" not in nex
        assert "method = " not in nex
        assert "algorithm = " not in nex
        assert "smoothing = " not in nex
        assert "plot=chronogram" not in nex

    def test_root_calibration_does_not_crash_missing_leaf_pair(self, temp_dir):
        """B-20：根校准 mrca_leaf_pair=None 时不得对 None 下标。"""
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        method = self._make(temp_dir)
        root = CalibrationPoint(
            name="Root",
            age_constraint=FixedAgeConstraint(fixed_age=500.0),
            resolved_taxa=["human", "chimp", "mouse", "rat"],
            mrca_leaf_pair=None,
            is_root_node=True,
        )
        no_pair = CalibrationPoint(
            name="Loose",
            age_constraint=FixedAgeConstraint(fixed_age=10.0),
            resolved_taxa=[],
            mrca_leaf_pair=None,
            is_root_node=False,
        )
        # 不得抛 TypeError/IndexError；缺叶对的校准应被跳过
        method.prepare_inputs(self._mock_tree(), [root, no_pair], alignment_path=None)
        nex = (method.work_dir / "r8s_input.nex").read_text()
        assert "MRCA Root" in nex
        assert "MRCA Loose" not in nex

    def test_execute_subprocess_invokes_batch_flags(self, temp_dir):
        method = self._make(temp_dir)
        nex = method.work_dir / "r8s_input.nex"
        nex.write_text("#NEXUS\n")
        method._input_files = {"nexus": nex}
        captured = {}

        def fake_run(cmd):
            captured["cmd"] = list(cmd)
            result = Mock()
            result.returncode = 0
            result.stdout = ""
            result.stderr = ""
            return result

        with patch("phylodater.adapters.r8s_pyr8s_method.ProcessRunner") as PR:
            PR.return_value.run.side_effect = fake_run
            assert method._execute_with_subprocess() is True

        # 手册：r8s -f <file> -b （批处理，避免等待键盘输入）
        assert captured["cmd"] == ["r8s", "-f", str(nex), "-b"]


class TestR8sOriginalSubprocessParsing:
    """B-25：从 chrono_description 按拓扑读年龄；空/冲突一律显式失败，不再以 PASS 返回空。"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="r8s_parse_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def _make(self, temp_dir, stdout_text):
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        method = R8sPyr8sMethod(
            R8sConfig(method="PL"), temp_dir, common_config=CommonConfig(nthreads=2)
        )
        method._is_pyr8s = False
        method._use_api = False
        root = CalibrationPoint(
            name="Root",
            age_constraint=FixedAgeConstraint(fixed_age=500.0),
            resolved_taxa=["human", "chimp", "mouse", "rat"],
            mrca_leaf_pair=None,
            is_root_node=True,
        )
        prim = CalibrationPoint(
            name="Primates",
            age_constraint=FixedAgeConstraint(fixed_age=80.0),
            resolved_taxa=["human", "chimp"],
            mrca_leaf_pair=("human", "chimp"),
            is_root_node=False,
        )
        method._calibrations = [root, prim]
        method._tip_name_map = {t: t for t in ["human", "chimp", "mouse", "rat"]}
        outp = temp_dir / "r8s_output.txt"
        outp.write_text(stdout_text)
        method._input_files = {"output": outp}
        return method

    def test_reads_true_age_from_chrono_tree_not_first_number(self, temp_dir):
        stdout = (
            "r8s> divtime\n"
            "showage shownamed=yes:\n"
            "Name\tType\tAge\tRate\n"
            "--------------------------------\n"
            "Primates\t1\t80.000000\t1.2345e-02\n"
            "Root\t1\t500.000000\t1.2000e-02\n"
            "\n"
            'tree "chronogram" = '
            "((human:80.0,chimp:80.0)Primates:420.0,"
            "(mouse:200.0,rat:200.0)Rodent:300.0)Root:0.0;\n"
        )
        method = self._make(temp_dir, stdout)
        result = method.parse_results()
        # 关键：80 Ma 必须读成 80，绝不再被行内第一个数字 1 顶替
        assert abs(result.node_ages["Primates"].mean_age - 80.0) < 1e-6
        assert abs(result.node_ages["Root"].mean_age - 500.0) < 1e-6
        assert "Primates" in result.dated_tree_newick

    def test_empty_output_raises_instead_of_passing(self, temp_dir):
        from phylodater.core.exceptions import ResultParsingError

        method = self._make(temp_dir, "r8s> nothing useful\nno importable tree here\n")
        with pytest.raises(ResultParsingError):
            method.parse_results()

    def test_validate_node_ages_raises_on_empty(self, temp_dir):
        from phylodater.core.exceptions import ResultParsingError

        method = self._make(temp_dir, "x")
        with pytest.raises(ResultParsingError):
            method._validate_node_ages({})

    def test_showage_fallback_exact_field_match(self, temp_dir):
        method = self._make(temp_dir, "x")
        ages = method._parse_showage_named_ages(
            "Name Type Age Rate\nPrimates 1 80.000000 1.2345e-02\n"
        )
        assert abs(ages["Primates"].mean_age - 80.0) < 1e-9

    def test_showage_fallback_rejects_conflicting_duplicates(self, temp_dir):
        from phylodater.core.exceptions import ResultParsingError

        method = self._make(temp_dir, "x")
        with pytest.raises(ResultParsingError):
            method._parse_showage_named_ages(
                "Primates 1 80.000000 1.2e-02\nPrimates 1 99.000000 1.2e-02\n"
            )


class TestC41DroppedInternalNodes:
    """C-41：年龄算成 <= 0 的内部节点不得**静默**从结果里消失。

    旧实现是 ``if age <= 0: continue``——比较表少一行、日志却报成功，无人知道少了谁。
    现在仍然不把非正年龄当成年龄印出来，但必须点名披露（按标签/后代叶集合，而不是
    遍历序号），并且随结果落盘成 ``DatingResult.warnings``。
    """

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="r8s_c41_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def _make(self, temp_dir):
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        method = R8sPyr8sMethod(
            R8sConfig(method="PL"), temp_dir, common_config=CommonConfig(nthreads=2)
        )
        method._is_pyr8s = False
        method._use_api = False
        method._calibrations = [
            CalibrationPoint(
                name="Root",
                age_constraint=FixedAgeConstraint(fixed_age=500.0),
                resolved_taxa=["human", "chimp", "mouse", "rat"],
                mrca_leaf_pair=None,
                is_root_node=True,
            ),
            CalibrationPoint(
                name="Primates",
                age_constraint=FixedAgeConstraint(fixed_age=80.0),
                resolved_taxa=["human", "chimp"],
                mrca_leaf_pair=("human", "chimp"),
                is_root_node=False,
            ),
        ]
        method._tip_name_map = {t: t for t in ["human", "chimp", "mouse", "rat"]}
        return method

    @staticmethod
    def _chronogram(newick):
        from dendropy import Tree

        return Tree.get(
            data=newick,
            schema="newick",
            preserve_underscores=True,
            rooting="force-rooted",
        )

    #: Rodent 这一支的两条枝长都是 0 ⇒ 该未命名内部节点的年龄算成 0
    ZERO_CLADE_NWK = (
        "((human:80.0,chimp:80.0)Primates:420.0,(mouse:0.0,rat:0.0)Rodent:500.0)"
        "Root:0.0;"
    )

    def test_zero_age_node_is_named_in_a_warning(self, temp_dir):
        from phylodater.core.exceptions import SemanticDegradationWarning

        method = self._make(temp_dir)
        node_ages = method._extract_ages_from_chronogram(
            self._chronogram(self.ZERO_CLADE_NWK)
        )

        # 校准节点照旧产出
        assert {"Root", "Primates"} <= set(node_ages)
        # 零年龄节点不产出"假年龄"，但必须被点名
        assert not any(estimate.mean_age <= 0 for estimate in node_ages.values())
        assert len(method._warnings) == 1
        message = str(method._warnings[0])
        assert isinstance(method._warnings[0], SemanticDegradationWarning)
        assert "excluded from the results" in message
        assert "1 unnamed internal node(s)" in message
        assert "label='Rodent'" in message  # 有标签就给标签
        assert "mouse" in message and "rat" in message  # 并且给后代叶集合
        assert "age=0" in message

    def test_unlabelled_zero_age_node_is_named_by_its_leaf_set(self, temp_dir):
        """未命名内部节点只能靠后代叶集合点名——位置序号 ``internal_node_1`` 不算。"""
        method = self._make(temp_dir)
        method._extract_ages_from_chronogram(
            self._chronogram(
                "((human:80.0,chimp:80.0):420.0,(mouse:0.0,rat:0.0):500.0);"
            )
        )

        assert len(method._warnings) == 1
        message = str(method._warnings[0])
        assert "leaves=[mouse, rat]" in message
        assert "internal_node_" not in message

    def test_warning_reaches_parsed_result(self, temp_dir):
        """只留在日志里等于没有：警告必须挂到 DatingResult 上随产物落盘。"""
        method = self._make(temp_dir)
        output = temp_dir / "r8s_output.txt"
        output.write_text(
            'tree "chronogram" = ' + self.ZERO_CLADE_NWK + "\n", encoding="utf-8"
        )
        method._input_files = {"output": output}

        result = method.parse_results()

        assert [str(w) for w in result.warnings] == [str(w) for w in method._warnings]
        assert any("excluded from the results" in str(w) for w in result.warnings)
        serialized = json.dumps(result.to_dict())
        assert "excluded from the results" in serialized

    def test_all_positive_ages_produce_no_warning(self, temp_dir):
        method = self._make(temp_dir)
        method._extract_ages_from_chronogram(
            self._chronogram(
                "((human:80.0,chimp:80.0)Primates:420.0,"
                "(mouse:200.0,rat:200.0)Rodent:300.0)Root:0.0;"
            )
        )
        assert method._warnings == []

    def test_whole_extraction_failure_is_no_longer_debug_only(self, temp_dir):
        """整段提取失败 = 所有未命名节点都没有年龄，同样不能只记 debug。"""
        method = self._make(temp_dir)
        method._calibrations = []

        ages = method._extract_ages_from_chronogram(object())  # 不得抛出
        assert ages == {}
        assert len(method._warnings) == 1
        assert "could not extract unnamed internal node ages" in str(
            method._warnings[0]
        )

    def test_warnings_are_reset_per_run(self, temp_dir):
        method = self._make(temp_dir)
        tree = Mock(spec=PhylogeneticTree)
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        tree.newick = (
            "((human:0.1,chimp:0.1)Primates:0.2,(mouse:0.3,rat:0.3)Rodent:0.1);"
        )
        tree.with_renamed_leaves = Mock(return_value=tree)
        tree.without_internal_labels = Mock(return_value=tree)
        method._warnings.append(SemanticDegradationWarning("上一次运行的残留"))

        method.prepare_inputs(tree, method._calibrations, alignment_path=None)

        assert method._warnings == []


class TestPyr8sAdvisoryRanges:
    """pyr8s 后端不强制 CONSTRAIN 区间 —— 必须预先披露并在事后对账。

    已用 pyr8s 0.3.1 直接实验确认（scripts/probe_pyr8s_constraints.py）：内部节点上的
    ``constrain taxon=X min_age=68; max_age=82;`` 会被 NPRS 忽略（实测返回 54.6 Ma），
    而 ``fixage`` 会被精确执行。本项目"绝不静默丢约束"的纪律要求：写输入文件前说明
    区间只是参考值，解析结果后逐条对账并把未满足的校准点名写进结果级警告。
    """

    @pytest.fixture
    def temp_dir(self):
        path = Path(tempfile.mkdtemp(prefix="r8s_advisory_"))
        yield path
        shutil.rmtree(path, ignore_errors=True)

    def _pyr8s_backend(self, temp_dir):
        method = R8sPyr8sMethod(
            R8sConfig(method="NPRS"), temp_dir, common_config=CommonConfig()
        )
        method._is_pyr8s = True
        method._use_api = True
        return method

    def test_range_calibrations_are_disclosed_before_the_run(self, temp_dir):
        method = self._pyr8s_backend(temp_dir)
        method._warn_pyr8s_advisory_ranges(
            [
                CalibrationPoint(
                    name="Primates",
                    mrca_leaf_pair=("human", "chimp"),
                    age_constraint=UniformAgeConstraint(min_age=68.0, max_age=82.0),
                ),
                CalibrationPoint(
                    name="Rodentia",
                    mrca_leaf_pair=("mouse", "rat"),
                    age_constraint=FixedAgeConstraint(fixed_age=12.0),
                ),
            ]
        )
        text = " ".join(str(w) for w in method._warnings)
        assert "ADVISORY" in text
        assert "Primates" in text and "Rodentia" not in text  # 只点名非点校准

    def test_point_only_calibrations_produce_no_disclosure(self, temp_dir):
        method = self._pyr8s_backend(temp_dir)
        method._warn_pyr8s_advisory_ranges(
            [
                CalibrationPoint(
                    name="Rodentia",
                    mrca_leaf_pair=("mouse", "rat"),
                    age_constraint=FixedAgeConstraint(fixed_age=12.0),
                )
            ]
        )
        assert method._warnings == []

    def test_violated_range_is_reported_after_the_run(self, temp_dir):
        from phylodater.models import CIType, NodeAgeEstimate

        method = self._pyr8s_backend(temp_dir)
        method._calibrations = [
            CalibrationPoint(
                name="Primates",
                mrca_leaf_pair=("human", "chimp"),
                age_constraint=UniformAgeConstraint(min_age=68.0, max_age=82.0),
            ),
            CalibrationPoint(
                name="Rodentia",
                mrca_leaf_pair=("mouse", "rat"),
                age_constraint=FixedAgeConstraint(fixed_age=12.0),
            ),
        ]
        method._verify_requested_constraints(
            {
                "Primates": NodeAgeEstimate(
                    mean_age=54.64, ci_type=CIType.NONE
                ),  # < min 68 → 未满足
                "Rodentia": NodeAgeEstimate(mean_age=12.0, ci_type=CIType.NONE),
            }
        )
        text = " ".join(str(w) for w in method._warnings)
        assert "'Primates' was NOT satisfied" in text
        assert "54.64 Ma" in text and "68.0" in text
        assert "Rodentia" not in text  # 被满足的校准不该被点名

    def test_r8s_binary_backend_does_not_claim_ranges_are_advisory(self, temp_dir):
        """原版 r8s 二进制路径不作该声明：结论只对 pyr8s 的 NPRS 实测有效。"""
        method = R8sPyr8sMethod(
            R8sConfig(method="PL"), temp_dir, common_config=CommonConfig()
        )
        method._is_pyr8s = False
        method._warn_pyr8s_advisory_ranges(
            [
                CalibrationPoint(
                    name="Primates",
                    mrca_leaf_pair=("human", "chimp"),
                    age_constraint=UniformAgeConstraint(min_age=68.0, max_age=82.0),
                )
            ]
        )
        assert method._warnings == []
