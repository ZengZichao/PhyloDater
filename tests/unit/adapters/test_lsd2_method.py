"""
LSD2Method 适配器单元测试

测试 LSD2 (IQ-TREE2) 定年方法的各项功能：
- 配置 / 环境检测 / 工作目录
- **输入生成路径**（constraints.date 文件的逐行内容）：区间约束、单边约束的
  NA 写法、降级警告、叶节点真实采样日期的保留（B-16 / B-17 回归防护）
"""

import shutil
import tempfile
import warnings
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.adapters.lsd2_method import LSD2Method
from phylodater.core.exceptions import CalibrationError, SemanticDegradationWarning
from phylodater.infrastructure.configuration import CommonConfig, LSD2Config
from phylodater.models import (
    CalibrationPoint,
    FixedAgeConstraint,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    PhylogeneticTree,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)

# 一棵根分裂为 (A,B) | (C,D) 的树；dendropy 前序叶名即 A,B,C,D
TEST_NEWICK = "((A:0.1,B:0.1)AB:0.2,(C:0.3,D:0.3)CD:0.1);"


def _cal(name, constraint, pair=None, is_root=False):
    return CalibrationPoint(
        name=name,
        age_constraint=constraint,
        resolved_taxa=list(pair) if pair else [],
        mrca_leaf_pair=pair,
        is_root_node=is_root,
    )


def _parse_lsd2_value(value: str):
    """把 LSD2 date 文件的取值解析为 (最小年龄 Ma, 最大年龄 Ma)，NA -> None。

    ``-80`` -> (80, 80)（等式）；``-90:-70`` -> (70, 90)（区间）；
    ``-60:NA`` -> (None, 60)（只设上界，年龄 ≤ 60）；``NA:-50`` -> (50, None)
    （只设下界，年龄 ≥ 50）。该单边语义已用 IQ-TREE 2.3.6 / 3.1.3 实测回读
    校准节点年龄确认（见 models/constraints.py 的 LSD2 渲染注释）。
    按语义断言而不是逐字符比对渲染，既不依赖浮点格式化细节，又能抓住
    "不等式被塌成等式"这一 B-16 失效形态。
    """
    sides = [side.strip() for side in value.strip().split(":")]
    if len(sides) == 1:
        point = -float(sides[0])
        return (point, point)

    def one(token):
        if token.upper() == "NA":
            return None
        return -float(token)

    # value = -age，``v1:v2`` 中 v1（更负）对应更老端、v2 对应更年轻端；
    # 因此最小年龄来自右端 v2、最大年龄来自左端 v1（NA 侧无界）。
    lower, upper = one(sides[1]), one(sides[0])
    if lower is not None and upper is not None:
        return (min(lower, upper), max(lower, upper))
    return (lower, upper)


class TestLSD2Method:
    """LSD2Method 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="lsd2_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def lsd2_config(self):
        """创建 LSD2 配置"""
        return LSD2Config(model="GTR+G")

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

    def test_init_creates_work_directory(self, temp_dir, lsd2_config, common_config):
        """测试初始化创建工作目录"""
        method = LSD2Method(lsd2_config, temp_dir, common_config=common_config)

        assert method.work_dir.exists()
        assert method.method_name == "lsd2"

    def test_method_name_property(self, temp_dir, lsd2_config, common_config):
        """测试方法名称属性"""
        method = LSD2Method(lsd2_config, temp_dir, common_config=common_config)
        assert method.method_name == "lsd2"

    def test_common_config_usage(self, temp_dir, lsd2_config, common_config):
        """测试 common_config 使用"""
        method = LSD2Method(lsd2_config, temp_dir, common_config=common_config)
        assert method.common_config.nthreads == 4

    @patch("phylodater.adapters.lsd2_method.ProcessRunner")
    def test_validate_environment_with_iqtree(
        self, mock_runner_class, temp_dir, lsd2_config, common_config
    ):
        """测试环境验证 - IQ-TREE 存在"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = "2.2.0"
        mock_runner_class.return_value = mock_runner

        method = LSD2Method(lsd2_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is True
        mock_runner.check_executable.assert_called()

    @patch("phylodater.adapters.lsd2_method.ProcessRunner")
    def test_validate_environment_without_iqtree(
        self, mock_runner_class, temp_dir, lsd2_config, common_config
    ):
        """测试环境验证 - IQ-TREE 不存在"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = False
        mock_runner_class.return_value = mock_runner

        method = LSD2Method(lsd2_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is False

    def test_model_parameter(self, temp_dir, lsd2_config, common_config):
        """测试模型参数"""
        method = LSD2Method(lsd2_config, temp_dir, common_config=common_config)
        assert method.config.model == "GTR+G"

    def test_nthreads_from_common_config(self, temp_dir, lsd2_config, common_config):
        """测试线程数从 common_config 获取"""
        method = LSD2Method(lsd2_config, temp_dir, common_config=common_config)
        assert method.common_config.nthreads == 4

    def test_work_dir_isolation(self, temp_dir, lsd2_config, common_config):
        """测试工作目录隔离"""
        method1 = LSD2Method(lsd2_config, temp_dir, common_config=common_config)
        method2 = LSD2Method(
            lsd2_config, temp_dir / "other", common_config=common_config
        )

        assert method1.work_dir != method2.work_dir
        assert method1.work_dir.exists()
        assert method2.work_dir.exists()


class TestLSD2Config:
    """LSD2 配置测试"""

    def test_default_config(self):
        """测试默认配置：model 默认不写死，由比对类型自动选择。"""
        config = LSD2Config()
        assert config.model is None

    def test_model_selection(self):
        """测试模型选择"""
        dna_config = LSD2Config(model="GTR+G")
        protein_config = LSD2Config(model="LG+G")

        assert dna_config.model == "GTR+G"
        assert protein_config.model == "LG+G"

    def test_date_ci_validation(self):
        """测试 date_ci 验证"""
        with pytest.raises(ValueError):
            LSD2Config(date_ci=-1)

    def test_clock_sd_validation(self):
        """测试 clock_sd 验证"""
        with pytest.raises(ValueError):
            LSD2Config(clock_sd=0)


class TestLSD2ModelResolution:
    """未指定 model 时按比对类型选择；类型明显不符时直接报错。

    旧默认 ``LG+G`` 是氨基酸模型，而主用例递给 IQ-TREE2 的是核酸比对，IQ-TREE2
    只会报一句难以理解的 ``ERROR: File not found LG``，整个方法直接失败。
    """

    @pytest.fixture
    def temp_dir(self, tmp_path):
        return tmp_path / "lsd2_model"

    def _method(self, temp_dir, model=None):
        temp_dir.mkdir(parents=True, exist_ok=True)
        return LSD2Method(
            LSD2Config(model=model), temp_dir, common_config=CommonConfig()
        )

    def _write_alignment(self, temp_dir, residues, name="aln.fasta"):
        path = temp_dir / name
        path.write_text(f">A\n{residues}\n>B\n{residues}\n", encoding="utf-8")
        return path

    def test_dna_alignment_resolves_to_gtr_plus_g(self, temp_dir):
        method = self._method(temp_dir)
        aln = self._write_alignment(temp_dir, "ACGTACGTACGTACGTACGTACGTACGTACGT")
        assert method._detect_alignment_seq_type(aln) == "nucleotide"
        assert method._resolve_model(aln) == "GTR+G"

    def test_protein_alignment_resolves_to_lg_plus_g(self, temp_dir):
        method = self._method(temp_dir)
        aln = self._write_alignment(temp_dir, "MQIFTQWDELIS PGVRNTIAp".replace(" ", ""))
        assert method._detect_alignment_seq_type(aln) == "protein"
        assert method._resolve_model(aln) == "LG+G"

    def test_explicit_model_is_kept(self, temp_dir):
        method = self._method(temp_dir, model="HKY+I+G")
        aln = self._write_alignment(temp_dir, "ACGTACGTACGTACGTACGTACGTACGTACGT")
        assert method._resolve_model(aln) == "HKY+I+G"

    def test_amino_acid_model_on_dna_alignment_fails_loudly(self, temp_dir):
        """IQ-TREE2 会把它当成文件路径，所以 PhyloDater 必须提前报错。"""
        from phylodater.core.exceptions import ConfigurationError

        method = self._method(temp_dir, model="LG+G")
        aln = self._write_alignment(temp_dir, "ACGTACGTACGTACGTACGTACGTACGTACGT")
        with pytest.raises(ConfigurationError, match="amino-acid model"):
            method._resolve_model(aln)

    def test_missing_alignment_falls_back_to_nucleotide(self, temp_dir):
        """比对不可读时不能静默丢失：按主用例（核酸）处理并记录日志。"""
        method = self._method(temp_dir)
        assert method._detect_alignment_seq_type(temp_dir / "absent.fasta") is None
        assert method._resolve_model(temp_dir / "absent.fasta") == "GTR+G"


class TestLSD2InputGeneration:
    """B-16 / B-17 回归防护：逐行断言 constraints.date 的真实内容。

    此前该路径由适配器内自带的 _age_to_str 处理，把均匀区间、软下界、上限三种
    **不等式**约束一律塌成点值（等式约束）且不发任何降级警告；date 文件还把全部
    叶节点无条件写成 0 Ma，抹掉异时采样的采样日期。
    """

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="lsd2_inputs_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def method(self, temp_dir):
        return LSD2Method(
            LSD2Config(model="GTR+G"), temp_dir, common_config=CommonConfig(nthreads=4)
        )

    @pytest.fixture
    def tree(self):
        return PhylogeneticTree(TEST_NEWICK)

    @staticmethod
    def _date_lines(method):
        text = (method.work_dir / "constraints.date").read_text()
        return [line for line in text.splitlines() if line.strip()]

    @staticmethod
    def _constraint_lines(method):
        """date 文件里的 MRCA/约束行（键含逗号；叶日期行的键是单个分类名）"""
        lines = TestLSD2InputGeneration._date_lines(method)
        return [line for line in lines if "," in line.split("\t")[0]]

    @staticmethod
    def _tip_lines(method):
        """date 文件里的叶节点日期行"""
        lines = TestLSD2InputGeneration._date_lines(method)
        return [line for line in lines if "," not in line.split("\t")[0]]

    @staticmethod
    def _constraint_map(method):
        """{MRCA 叶对键: (最小年龄, 最大年龄)}，NA 记为 None"""
        mapping = {}
        for line in TestLSD2InputGeneration._constraint_lines(method):
            key, _, value = line.partition("\t")
            mapping[key] = _parse_lsd2_value(value)
        return mapping

    # ---- B-16：不等式必须按区间 / 单边 NA 语法写出 ------------------------

    def test_fixed_age_stays_an_exact_point_without_warning(self, method, tree):
        """固定点在 LSD2 上是精确表达：等式取值，且不发约束降级警告"""
        method.prepare_inputs(tree, [_cal("F", FixedAgeConstraint(80.0), ("A", "B"))])

        assert self._constraint_map(method) == {"A,B": (80.0, 80.0)}
        assert not [
            w for w in method._warnings if "cannot represent" in str(w)
        ], method._warnings

    def test_uniform_interval_is_not_collapsed_to_a_point(self, method, tree):
        """UniformAgeConstraint(70,90) 必须写成区间，不得塌成中点 80 的等式"""
        method.prepare_inputs(
            tree, [_cal("U", UniformAgeConstraint(70.0, 90.0), ("A", "B"))]
        )

        ((key, bounds),) = self._constraint_map(method).items()
        assert key == "A,B"
        assert bounds == (70.0, 90.0)
        assert not [
            w for w in method._warnings if "cannot represent" in str(w)
        ], method._warnings

    def test_soft_lower_bound_uses_single_sided_na_syntax(self, method, tree):
        """SoftLowerBoundConstraint(60) -> 只设下界 (60, None)，不是等式 (60, 60)"""
        method.prepare_inputs(
            tree, [_cal("L", SoftLowerBoundConstraint(60.0), ("C", "D"))]
        )

        assert self._constraint_map(method) == {"C,D": (60.0, None)}
        # 语法本身也必须带上缺失界记号
        assert "NA" in self._constraint_lines(method)[0]

    def test_maximum_age_uses_single_sided_na_syntax(self, method, tree):
        """MaximumAgeConstraint(50) -> 只设上界 (None, 50)，不是等式 (50, 50)"""
        method.prepare_inputs(tree, [_cal("M", MaximumAgeConstraint(50.0), ("A", "C"))])

        assert self._constraint_map(method) == {"A,C": (None, 50.0)}
        assert "NA" in self._constraint_lines(method)[0]

    def test_mixed_inequality_constraints_keep_their_semantics(self, method, tree):
        """三种不等式同时进去时语义互不混淆（旧实现会写成 -80/-60/-50 三个等式）"""
        method.prepare_inputs(
            tree,
            [
                _cal("U", UniformAgeConstraint(70.0, 90.0), ("A", "B")),
                _cal("L", SoftLowerBoundConstraint(60.0), ("C", "D")),
                _cal("M", MaximumAgeConstraint(50.0), ("A", "C")),
            ],
        )

        assert self._constraint_map(method) == {
            "A,B": (70.0, 90.0),
            "C,D": (60.0, None),
            "A,C": (None, 50.0),
        }

    def test_gamma_prior_degrades_to_interval_with_warning(self, method, tree):
        """概率先验降级为 95% 硬区间时必须留痕"""
        method.prepare_inputs(
            tree,
            [
                _cal(
                    "G",
                    GammaPriorConstraint(alpha=2.0, beta=0.1, offset=10.0),
                    ("A", "B"),
                )
            ],
        )

        ((key, bounds),) = self._constraint_map(method).items()
        assert key == "A,B"
        assert bounds[0] is not None and bounds[1] is not None
        assert bounds[0] < bounds[1]
        assert any("Gamma" in str(w) for w in method._warnings)

    # ---- B-16：无法精确表达时必须有降级警告 ------------------------------

    def test_degradation_warnings_are_recorded_and_warned(self, temp_dir):
        """软下界 / 上限 / 软区间各发一条 SemanticDegradationWarning（日志 + 结果级）"""
        method = LSD2Method(
            LSD2Config(model="LG+G"), temp_dir, common_config=CommonConfig()
        )
        # 所有叶都带日期 -> 不产生"按现在处理"的警告，只剩三条降级警告
        tree = PhylogeneticTree(
            "((A|-0.1:0.1,B|-0.2:0.1)AB:0.2,(C|-0.3:0.3,D|-0.4:0.3)CD:0.1);"
        )
        cals = [
            _cal("L", SoftLowerBoundConstraint(60.0), ("C", "D")),
            _cal("M", MaximumAgeConstraint(50.0), ("A", "C")),
            _cal("S", SoftBoundsConstraint(100.0, 120.0), ("A", "B")),
        ]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            method.prepare_inputs(tree, cals)

        messages = [str(w) for w in method._warnings]
        assert len(messages) == 3, messages
        assert all(isinstance(w, SemanticDegradationWarning) for w in method._warnings)
        assert any("soft lower bound" in m for m in messages), messages
        assert any("soft upper bound" in m for m in messages), messages
        assert any("soft bounds" in m for m in messages), messages
        # 同时通过 warnings 模块发出，供日志层/报告层捕获
        raised = [
            w for w in caught if issubclass(w.category, SemanticDegradationWarning)
        ]
        assert len(raised) >= 3, raised

    def test_warnings_reach_the_dating_result(self, temp_dir):
        """降级警告写入 DatingResult.warnings（与 pathd8 的既有做法一致）"""
        method = LSD2Method(
            LSD2Config(model="LG+G"), temp_dir, common_config=CommonConfig()
        )
        tree = PhylogeneticTree(TEST_NEWICK)
        method.prepare_inputs(
            tree, [_cal("L", SoftLowerBoundConstraint(60.0), ("C", "D"))]
        )
        timetree = method.work_dir / "lsd2_run.timetree.nex"
        timetree.write_text(
            "#NEXUS\nbegin trees;\n"
            '  tree tree1 = [&R] ((A:0.1[&date="-0.05"],B:0.1[&date="-0.05"])'
            '[&date="-80.0"]:0.2,(C:0.3[&date="-0.05"],D:0.3[&date="-0.05"])'
            '[&date="-60.0"]:0.1)[&date="-120.0"];\nend;\n'
        )
        method._input_files["timetree_nex"] = timetree

        result = method.parse_results()
        assert any("soft lower bound" in str(w) for w in result.warnings)
        assert (
            result.node_ages
        ), "node ages should be decoded from [&date=...] annotations"

    # ---- B-17：叶节点采样日期 ---------------------------------------------

    def test_tip_without_dates_written_as_zero_but_warns(self, method, tree):
        """全树无日期时仍写 0（保证 date 文件可用），但必须告警说明按"现在"处理"""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            method.prepare_inputs(
                tree, [_cal("U", UniformAgeConstraint(70.0, 90.0), ("A", "B"))]
            )

        tip_lines = self._tip_lines(method)
        assert tip_lines == ["A\t0", "B\t0", "C\t0", "D\t0"]
        assert any("contemporaneous" in str(w) for w in method._warnings)
        assert any("contemporaneous" in str(w.message) for w in caught)

    def test_heterochronous_tip_dates_are_preserved(self, temp_dir):
        """分类名里带的采样日期（A|-0.045 / C[&&DATE=-0.12]）必须原样写入 date 文件"""
        method = LSD2Method(
            LSD2Config(model="LG+G"), temp_dir, common_config=CommonConfig()
        )
        tree = PhylogeneticTree(
            "((A|-0.045:0.1,B:0.1)AB:0.2,(C[&&DATE=-0.12]:0.3,D:0.3)CD:0.1);"
        )

        method.prepare_inputs(
            tree, [_cal("U", UniformAgeConstraint(70.0, 90.0), ("A", "B"))]
        )

        tip_lines = {
            line.split("\t")[0]: line.split("\t")[1] for line in self._tip_lines(method)
        }
        assert tip_lines["A|-0.045"] == "-0.045"
        assert tip_lines["C"] == "-0.12"
        # 没有日期的叶仍为 0，并被告警点名
        assert tip_lines["B"] == "0" and tip_lines["D"] == "0"
        assert any(
            "2/4 tip(s) carry no sampling date" in str(w) for w in method._warnings
        )

    def test_calendar_tip_dates_fail_loudly_instead_of_unit_mismatch(self, temp_dir):
        """日历式日期与本适配器输出的 Ma 轴节点约束不可混用：显式报错而不是猜"""
        method = LSD2Method(
            LSD2Config(model="LG+G"), temp_dir, common_config=CommonConfig()
        )
        tree = PhylogeneticTree(
            "((hCoV/Wuhan|2019-12-31:0.1,B:0.1)AB:0.2,(C:0.3,D:0.3)CD:0.1);"
        )

        with pytest.raises(CalibrationError, match="CALENDAR sampling dates"):
            method.prepare_inputs(
                tree, [_cal("U", UniformAgeConstraint(70.0, 90.0), ("A", "B"))]
            )

    def test_number_in_taxon_name_is_not_mistaken_for_a_date(self, method):
        """登录号尾随的裸整数不当作日期（避免把名字里的数字写成采样时间）"""
        tree = PhylogeneticTree("((GCA_0001:0.1,B:0.1)AB:0.2,(C:0.3,D:0.3)CD:0.1);")
        numeric, calendar = method._extract_tip_sampling_dates(tree)

        assert numeric == {} and calendar == {}

    def test_tree_file_is_written_without_internal_labels(self, method, tree):
        """输入树落盘且剥掉内部标签（避免约束按名字挂到错误的旧标签上）"""
        method.prepare_inputs(
            tree, [_cal("U", UniformAgeConstraint(70.0, 90.0), ("A", "B"))]
        )

        written = (method.work_dir / "input_tree.nwk").read_text()
        assert "AB" not in written and "CD" not in written
        assert written.strip().endswith(";")

    def test_root_calibration_uses_a_leaf_pair_across_the_root(self, method, tree):
        """根校准用跨根两侧的叶对表达（dendropy 定位，不依赖 ete3）"""
        method.prepare_inputs(
            tree,
            [
                _cal("U", UniformAgeConstraint(70.0, 90.0), ("A", "B")),
                _cal("Root", FixedAgeConstraint(200.0), is_root=True),
            ],
        )

        assert self._constraint_map(method) == {
            "A,C": (200.0, 200.0),
            "A,B": (70.0, 90.0),
        }

    # ---- 丢点必须留痕 / 零产出必须失败 -----------------------------------

    def test_calibration_without_mrca_pair_warns_before_skipping(self, method, tree):
        """缺 mrca_leaf_pair 的校准被跳过时必须有 warning（此前静默丢弃）"""
        method.logger = Mock()
        method.prepare_inputs(
            tree,
            [
                _cal("Kept", FixedAgeConstraint(80.0), ("A", "B")),
                _cal("Dropped", FixedAgeConstraint(30.0), None),
            ],
        )

        assert len(self._constraint_lines(method)) == 1
        assert "Dropped" not in (method.work_dir / "constraints.date").read_text()
        warnings_logged = [
            str(call.args[0]) for call in method.logger.warning.call_args_list
        ]
        assert any(
            "Dropped" in m and "no mrca_leaf_pair" in m for m in warnings_logged
        ), warnings_logged

    def test_missing_age_constraint_raises_instead_of_silent_drop(self, method, tree):
        """没有 age_constraint 的校准点显式失败，而不是少写一行"""
        with pytest.raises(CalibrationError, match="no age_constraint"):
            method.prepare_inputs(tree, [_cal("Empty", None, ("A", "B"))])

    def test_zero_internal_calibrations_still_raises(self, method, tree):
        """只有尾日期、无内部节点校准：保留原有的响亮失败"""
        with pytest.raises(CalibrationError, match="at least one internal node"):
            method.prepare_inputs(
                tree, [_cal("Root", FixedAgeConstraint(90.0), None, True)]
            )

    def test_repeated_prepare_inputs_does_not_accumulate_warnings(self, method, tree):
        """每次 prepare_inputs 重置警告列表，避免同一降级被记录多次"""
        cals = [_cal("L", SoftLowerBoundConstraint(60.0), ("C", "D"))]
        method.prepare_inputs(tree, cals)
        first = len(method._warnings)
        method.prepare_inputs(tree, cals)

        assert len(method._warnings) == first
