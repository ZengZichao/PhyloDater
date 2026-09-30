"""
TreePLMethod 适配器单元测试

测试 TreePL 定年方法的各项功能
"""

import shutil
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.adapters.treepl_method import TreePLMethod
from phylodater.core.exceptions import CalibrationError
from phylodater.infrastructure.configuration import CommonConfig, TreePLConfig
from phylodater.infrastructure.safe_io import safe_writer
from phylodater.models import (
    CalibrationPoint,
    FixedAgeConstraint,
    MaximumAgeConstraint,
    PhylogeneticTree,
    SoftBoundsConstraint,
    UniformAgeConstraint,
)


class TestTreePLMethod:
    """TreePLMethod 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="treepl_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def treepl_config(self):
        """创建 TreePL 配置"""
        return TreePLConfig(
            initial_smooth=100.0,
            cv=True,
            randomcv=True,
            run_prime=True,
        )

    @pytest.fixture
    def common_config(self):
        """创建通用配置"""
        return CommonConfig(nthreads=4, seed=42)

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
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

    def test_init_creates_work_directory(self, temp_dir, treepl_config, common_config):
        """测试初始化创建工作目录"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        assert method.work_dir.exists()
        assert method.method_name == "treepl"

    def test_method_name_property(self, temp_dir, treepl_config, common_config):
        """测试方法名称属性"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)
        assert method.method_name == "treepl"

    def test_init_sets_random_seed(self, temp_dir, treepl_config):
        """测试初始化设置随机种子"""
        common_config = CommonConfig(seed=None)
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        # 种子应该被自动设置
        assert method.common_config.seed is not None
        assert isinstance(method.common_config.seed, int)

    def test_init_preserves_user_seed(self, temp_dir, treepl_config, common_config):
        """测试初始化保留用户设置的种子"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        assert method.common_config.seed == 42

    def test_actual_multiplier_default(self, temp_dir, treepl_config, common_config):
        """测试实际乘数默认值"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        assert method._actual_multiplier == 1.0
        assert method._tree_modified is False

    @patch("phylodater.adapters.treepl_method.ProcessRunner")
    def test_validate_environment(
        self, mock_runner_class, temp_dir, treepl_config, common_config
    ):
        """测试环境验证"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = "1.0"
        mock_runner_class.return_value = mock_runner

        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is True

    @patch("phylodater.adapters.treepl_method.ProcessRunner")
    def test_validate_environment_not_found(
        self, mock_runner_class, temp_dir, treepl_config, common_config
    ):
        """测试环境验证 - 未找到"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = False
        mock_runner_class.return_value = mock_runner

        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is False

    def test_write_common_params_nthreads(self, temp_dir, treepl_config, common_config):
        """测试通用参数写入线程数"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        # 创建临时文件测试参数写入
        test_file = temp_dir / "test.conf"
        with safe_writer(test_file) as f:
            method._write_common_params(f, Path("tree.nwk"))

        with open(test_file, "r") as f:
            content = f.read()

        assert "nthreads = 4" in content

    def test_write_common_params_seed(self, temp_dir, treepl_config, common_config):
        """测试通用参数写入种子"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        test_file = temp_dir / "test.conf"
        with safe_writer(test_file) as f:
            method._write_common_params(f, Path("tree.nwk"))

        with open(test_file, "r") as f:
            content = f.read()

        assert "seed = 42" in content

    def test_write_common_params_verbose(self, temp_dir, treepl_config):
        """测试通用参数写入详细输出"""
        common_config = CommonConfig(verbose=True)
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        test_file = temp_dir / "test.conf"
        with safe_writer(test_file) as f:
            method._write_common_params(f, Path("tree.nwk"))

        with open(test_file, "r") as f:
            content = f.read()

        assert "verbose" in content

    def test_write_treepl_constraint_fixed(
        self, temp_dir, treepl_config, common_config
    ):
        """测试固定约束写入"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        cal = Mock(spec=CalibrationPoint)
        cal.name = "TestCal"
        cal.age_constraint = FixedAgeConstraint(fixed_age=10.0)

        test_file = temp_dir / "test.conf"
        with safe_writer(test_file) as f:
            method._write_treepl_constraint(f, cal)

        with open(test_file, "r") as f:
            content = f.read()

        assert "min = TestCal 10.0" in content
        assert "max = TestCal 10.0" in content

    def test_write_treepl_constraint_uniform(
        self, temp_dir, treepl_config, common_config
    ):
        """测试均匀分布约束写入"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        cal = Mock(spec=CalibrationPoint)
        cal.name = "TestCal"
        cal.age_constraint = UniformAgeConstraint(min_age=6.0, max_age=8.0)

        test_file = temp_dir / "test.conf"
        with safe_writer(test_file) as f:
            method._write_treepl_constraint(f, cal)

        with open(test_file, "r") as f:
            content = f.read()

        assert "min = TestCal 6.0" in content
        assert "max = TestCal 8.0" in content

    def test_write_treepl_constraint_soft_bounds(
        self, temp_dir, treepl_config, common_config
    ):
        """测试软边界约束写入（降级为均匀分布）"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        cal = Mock(spec=CalibrationPoint)
        cal.name = "TestCal"
        cal.age_constraint = SoftBoundsConstraint(
            min_age=5.0, max_age=10.0, tail_prob=0.025
        )

        test_file = temp_dir / "test.conf"
        with safe_writer(test_file) as f:
            method._write_treepl_constraint(f, cal)

        with open(test_file, "r") as f:
            content = f.read()

        assert "min = TestCal 5.0" in content
        assert "max = TestCal 10.0" in content

    def test_parse_prime_output_success(self, temp_dir, treepl_config, common_config):
        """测试 Prime 输出解析成功"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        output = """
        PLACE THE LINES BELOW IN THE CONFIGURATION FILE
        opt = 5
        optad = 3
        """

        result = method._parse_prime_output(output)

        assert result.success is True
        assert result.opt == 5
        assert result.optad == 3

    def test_parse_prime_output_failure(self, temp_dir, treepl_config, common_config):
        """测试 Prime 输出解析失败"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        output = "Some random output without opt values"

        result = method._parse_prime_output(output)

        assert result.success is False

    def test_parse_cv_file(self, temp_dir, treepl_config, common_config):
        """测试 CV 文件解析"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        content = """
        smoothing    chi-square
        1000.0       123.45
        100.0        89.12
        10.0         95.67
        """

        result = method._parse_cv_file(content)

        assert result == 100.0  # chi-square 最小的平滑值

    def test_parse_cv_file_empty(self, temp_dir, treepl_config, common_config):
        """测试 CV 文件解析 - 空文件"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        content = ""

        result = method._parse_cv_file(content)

        assert result is None

    def test_sanitize_name(self, temp_dir, treepl_config, common_config):
        """测试名称安全化"""
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        assert method._sanitize_name("Test Name") == "Test_Name"
        assert method._sanitize_name("Test@Name#123") == "Test_Name_123"
        assert method._sanitize_name("A" * 40) == "A" * 30  # 截断到 30 字符

    def test_actual_multiplier_used_in_adjustment(
        self, temp_dir, treepl_config, common_config
    ):
        """测试回退调整确实按 _actual_multiplier 走，并且失败时会说话

        原先这个测试写着"不实际调用，只验证方法存在"然后 `pass`，
        无论实现如何都会通过（审阅报告 §四 批评的自我实现型测试）。
        """
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)
        method.logger = Mock()

        method._tree_modified = True
        method._actual_multiplier = 1000.0
        method._calibrations = []

        newick = "((human:0.001,chimp:0.001):0.002,(mouse:0.003,rat:0.003):0.001);"
        out = method._adjust_dated_tree_lengths(newick)

        # ete3 可用时应按比例除回；不可用（本工作区 Python 3.14 上 ete3 装不上）
        # 时必须**说出来**并原样返回，而不是静默给出未回退的树。
        try:
            from ete3 import Tree as _Ete3Tree  # 仅探测 ete3 在本解释器上能否导入

            ete3_available = _Ete3Tree is not None
        except Exception:
            ete3_available = False

        if ete3_available:
            assert out != newick
            # ete3/本实现可能输出 1e-06 或 0.000001，两种都要接受
            assert ("0.000001" in out) or ("1e-06" in out) or ("1e-06" in out.lower())
        else:
            assert out == newick
            assert method.logger.warning.called

    def test_check_tiny_branches_disabled(self, temp_dir, treepl_config, common_config):
        """测试禁用 tiny branch 检查"""
        treepl_config.check_tiny_branches = False
        method = TreePLMethod(treepl_config, temp_dir, common_config=common_config)

        mock_tree = Mock(spec=PhylogeneticTree)
        result = method._check_and_fix_tiny_branches(mock_tree)

        # 应该返回原始树
        assert result is mock_tree
        assert method._tree_modified is False


class TestTreePLConfig:
    """TreePL 配置测试"""

    def test_default_config(self):
        """测试默认配置"""
        config = TreePLConfig()

        assert config.initial_smooth == 100.0
        assert config.cv is True
        assert config.randomcv is True
        assert config.run_prime is True

    def test_custom_config(self):
        """测试自定义配置"""
        config = TreePLConfig(
            initial_smooth=50.0,
            cvstart=500.0,
            cvstop=0.01,
        )

        assert config.initial_smooth == 50.0
        assert config.cvstart == 500.0
        assert config.cvstop == 0.01

    def test_validation_cvstop_must_be_less_than_cvstart(self):
        """测试验证 cvstop 必须小于 cvstart"""
        with pytest.raises(ValueError, match="cvstop.*must be.*cvstart"):
            TreePLConfig(cvstart=100.0, cvstop=200.0)

    def test_validation_initial_smooth_must_be_positive(self):
        """测试验证 initial_smooth 必须为正"""
        with pytest.raises(ValueError, match="initial_smooth must be > 0"):
            TreePLConfig(initial_smooth=-1.0)

    def test_validation_cviter_must_be_positive(self):
        """测试验证 cviter 必须为正"""
        with pytest.raises(ValueError, match="cviter must be >= 1"):
            TreePLConfig(cviter=0)


TREEPL_NEWICK = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
TREEPL_TIPS = ["human", "chimp", "mouse", "rat"]


def _treepl_tree():
    tree = Mock(spec=PhylogeneticTree)
    tree.newick = TREEPL_NEWICK
    tree.tip_names = list(TREEPL_TIPS)
    tree.num_tips = len(TREEPL_TIPS)
    tree.without_internal_labels = Mock(return_value=tree)
    # _build_name_mapping 用 ACCESSION 模式，总会产出一套短名，
    # 于是 prepare_inputs 会走 with_renamed_leaves；这里让它返回同一个
    # mock，好让 tip_names 仍然可迭代。
    tree.with_renamed_leaves = Mock(return_value=tree)
    tree.write = Mock()
    return tree


def _treepl_cal(name, constraint, pair=("human", "chimp"), root=False, taxa=None):
    cal = Mock(spec=CalibrationPoint)
    cal.name = name
    cal.age_constraint = constraint
    cal.mrca_leaf_pair = pair
    cal.is_root_node = root
    cal.resolved_taxa = taxa if taxa is not None else (list(pair) if pair else [])
    return cal


class TestTreePLWrittenConfigGuard:
    """treePL 的前置守卫要看**写出去的配置文件**，不看内存里的校准列表

    treePL 在多种错误路径上打印信息后 exit(0)（审阅报告 §六 末"共性风险"），
    所以"配置里一条约束都没有"必须在本机就拦下，而不是等上游静默跑完。
    """

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="treepl_guard_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def method(self, temp_dir):
        return TreePLMethod(
            TreePLConfig(), temp_dir, common_config=CommonConfig(nthreads=2, seed=7)
        )

    def test_config_with_no_constraint_is_refused(self, method):
        """校准既无 mrca_leaf_pair 也无 resolved_taxa → 配置里零约束 → 拒绝"""
        with pytest.raises(CalibrationError, match="no usable calibration"):
            method.prepare_inputs(
                _treepl_tree(),
                [
                    _treepl_cal(
                        "Homeless",
                        UniformAgeConstraint(min_age=6.0, max_age=8.0),
                        pair=None,
                        taxa=[],
                    )
                ],
            )

    def test_unresolvable_calibration_is_named(self, method):
        """掉下去的那条校准必须点名，并且不能拖垮整次运行"""
        method.prepare_inputs(
            _treepl_tree(),
            [
                _treepl_cal("Primates", UniformAgeConstraint(min_age=6.0, max_age=8.0)),
                _treepl_cal(
                    "Homeless",
                    UniformAgeConstraint(min_age=6.0, max_age=8.0),
                    pair=None,
                    taxa=[],
                ),
            ],
        )
        prime = method._input_files["prime_config"].read_text()
        assert "mrca = Primates" in prime
        assert "min = Primates 6.0" in prime
        assert "Homeless" not in prime  # 确实没能落盘

        recorded = " ".join(str(w) for w in method._warnings)
        assert "Calibration 'Homeless' was NOT written" in recorded
        assert "2 calibration(s) supplied, 1 of them reached" in recorded

    def test_root_calibration_still_gets_a_full_tip_mrca(self, method):
        """回归护栏（审阅报告 §六.15(b)）：根校准要用全部叶节点定义 mrca"""
        method.prepare_inputs(
            _treepl_tree(),
            [
                _treepl_cal(
                    "Root",
                    MaximumAgeConstraint(max_age=500.0),
                    pair=None,
                    root=True,
                    taxa=[],
                )
            ],
        )
        prime = method._input_files["prime_config"].read_text()
        mrca_lines = [
            line for line in prime.splitlines() if line.startswith("mrca = Root ")
        ]
        assert len(mrca_lines) == 1
        # "mrca = Root <tip1> ... <tipN>"：全部 4 片叶子（经短名映射后
        # 仍与树里的 tip 名一致）都参与定义根节点
        assert len(mrca_lines[0].split()) == 3 + len(TREEPL_TIPS)
        assert "max = Root 500.0" in prime

    def test_warning_is_not_duplicated_across_the_three_stages(self, method):
        """prime/CV/final 都写一遍约束，同一问题只报一次"""
        method.prepare_inputs(
            _treepl_tree(),
            [
                _treepl_cal("Primates", UniformAgeConstraint(min_age=6.0, max_age=8.0)),
                _treepl_cal(
                    "Homeless",
                    UniformAgeConstraint(min_age=6.0, max_age=8.0),
                    pair=None,
                    taxa=[],
                ),
            ],
        )
        method._generate_final_config(
            method._input_files["tree"], method.config.initial_smooth
        )
        msgs = [str(w) for w in method._warnings if "was NOT written" in str(w)]
        assert len(msgs) == 1


class TestTreePLNoSilentSwallow:
    """C-8：本文件里过去有 5 处 `except Exception: pass`，吞掉的问题都得留痕"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="treepl_swallow_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def method(self, temp_dir):
        return TreePLMethod(
            TreePLConfig(), temp_dir, common_config=CommonConfig(nthreads=2)
        )

    def test_version_probe_failure_is_logged_not_swallowed(self, method):
        runner = Mock()
        runner.get_version.side_effect = OSError("probe blew up")
        method.logger = Mock()

        assert method._resolve_version_for_provenance(runner) == "unknown"
        assert method.logger.debug.called

    def test_version_probe_returns_the_real_string(self, method):
        runner = Mock()
        runner.get_version.return_value = "treePL version 1.0"
        assert method._resolve_version_for_provenance(runner) == "treePL version 1.0"

    def test_binary_version_scan_failure_is_logged(self, method):
        method.logger = Mock()
        # 不存在的路径：shutil.which 返回 None → 直接返回 None，不该抛错
        assert method._extract_treepl_version("definitely-not-a-real-binary") is None

    def test_internal_node_annotation_failure_is_reported(self, method):
        """标注失败时输出树只是少一个标签，但必须说出来（过去是 pass）"""
        method.logger = Mock()
        method._calibrations = [
            _treepl_cal(
                "Ghost", FixedAgeConstraint(fixed_age=10.0), pair=("nope", "nada")
            )
        ]
        newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
        out = method._annotate_internal_nodes(newick, None)

        assert out  # 仍然返回一棵树（不中止运行）
        assert method.logger.warning.called
