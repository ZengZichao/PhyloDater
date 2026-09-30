"""
WLogDateMethod 适配器单元测试

测试 wLogDate 定年方法的各项功能：
- CLI 参数构建（-i, -t, -o 正确性）
- 向后时间取反逻辑
- 输出文件读取（非 stdout）
- 约束降级警告
- 配置访问模式
"""

import shutil
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.adapters.wlogdate_method import WLogDateMethod
from phylodater.core.exceptions import CalibrationError
from phylodater.infrastructure.configuration import ToolConfig
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


class TestWLogDateMethod:
    """WLogDateMethod 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="wlogdate_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def tool_config(self):
        """创建 ToolConfig（使用默认 WLogDateConfig）"""
        return ToolConfig()

    @pytest.fixture
    def tool_config_backward(self):
        """创建启用向后时间的 ToolConfig"""
        config = ToolConfig()
        config.wlogdate.backward_time = True
        config.wlogdate.num_replicates = 5
        config.wlogdate.sequence_length = 2000
        config.common.seed = 42
        return config

    @pytest.fixture
    def tool_config_forward(self):
        """创建前向时间的 ToolConfig"""
        config = ToolConfig()
        config.wlogdate.backward_time = False
        return config

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.num_tips = 4
        tree.write = Mock()
        tree.without_internal_labels = Mock(return_value=tree)
        return tree

    @pytest.fixture
    def mock_calibrations(self):
        """创建模拟校准点列表"""
        cal1 = Mock(spec=CalibrationPoint)
        cal1.name = "Primates"
        cal1.resolved_taxa = ["human", "chimp"]
        cal1.mrca_leaf_pair = ("human", "chimp")
        cal1.is_root_node = False
        cal1.age_constraint = UniformAgeConstraint(min_age=6.0, max_age=8.0)

        cal2 = Mock(spec=CalibrationPoint)
        cal2.name = "Root"
        cal2.resolved_taxa = ["human", "chimp", "mouse", "rat"]
        cal2.mrca_leaf_pair = ("human", "rat")
        cal2.is_root_node = True
        cal2.age_constraint = FixedAgeConstraint(fixed_age=90.0)

        return [cal1, cal2]

    # ---- 初始化测试 ----

    def test_init_creates_work_directory(self, temp_dir, tool_config):
        """测试初始化创建工作目录"""
        method = WLogDateMethod(tool_config, temp_dir)

        assert method.work_dir.exists()
        assert method.method_name == "wlogdate"

    def test_config_access_pattern(self, temp_dir, tool_config):
        """测试配置通过 self.config.wlogdate 访问"""
        method = WLogDateMethod(tool_config, temp_dir)

        # self.config 是 ToolConfig
        assert method.config is tool_config
        # self._wld_config 是 WLogDateConfig
        assert method._wld_config is tool_config.wlogdate
        assert method._wld_config.wlogdate_bin == "launch_wLogDate.py"
        assert method._wld_config.backward_time is True

    def test_method_name_property(self, temp_dir, tool_config):
        """测试方法名称属性"""
        method = WLogDateMethod(tool_config, temp_dir)
        assert method.method_name == "wlogdate"

    # ---- 环境验证测试 ----

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_validate_environment_success(
        self, mock_runner_class, temp_dir, tool_config
    ):
        """测试环境验证成功"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = "wLogDate 1.0.4"
        mock_runner_class.return_value = mock_runner

        method = WLogDateMethod(tool_config, temp_dir)
        result = method.validate_environment()

        assert result is True
        mock_runner.check_executable.assert_called_once_with("launch_wLogDate.py")

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_validate_environment_not_found(
        self, mock_runner_class, temp_dir, tool_config
    ):
        """测试环境验证失败"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = False
        mock_runner_class.return_value = mock_runner

        method = WLogDateMethod(tool_config, temp_dir)
        result = method.validate_environment()

        assert result is False

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_validate_environment_custom_path(
        self, mock_runner_class, temp_dir, tool_config
    ):
        """测试通过 software_paths 自定义路径"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = None
        mock_runner_class.return_value = mock_runner

        method = WLogDateMethod(tool_config, temp_dir)
        # 模拟 software_paths 中的自定义路径
        method._software_paths.wlogdate_bin = "/custom/path/wlogdate"
        result = method.validate_environment()

        assert result is True
        mock_runner.check_executable.assert_called_once_with("/custom/path/wlogdate")

    # ---- 向后时间取反测试 ----

    def test_prepare_age_backward_time(self, temp_dir, tool_config_backward):
        """测试向后时间模式下年龄不再由 PhyloDater 取反（wLogDate CLI 自行处理）"""
        method = WLogDateMethod(tool_config_backward, temp_dir)

        # 正向 Ma 值应保持不变；wLogDate 在 -b 模式下内部取反
        assert method._prepare_age(65.0) == 65.0
        assert method._prepare_age(10.0) == 10.0
        assert method._prepare_age(0.5) == 0.5

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_root_time_unit_consistency(
        self, mock_runner_class, temp_dir, mock_tree, mock_calibrations
    ):
        """根校准推导的 -r 值应不随 time_unit（Ma/Ga）改变而漂移。

        校准年龄以 Ma 存储；prepare_inputs 会按 time_unit 归一化后再存入
        root_time，execute 阶段的 _prepare_age（即 _convert_time_unit）只做
        一次 Ga→Ma 转换。无论配置为 Ma 还是 Ga，传入 wLogDate 的 -r 内部值
        都应相同（避免旧版 Ga 下把已是 Ma 的校准年龄再 ×1000 造成 1000× 偏差）。
        """
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner

        newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"

        def run_and_get_r(cfg):
            method = WLogDateMethod(cfg, temp_dir)
            method.prepare_inputs(mock_tree, mock_calibrations)
            output_file = method.work_dir / "wlogdate_output.tre"
            output_file.write_text(newick)
            method.execute()
            cmd = mock_runner.run.call_args[0][0]
            return float(cmd[cmd.index("-r") + 1])

        # Ma 配置（默认，WLogDateConfig 无 time_unit 字段 → 视为 Ma）
        cfg_ma = ToolConfig()
        cfg_ma.wlogdate.backward_time = True
        r_ma = run_and_get_r(cfg_ma)

        # Ga 配置（monkeypatch time_unit）
        cfg_ga = ToolConfig()
        cfg_ga.wlogdate.backward_time = True
        method_ga = WLogDateMethod(cfg_ga, temp_dir)
        method_ga._wld_config.time_unit = "Ga"
        method_ga.prepare_inputs(mock_tree, mock_calibrations)
        ga_out = method_ga.work_dir / "wlogdate_output.tre"
        ga_out.write_text(newick)
        method_ga.execute()
        cmd_ga = mock_runner.run.call_args[0][0]
        r_ga = float(cmd_ga[cmd_ga.index("-r") + 1])

        # 根校准为 FixedAgeConstraint(fixed_age=90.0)
        assert r_ma == pytest.approx(90.0)
        assert r_ga == pytest.approx(90.0)
        assert r_ma == pytest.approx(r_ga)

    def test_prepare_age_forward_time(self, temp_dir, tool_config_forward):
        """测试前向时间模式下年龄不变"""
        method = WLogDateMethod(tool_config_forward, temp_dir)

        assert method._prepare_age(65.0) == 65.0
        assert method._prepare_age(10.0) == 10.0

    def test_convert_time_unit_ga(self, temp_dir, tool_config):
        """测试 Ga -> Ma 单位转换"""
        method = WLogDateMethod(tool_config, temp_dir)
        method._wld_config.time_unit = "Ga"

        assert method._convert_time_unit(1.0) == 1000.0
        assert method._convert_time_unit(0.065) == 65.0

    def test_convert_time_unit_ma(self, temp_dir, tool_config):
        """测试 Ma 单位不变"""
        method = WLogDateMethod(tool_config, temp_dir)

        assert method._convert_time_unit(65.0) == 65.0

    # ---- 约束降级测试 ----

    def test_extract_age_fixed(self, temp_dir, tool_config):
        """测试 FixedAgeConstraint 提取"""
        method = WLogDateMethod(tool_config, temp_dir)
        constraint = FixedAgeConstraint(fixed_age=65.0)
        assert method._extract_age(constraint, "test") == 65.0

    def test_extract_age_uniform_degrades_to_midpoint(self, temp_dir, tool_config):
        """测试 UniformAgeConstraint 降级为中点"""
        method = WLogDateMethod(tool_config, temp_dir)
        constraint = UniformAgeConstraint(min_age=10.0, max_age=20.0)
        assert method._extract_age(constraint, "test") == 15.0

    def test_extract_age_soft_lower_bound(self, temp_dir, tool_config):
        """测试 SoftLowerBoundConstraint 降级为最小值"""
        method = WLogDateMethod(tool_config, temp_dir)
        constraint = SoftLowerBoundConstraint(min_age=50.0)
        assert method._extract_age(constraint, "test") == 50.0

    def test_extract_age_maximum(self, temp_dir, tool_config):
        """测试 MaximumAgeConstraint 与 UniformAgeConstraint 统一降级为区间中点"""
        method = WLogDateMethod(tool_config, temp_dir)
        constraint = MaximumAgeConstraint(max_age=100.0)
        assert method._extract_age(constraint, "test") == 50.0

    def test_extract_age_soft_bounds(self, temp_dir, tool_config):
        """测试 SoftBoundsConstraint 降级为中点"""
        method = WLogDateMethod(tool_config, temp_dir)
        constraint = SoftBoundsConstraint(min_age=10.0, max_age=20.0)
        assert method._extract_age(constraint, "test") == 15.0

    def test_extract_age_gamma(self, temp_dir, tool_config):
        """测试 GammaPriorConstraint 降级为均值"""
        method = WLogDateMethod(tool_config, temp_dir)
        constraint = GammaPriorConstraint(alpha=2.0, beta=0.1, offset=10.0)
        # mean = alpha/beta + offset = 2.0/0.1 + 10.0 = 30.0
        assert method._extract_age(constraint, "test") == 30.0

    # ---- 输出解析测试 ----

    # ---- 输出解析测试 ----

    def test_parse_tree_ages_matches_labels(self, temp_dir, tool_config):
        """测试按节点标签名精确匹配校准点年龄"""
        method = WLogDateMethod(tool_config, temp_dir)

        cal_root = Mock(spec=CalibrationPoint)
        cal_root.name = "Root"
        cal_root.mrca_leaf_pair = ("human", "rat")
        cal_root.is_root_node = True

        cal_primates = Mock(spec=CalibrationPoint)
        cal_primates.name = "Primates"
        cal_primates.mrca_leaf_pair = ("human", "chimp")
        cal_primates.is_root_node = False

        method._calibrations = [cal_root, cal_primates]

        newick = (
            "((human:0.1,chimp:0.1)Primates[t=30.0,mu=0.05]:0.2,"
            "(mouse:0.3,rat:0.3):0.15)Root[t=65.0,mu=0.03];"
        )
        node_ages = method._parse_tree_ages(newick)

        assert "Root" in node_ages
        assert node_ages["Root"].mean_age == 65.0
        assert "Primates" in node_ages
        assert node_ages["Primates"].mean_age == 30.0

    def test_parse_tree_ages_extracts_unnamed_internal_nodes(
        self, temp_dir, tool_config
    ):
        """测试提取未命名内部节点"""
        method = WLogDateMethod(tool_config, temp_dir)

        cal_root = Mock(spec=CalibrationPoint)
        cal_root.name = "Root"
        cal_root.mrca_leaf_pair = ("human", "rat")

        cal_primates = Mock(spec=CalibrationPoint)
        cal_primates.name = "Primates"
        cal_primates.mrca_leaf_pair = ("human", "chimp")

        method._calibrations = [cal_root, cal_primates]

        newick = (
            "((human:0.1,chimp:0.1)Primates[t=30.0,mu=0.05]:0.2,"
            "(mouse:0.3,rat:0.3)[t=50.0,mu=0.04]:0.15)Root[t=65.0,mu=0.03];"
        )
        node_ages = method._parse_tree_ages(newick)

        assert len(node_ages) == 3
        assert node_ages["Root"].mean_age == 65.0
        assert node_ages["Primates"].mean_age == 30.0
        assert "internal_node_1" in node_ages
        assert node_ages["internal_node_1"].mean_age == 50.0

    def test_parse_tree_ages_backward_time_positive(
        self, temp_dir, tool_config_backward
    ):
        """测试 wLogDate 向后时间输出为正，解析时不再取反"""
        method = WLogDateMethod(tool_config_backward, temp_dir)

        cal_root = Mock(spec=CalibrationPoint)
        cal_root.name = "Root"
        cal_root.mrca_leaf_pair = ("human", "rat")
        cal_root.is_root_node = True
        method._calibrations = [cal_root]

        newick = (
            "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.15)Root[t=65.0,mu=0.03];"
        )
        node_ages = method._parse_tree_ages(newick)

        assert node_ages["Root"].mean_age == 65.0

    def test_regex_matches_node_annotation(self, temp_dir, tool_config):
        """测试正则能正确识别带标签、无标签和根节点注释"""
        import re

        pattern = r"\)([^\[\]()]*?)\[t=(-?[\d.]+(?:[eE][-+]?\d+)?)(?:,mu=[\d.]+(?:[eE][-+]?\d+)?)?\](?=[;:])"
        newick = (
            "((human:0.1,chimp:0.1)Primates[t=30.0,mu=0.05]:0.2,"
            "(mouse:0.3,rat:0.3)[t=50.0,mu=0.04]:0.15)Root[t=65.0,mu=0.03];"
        )
        matches = [
            (m.group(1).strip(), float(m.group(2)))
            for m in re.finditer(pattern, newick)
        ]
        assert matches == [
            ("Primates", 30.0),
            ("", 50.0),
            ("Root", 65.0),
        ]

    # ---- CLI 参数构建测试 ----

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_execute_cli_flags(
        self,
        mock_runner_class,
        temp_dir,
        tool_config_backward,
        mock_tree,
        mock_calibrations,
    ):
        """测试 CLI 参数包含正确的 -i, -t, -o 标志"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner

        method = WLogDateMethod(tool_config_backward, temp_dir)
        method.prepare_inputs(mock_tree, mock_calibrations)

        # 创建模拟输出文件
        output_file = method.work_dir / "wlogdate_output.tre"
        output_file.write_text("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")

        method.execute()

        # 获取实际调用的命令
        call_args = mock_runner.run.call_args
        cmd = call_args[0][0]

        # 验证关键参数
        assert "-i" in cmd
        assert "-t" in cmd
        assert "-o" in cmd
        assert "-b" in cmd  # backward_time=True
        assert "-p" in cmd
        assert "5" in cmd  # num_replicates=5

        # 验证 -i 在 -t 之前（输入树在前）
        i_idx = cmd.index("-i")
        t_idx = cmd.index("-t")
        o_idx = cmd.index("-o")
        assert i_idx < t_idx < o_idx

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_execute_algorithm_params(
        self,
        mock_runner_class,
        temp_dir,
        tool_config_backward,
        mock_tree,
        mock_calibrations,
    ):
        """测试算法参数传递"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner

        method = WLogDateMethod(tool_config_backward, temp_dir)
        method.prepare_inputs(mock_tree, mock_calibrations)

        output_file = method.work_dir / "wlogdate_output.tre"
        output_file.write_text("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")

        method.execute()

        cmd = mock_runner.run.call_args[0][0]

        # 验证算法参数
        assert "-l" in cmd  # sequence_length
        assert "2000" in cmd
        assert "-s" in cmd  # seed (from common_config)
        assert "42" in cmd

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_execute_backward_time_passes_positive_root_leaf(
        self, mock_runner_class, temp_dir, mock_tree, mock_calibrations
    ):
        """测试向后时间模式下 -r 和 -f 传入正向年龄（wLogDate CLI 内部取反）"""
        config = ToolConfig()
        config.wlogdate.backward_time = True
        config.wlogdate.root_time = 90.0
        config.wlogdate.leaf_time = 0.0

        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner

        method = WLogDateMethod(config, temp_dir)
        method.prepare_inputs(mock_tree, mock_calibrations)

        output_file = method.work_dir / "wlogdate_output.tre"
        output_file.write_text("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")

        method.execute()

        cmd = mock_runner.run.call_args[0][0]

        # 验证 -r 和 -f 参数保持正向
        r_idx = cmd.index("-r")
        assert cmd[r_idx + 1] == "90.0"
        f_idx = cmd.index("-f")
        assert cmd[f_idx + 1] == "0.0"

    # ---- 输出文件读取测试 ----

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_execute_reads_output_file_not_stdout(
        self, mock_runner_class, temp_dir, tool_config, mock_tree, mock_calibrations
    ):
        """测试从输出文件读取而非 stdout"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="SHOULD_NOT_BE_USED", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner

        method = WLogDateMethod(tool_config, temp_dir)
        method.prepare_inputs(mock_tree, mock_calibrations)

        # 写入模拟输出文件
        output_file = method.work_dir / "wlogdate_output.tre"
        expected_newick = (
            "((human:0.1[chimp:0.1]:0.2):0.1,(mouse:0.3,rat:0.3):0.1)[t=65.0,mu=0.05];"
        )
        output_file.write_text(expected_newick)

        method.execute()
        result = method.parse_results()

        # 验证从文件读取，而非 stdout
        assert result.dated_tree_newick == expected_newick

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_execute_missing_output_file_raises_error(
        self, mock_runner_class, temp_dir, tool_config, mock_tree, mock_calibrations
    ):
        """测试输出文件缺失时抛出错误"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner

        method = WLogDateMethod(tool_config, temp_dir)
        method.prepare_inputs(mock_tree, mock_calibrations)

        # 不创建输出文件
        from phylodater.core.exceptions import ExecutionError

        with pytest.raises(ExecutionError, match="output file not found"):
            method.execute()

    # ---- 注册测试 ----

    def test_adapter_registered(self):
        """测试适配器已注册到 DatingMethodRegistry"""
        # 重新导入以触发注册
        import importlib

        import phylodater.adapters.wlogdate_method

        importlib.reload(phylodater.adapters.wlogdate_method)

        from phylodater.core.method_interface import DatingMethodRegistry

        assert DatingMethodRegistry.is_registered("wlogdate")

    # ---- 无校准点测试 ----

    def test_prepare_inputs_skips_empty_taxa(self, temp_dir, tool_config, mock_tree):
        """测试跳过无 resolved_taxa 的校准点:全部无效时 raise CalibrationError"""
        cal = Mock(spec=CalibrationPoint)
        cal.name = "Empty"
        cal.resolved_taxa = []
        cal.mrca_leaf_pair = None
        cal.age_constraint = FixedAgeConstraint(fixed_age=10.0)

        method = WLogDateMethod(tool_config, temp_dir)
        # 所有校准点都被跳过 → 无有效约束可写,应 raise 而非静默产出空文件
        with pytest.raises(CalibrationError):
            method.prepare_inputs(mock_tree, [cal])

    def test_prepare_inputs_skips_no_constraint(self, temp_dir, tool_config, mock_tree):
        """测试跳过无 age_constraint 的校准点:全部无效时 raise CalibrationError"""
        cal = Mock(spec=CalibrationPoint)
        cal.name = "NoConstraint"
        cal.resolved_taxa = ["human", "chimp"]
        cal.mrca_leaf_pair = ("human", "chimp")
        cal.age_constraint = None

        method = WLogDateMethod(tool_config, temp_dir)
        with pytest.raises(CalibrationError):
            method.prepare_inputs(mock_tree, [cal])


class TestWLogDateConfigSafety:
    """B-28 回归防护：适配器不得就地改写共享配置。

    此前 prepare_inputs 会把从根校准推导出的 root_time、以及从比对探测到的
    sequence_length 直接写回 config.wlogdate —— 那是跨方法、跨调用共享的实例，
    于是用户显式设置的根锚点会被静默顶掉，且副作用残留到下一次调用。
    """

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="wlogdate_safety_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @staticmethod
    def _cal(name, constraint, pair=None, is_root=False):
        return CalibrationPoint(
            name=name,
            age_constraint=constraint,
            resolved_taxa=list(pair) if pair else [],
            mrca_leaf_pair=pair,
            is_root_node=is_root,
        )

    @staticmethod
    def _tree_mock():
        tree = Mock(spec=PhylogeneticTree)
        tree.num_tips = 4
        tree.write = Mock()
        return tree

    def test_derived_root_time_is_not_written_back_to_config(self, temp_dir):
        """根校准派生的锚点只留在实例上，config.wlogdate.root_time 保持原样"""
        config = ToolConfig()
        assert config.wlogdate.root_time is None
        method = WLogDateMethod(config, temp_dir)

        method.prepare_inputs(
            self._tree_mock(),
            [
                self._cal("Root", FixedAgeConstraint(90.0), is_root=True),
                self._cal("C", UniformAgeConstraint(60.0, 80.0), ("human", "chimp")),
            ],
        )

        assert config.wlogdate.root_time is None, "shared config was mutated"
        assert method._effective_root_time == 90.0

    def test_no_cross_call_leak_of_root_time(self, temp_dir):
        """同一个共享配置上的第二次调用不会看到第一次的派生值"""
        config = ToolConfig()
        first = WLogDateMethod(config, temp_dir)
        first.prepare_inputs(
            self._tree_mock(),
            [
                self._cal("Root", FixedAgeConstraint(90.0), is_root=True),
                self._cal("C", FixedAgeConstraint(50.0), ("human", "chimp")),
            ],
        )
        assert first._effective_root_time == 90.0

        second = WLogDateMethod(config, temp_dir)
        second.prepare_inputs(
            self._tree_mock(),
            [self._cal("C", FixedAgeConstraint(50.0), ("human", "chimp"))],
        )

        assert second._effective_root_time is None

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_user_root_time_wins_and_conflict_is_warned(
        self, mock_runner_class, temp_dir
    ):
        """用户显式 root_time 与根校准冲突时：保留用户值 + 对照式警告"""
        config = ToolConfig()
        config.wlogdate.root_time = 200.0
        mock_runner = Mock()
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner
        method = WLogDateMethod(config, temp_dir)

        method.prepare_inputs(
            self._tree_mock(),
            [
                self._cal("Root", FixedAgeConstraint(90.0), is_root=True),
                self._cal("C", FixedAgeConstraint(50.0), ("human", "chimp")),
            ],
        )
        output_file = method.work_dir / "wlogdate_output.tre"
        output_file.write_text("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")
        method.execute()

        assert config.wlogdate.root_time == 200.0
        cmd = mock_runner.run.call_args[0][0]
        assert cmd[cmd.index("-r") + 1] == "200.0"
        messages = [str(w) for w in method._warnings]
        assert any(
            "conflicts" in m and "200.0" in m and "90.0" in m for m in messages
        ), messages

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_detected_sequence_length_not_written_back(
        self, mock_runner_class, temp_dir
    ):
        """从比对探测到的长度用局部值传给 -l，不写回 config.wlogdate"""
        config = ToolConfig()
        config.wlogdate.sequence_length = 500
        mock_runner = Mock()
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner
        method = WLogDateMethod(config, temp_dir)

        alignment = method.work_dir / "aln.fasta"
        alignment.write_text(">human\nACGT\n>chimp\nACGT\n")

        with patch("phylodater.infrastructure.AlignmentMetadataExtractor") as extractor:
            extractor.return_value.extract.return_value = Mock(sequence_length=44734)
            method.prepare_inputs(
                self._tree_mock(),
                [self._cal("C", FixedAgeConstraint(50.0), ("human", "chimp"))],
                alignment_path=alignment,
            )
        output_file = method.work_dir / "wlogdate_output.tre"
        output_file.write_text("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")
        method.execute()

        assert config.wlogdate.sequence_length == 500, "shared config was mutated"
        assert method._effective_sequence_length == 44734
        cmd = mock_runner.run.call_args[0][0]
        assert cmd[cmd.index("-l") + 1] == "44734"

    def test_root_time_and_config_units_are_still_consistent(self, temp_dir):
        """派生锚点在 Ga 配置下与 Ma 配置下给出同一个 -r（单位换算只做一次）"""
        cal = [
            self._cal("Root", FixedAgeConstraint(90.0), is_root=True),
            self._cal("C", FixedAgeConstraint(50.0), ("human", "chimp")),
        ]

        ma = WLogDateMethod(ToolConfig(), temp_dir)
        ma.prepare_inputs(self._tree_mock(), cal)
        assert ma._effective_root_time == pytest.approx(90.0)

        ga_config = ToolConfig()
        ga_config.wlogdate.time_unit = "Ga"
        ga = WLogDateMethod(ga_config, temp_dir)
        ga.prepare_inputs(self._tree_mock(), cal)
        assert ga._effective_root_time == pytest.approx(0.09)
        assert ga._prepare_age(ga._effective_root_time) == pytest.approx(90.0)


class TestWLogDateCliAndParsing:
    """C-53 / C-54 / C-55 / C-23 回归防护。"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="wlogdate_cli_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @staticmethod
    def _cal(name, pair=("human", "chimp"), constraint=None, is_root=False):
        return CalibrationPoint(
            name=name,
            age_constraint=constraint if constraint else FixedAgeConstraint(50.0),
            resolved_taxa=list(pair),
            mrca_leaf_pair=pair,
            is_root_node=is_root,
        )

    @staticmethod
    def _tree_mock():
        tree = Mock(spec=PhylogeneticTree)
        tree.num_tips = 4
        tree.write = Mock()
        return tree

    # ---- C-53：算法参数无条件传递 ----------------------------------------

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_algorithm_parameters_are_always_passed(self, mock_runner_class, temp_dir):
        """等于内置默认值时也必须显式传 -m/-u/-z（不再假定上游默认与本地一致）"""
        config = ToolConfig()  # max_iter=50000, pseudocount=0.01, zero_branch_len=1e-10
        mock_runner = Mock()
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner
        method = WLogDateMethod(config, temp_dir)
        method.prepare_inputs(self._tree_mock(), [self._cal("C")])
        (method.work_dir / "wlogdate_output.tre").write_text("((human:0.1):0.2);")

        method.execute()

        cmd = mock_runner.run.call_args[0][0]
        assert cmd[cmd.index("-m") + 1] == "50000"
        assert cmd[cmd.index("-u") + 1] == "0.01"
        assert cmd[cmd.index("-z") + 1] == "1e-10"

    @patch("phylodater.adapters.wlogdate_method.ProcessRunner")
    def test_effective_parameters_recorded_in_metadata(
        self, mock_runner_class, temp_dir
    ):
        """实际生效值写进 metadata（而非配置里的假设值）"""
        config = ToolConfig()
        mock_runner = Mock()
        mock_runner.run.return_value = Mock(
            returncode=0, stdout="", stderr="", execution_time=1.0
        )
        mock_runner_class.return_value = mock_runner
        method = WLogDateMethod(config, temp_dir)
        method.prepare_inputs(self._tree_mock(), [self._cal("C")])
        (method.work_dir / "wlogdate_output.tre").write_text(
            "((human:0.1,chimp:0.1)C[t=40.0,mu=0.05]:0.2,(mouse:0.3,rat:0.3):0.1)Root[t=65.0,mu=0.03];"
        )
        method.execute()
        method._execution_seconds = 12.5

        result = method.parse_results()

        assert result.execution_seconds == 12.5
        assert result.metadata["execution_time"] == 12.5
        assert result.metadata["max_iter"] == config.wlogdate.max_iter
        assert result.metadata["effective_sequence_length"] is not None

    # ---- C-55：负年龄不再静默取绝对值 ------------------------------------

    def test_negative_age_on_calibration_is_not_absoluted(self, temp_dir):
        """t=-100 必须留下负值 + 显式披露（旧实现折成 mean_age=100.0 且无警告）"""
        method = WLogDateMethod(ToolConfig(), temp_dir)
        method._calibrations = [
            self._cal("CladeAB", constraint=FixedAgeConstraint(100.0))
        ]
        newick = (
            "((human:0.1,chimp:0.1)CladeAB[t=-100.0,mu=0.05]:0.2,"
            "(mouse:0.3,rat:0.3):0.15)Root[t=200.0,mu=0.03];"
        )

        with patch.object(method.logger, "error") as error_log:
            node_ages = method._parse_tree_ages(newick)

        assert node_ages["CladeAB"].mean_age == -100.0
        assert error_log.called
        assert any("NEGATIVE age" in str(w) for w in method._warnings)
        assert not any(
            node.mean_age == 100.0 for node in node_ages.values()
        ), "negative age was silently folded to a plausible positive"

    def test_negative_age_on_unnamed_internal_node_is_not_absoluted(self, temp_dir):
        """未命名内部节点那一处（第三个 abs 站点）同样保留符号并披露"""
        method = WLogDateMethod(ToolConfig(), temp_dir)
        method._calibrations = []
        newick = (
            "((human:0.1,chimp:0.1)[t=-50.0,mu=0.05]:0.2,"
            "(mouse:0.3,rat:0.3):0.15)Root[t=65.0,mu=0.03];"
        )

        with patch.object(method.logger, "error") as error_log:
            node_ages = method._parse_tree_ages(newick)

        assert node_ages["internal_node_1"].mean_age == -50.0
        assert error_log.called

    def test_negative_age_warning_reaches_dating_result(self, temp_dir):
        """负年龄的披露进入 DatingResult.warnings，报告里看得见"""
        config = ToolConfig()
        method = WLogDateMethod(config, temp_dir)
        method._calibrations = [self._cal("CladeAB")]
        method._input_files["output_tree"] = method.work_dir / "out.tre"
        (method.work_dir / "out.tre").write_text(
            "((human:0.1,chimp:0.1)CladeAB[t=-100.0,mu=0.05]:0.2,"
            "(mouse:0.3,rat:0.3):0.15)Root[t=200.0,mu=0.03];"
        )

        with patch.object(method.logger, "error"):
            result = method.parse_results()

        assert any("NEGATIVE age" in str(w) for w in result.warnings)

    def test_positive_ages_produce_no_negative_age_warning(self, temp_dir):
        """正常（正向）输出不应产生任何噪声警告"""
        method = WLogDateMethod(ToolConfig(), temp_dir)
        method._calibrations = [self._cal("CladeAB")]
        newick = (
            "((human:0.1,chimp:0.1)CladeAB[t=30.0,mu=0.05]:0.2,"
            "(mouse:0.3,rat:0.3):0.15)Root[t=65.0,mu=0.03];"
        )

        node_ages = method._parse_tree_ages(newick)

        assert node_ages["CladeAB"].mean_age == 30.0
        assert method._warnings == []

    # ---- C-23：解析集中渲染结果时对制表符鲁棒 ----------------------------

    def test_age_parsed_from_last_tab_field(self, temp_dir):
        """taxa_str 里含制表符时仍能取回年龄（rsplit 语义），不再 ValueError"""
        method = WLogDateMethod(ToolConfig(), temp_dir)
        constraint = FixedAgeConstraint(65.0)

        formatted = constraint.to_software_format(
            "wlogdate", name="Clade", taxa_str="human\tchimp"
        )
        assert method._age_from_wlogdate_line(formatted, "Clade") == 65.0

    def test_malformed_rendering_raises_calibration_error(self, temp_dir):
        """渲染结果形状不符时显式报错，而不是静默丢点"""
        method = WLogDateMethod(ToolConfig(), temp_dir)

        with pytest.raises(CalibrationError, match="unexpected rendering"):
            method._age_from_wlogdate_line("no tab here", "Clade")
        with pytest.raises(CalibrationError, match="non-numeric"):
            method._age_from_wlogdate_line("Clade=human+chimp\tnot_a_number", "Clade")


class TestWLogDateRealOutputParsing:
    """回归防护：真实 wLogDate 上游输出与单元测试样例有两个关键差异，
    早期实现因此在实际集成测试上解析出 0 个节点年龄：

    1. 真实注解是 ', mu=...'（逗号后有空格），单元测试样例是 ',mu=...'（无空格）
       -> 旧正则一个注解都匹配不到（本仓库集成测试在 benchmark 数据上复现）；
    2. 真实树里的 )Name[t=..., mu=...] 注解会让 ete3 解析器直接抛
       'NoneType' object has no attribute 'add_child'，MRCA 回退应使用
       模块内的轻量解析器（_parse_wlogdate_newick + _wlogdate_lca）。
    """

    # 一次真实 wLogDate 运行的输出（6 条灵长类，HumanChimp 固定 6.0、
    # GreatApes 均匀 8-10 降级 9.0、Root 最大 15 降级 7.5；节点 t 由上游重估）
    REAL_OUTPUT = (
        "((human:6.0,(chimp:1.1231608379607356,bonobo:1.1231608379607356)"
        "I2[t=1.123160837960735, mu=0.004101016936477549]:4.876839162039265)"
        "HumanChimp[t=6.0, mu=0.03333333333333333]:1.5,"
        "(gorilla:5.652568701593696,"
        "(orangutan:3.779519909958468,sumatran:3.779519909958468)"
        "I4[t=3.779519909958468, mu=0.010677778437655865]:1.8730487916352272)"
        "I3[t=5.652568701593696, mu=0.010825842355952884]:1.8474312984063042)"
        "Root[t=7.5, mu=0.010351687777664391]:0.02;"
    )

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="wlogdate_real_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @staticmethod
    def _method(temp_dir):
        return WLogDateMethod(ToolConfig(), temp_dir)

    def test_regex_accepts_real_space_before_mu(self, temp_dir):
        """真实格式 ', mu='（带空格）必须全部匹配（旧正则返回 0 个匹配）"""
        from phylodater.adapters.wlogdate_method import _ANNOTATION_RE

        matches = [
            (m.group(1), float(m.group(2)))
            for m in _ANNOTATION_RE.finditer(self.REAL_OUTPUT)
        ]
        assert matches == [
            ("I2", 1.123160837960735),
            ("HumanChimp", 6.0),
            ("I4", 3.779519909958468),
            ("I3", 5.652568701593696),
            ("Root", 7.5),
        ]

    def test_parse_tree_ages_on_real_output(self, temp_dir):
        """真实输出 + benchmark 校准集 -> 3 个校准点 + 3 个未命名内部节点都有年龄"""
        method = self._method(temp_dir)
        method._calibrations = [
            CalibrationPoint(
                name="HumanChimp",
                age_constraint=FixedAgeConstraint(6.0),
                mrca_leaf_pair=("human", "chimp"),
                resolved_taxa=["human", "chimp"],
            ),
            CalibrationPoint(
                name="GreatApes",
                age_constraint=UniformAgeConstraint(8.0, 10.0),
                mrca_leaf_pair=("human", "gorilla"),
                resolved_taxa=["human", "gorilla"],
            ),
            CalibrationPoint(
                name="Root",
                age_constraint=MaximumAgeConstraint(15.0),
                is_root_node=True,
            ),
        ]

        node_ages = method._parse_tree_ages(self.REAL_OUTPUT)

        # 精确标签匹配
        assert node_ages["HumanChimp"].mean_age == pytest.approx(6.0)
        assert node_ages["Root"].mean_age == pytest.approx(7.5)
        # MRCA 回退：human+gorilla 的 MRCA 是 Root 节点（标签已被 Root 校准占用）
        assert node_ages["GreatApes"].mean_age == pytest.approx(7.5)
        # 未命名内部节点
        assert len(node_ages) == 6
        assert node_ages["internal_node_1"].mean_age == pytest.approx(1.123160837960735)
        assert node_ages["internal_node_2"].mean_age == pytest.approx(3.779519909958468)
        assert node_ages["internal_node_3"].mean_age == pytest.approx(5.652568701593696)
        # 无噪声警告
        assert method._warnings == []

    def test_mrca_fallback_without_any_label(self, temp_dir):
        """输出树完全没有校准节点名时，MRCA 回退仍能按 mrca_leaf_pair 取到年龄"""
        method = self._method(temp_dir)
        method._calibrations = [
            CalibrationPoint(
                name="Apes",
                age_constraint=FixedAgeConstraint(10.0),
                mrca_leaf_pair=("human", "gorilla"),
                resolved_taxa=["human", "gorilla"],
            ),
        ]
        # 同一棵拓扑，但把校准名改成 A/B/C（不匹配任何标签）
        renamed = self.REAL_OUTPUT.replace("HumanChimp", "N1").replace("Root", "N2")
        node_ages = method._parse_tree_ages(renamed)

        assert "Apes" in node_ages
        assert node_ages["Apes"].mean_age == pytest.approx(7.5)  # MRCA = 根节点

    def test_lightweight_parser_structure(self, temp_dir):
        """轻量解析器正确还原拓扑、叶子名与注解块"""
        from phylodater.adapters.wlogdate_method import (
            _parse_wlogdate_newick,
            _wlogdate_find_leaf,
            _wlogdate_lca,
        )

        root = _parse_wlogdate_newick(self.REAL_OUTPUT)
        assert root.name == "Root"
        assert root.annotation == "[t=7.5, mu=0.010351687777664391]"
        assert sorted(_leaf_names(root)) == [
            "bonobo",
            "chimp",
            "gorilla",
            "human",
            "orangutan",
            "sumatran",
        ]
        n_human = _wlogdate_find_leaf(root, "human")
        n_gorilla = _wlogdate_find_leaf(root, "gorilla")
        mrca = _wlogdate_lca(n_human, n_gorilla)
        assert mrca is root
        # 分支长度与注解共存时拓扑仍然正确
        n_chimp = _wlogdate_find_leaf(root, "chimp")
        n_bonobo = _wlogdate_find_leaf(root, "bonobo")
        assert _wlogdate_lca(n_chimp, n_bonobo).name == "I2"
        assert (
            _wlogdate_lca(
                _wlogdate_find_leaf(root, "orangutan"),
                _wlogdate_find_leaf(root, "sumatran"),
            ).name
            == "I4"
        )

    def test_parser_rejects_malformed_newick(self, temp_dir):
        """畸形树抛 ResultParsingError，调用方降级为仅标签匹配而非崩溃"""
        from phylodater.adapters.wlogdate_method import _parse_wlogdate_newick
        from phylodater.core.exceptions import ResultParsingError

        method = self._method(temp_dir)
        method._calibrations = [
            CalibrationPoint(
                name="X",
                age_constraint=FixedAgeConstraint(5.0),
                mrca_leaf_pair=("human", "chimp"),
                resolved_taxa=["human", "chimp"],
            ),
        ]
        with pytest.raises(ResultParsingError):
            _parse_wlogdate_newick("((human:0.1,chimp:0.1):0.2, (no closing)")

        # 方法级：畸形树不抛异常，只警告（保持结果可用）
        node_ages = method._parse_tree_ages("(unbalanced[t=5.0, mu=0.1]:0.2;")
        assert isinstance(node_ages, dict)


def _leaf_names(node) -> list:
    """测试辅助：收集子树全部叶节点名（供轻量解析器断言）。"""
    if not node.children:
        return [node.name]
    names = []
    for child in node.children:
        names.extend(_leaf_names(child))
    return names
