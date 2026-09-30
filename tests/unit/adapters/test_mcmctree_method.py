"""
MCMCTreeMethod 适配器单元测试

测试 MCMCTree 定年方法的各项功能
"""

import os
import shutil
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.adapters.mcmctree_method import MCMCTreeMethod
from phylodater.core.exceptions import CalibrationError
from phylodater.infrastructure.configuration import CommonConfig, MCMCTreeConfig
from phylodater.models import (
    CalibrationPoint,
    FixedAgeConstraint,
    PhylogeneticTree,
    UniformAgeConstraint,
)


@pytest.fixture
def temp_dir():
    """模块级临时目录 fixture（供本文件所有测试类使用）。"""
    temp_path = Path(tempfile.mkdtemp(prefix="mcmctree_test_"))
    yield temp_path
    shutil.rmtree(temp_path, ignore_errors=True)


class TestMCMCTreeMethod:
    """MCMCTreeMethod 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="mcmctree_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def mcmctree_config(self):
        """创建 MCMCTree 配置"""
        return MCMCTreeConfig(
            clock=2,
            burnin=1000,
            nsample=2000,
            num_runs=1,
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
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        tree.num_tips = 4
        tree.with_calibration_annotations = Mock(return_value=tree)
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

        calib2 = Mock(spec=CalibrationPoint)
        calib2.name = "Root"
        calib2.is_root_node = True
        calib2.age_constraint = FixedAgeConstraint(fixed_age=100.0)

        return [calib1, calib2]

    def test_init_creates_work_directory(
        self, temp_dir, mcmctree_config, common_config
    ):
        """测试初始化创建工作目录"""
        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)

        assert method.work_dir.exists()
        assert method.method_name == "mcmctree"

    def test_method_name_property(self, temp_dir, mcmctree_config, common_config):
        """测试方法名称属性"""
        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)
        assert method.method_name == "mcmctree"

    def test_common_config_usage(self, temp_dir, mcmctree_config, common_config):
        """测试 common_config 使用"""
        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)
        assert method.common_config.nthreads == 4
        assert method.common_config.seed == 42

    @patch("phylodater.adapters.mcmctree_method.ProcessRunner")
    def test_validate_environment_with_executable(
        self, mock_runner_class, temp_dir, mcmctree_config, common_config
    ):
        """测试环境验证 - MCMCTree 存在且版本匹配"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = (
            "MCMCTREE in paml version 4.10.8, 15 May 2025"
        )
        mock_runner_class.return_value = mock_runner

        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is True

    @patch("phylodater.adapters.mcmctree_method.ProcessRunner")
    def test_validate_environment_version_mismatch(
        self, mock_runner_class, temp_dir, mcmctree_config, common_config
    ):
        """测试环境验证 - 跨版本系列（4.8 vs 4.10）必须阻断"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = (
            "MCMCTREE in paml version 4.8, 15 May 2025"
        )
        mock_runner_class.return_value = mock_runner

        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is False

    @patch("phylodater.adapters.mcmctree_method.ProcessRunner")
    def test_validate_environment_same_series_other_patch_is_accepted(
        self, mock_runner_class, temp_dir, mcmctree_config, common_config
    ):
        """同一发布系列内的补丁号差异不得阻断运行。

        控制文件格式/超先语法的差异发生在 4.8 与 4.10.x 之间，而不是 4.10.7 /
        4.10.8 / 4.10.10 之间。旧实现要求逐字相等，于是按文档用
        ``conda install -c bioconda paml`` 装到的 4.10.10 会被直接拒跑。
        """
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = (
            "MCMCTREE in paml version 4.10.10, 29 Jan 2026"
        )
        mock_runner_class.return_value = mock_runner

        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)
        assert method.validate_environment() is True
        assert method._detected_paml_version == "4.10.10"
        # 引擎版本必须被记录，随包发布的结果才有可溯源性
        assert "4.10.10" in (method.software_version or "")

    def test_paml_versions_compatible_series_semantics(self):
        """版本兼容判据：主.次版本系列相同即可，不同则阻断。"""
        from phylodater.adapters.mcmctree_method import _paml_versions_compatible

        assert _paml_versions_compatible("4.10.8", "4.10.10") is True
        assert _paml_versions_compatible("4.10.8", "4.10.7") is True
        assert _paml_versions_compatible("4.10.8", "4.8") is False
        assert _paml_versions_compatible("4.8", "4.9.10") is False
        # 无法解析时一律不阻断（交由上层日志提示）
        assert _paml_versions_compatible("", "4.10.10") is True
        assert _paml_versions_compatible("4.10.8", "unknown") is True

    @patch("phylodater.adapters.mcmctree_method.ProcessRunner")
    def test_validate_environment_without_executable(
        self, mock_runner_class, temp_dir, mcmctree_config, common_config
    ):
        """测试环境验证 - MCMCTree 不存在"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = False
        mock_runner_class.return_value = mock_runner

        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is False

    def test_get_executable_uses_paml_path(
        self, temp_dir, mcmctree_config, common_config
    ):
        """测试 _get_executable 优先从 paml_path/bin 拼接可执行文件"""
        from phylodater.adapters.mcmctree_method import _has_path_separator
        from phylodater.infrastructure.configuration import SoftwarePaths

        with tempfile.TemporaryDirectory() as paml_root:
            bin_dir = Path(paml_root) / "bin"
            bin_dir.mkdir()
            mcmctree_bin = bin_dir / "mcmctree"
            mcmctree_bin.write_text("#!/bin/sh\necho test")
            mcmctree_bin.chmod(0o755)

            software_paths = SoftwarePaths(paml_path=paml_root)
            method = MCMCTreeMethod(
                mcmctree_config,
                temp_dir,
                software_paths=software_paths,
                common_config=common_config,
            )
            exe = method._get_executable("mcmctree")
            assert _has_path_separator(exe)
            assert Path(exe).resolve() == mcmctree_bin.resolve()

    def test_get_executable_explicit_full_path_wins(
        self, temp_dir, mcmctree_config, common_config
    ):
        """测试 _get_executable 中显式完整路径优先于 paml_path"""
        from phylodater.infrastructure.configuration import SoftwarePaths

        with (
            tempfile.TemporaryDirectory() as paml_root,
            tempfile.TemporaryDirectory() as other_root,
        ):
            (Path(paml_root) / "bin").mkdir()
            default_bin = Path(paml_root) / "bin" / "mcmctree"
            default_bin.write_text("#!/bin/sh\necho default")
            default_bin.chmod(0o755)

            explicit_bin = Path(other_root) / "mcmctree"
            explicit_bin.write_text("#!/bin/sh\necho explicit")
            explicit_bin.chmod(0o755)

            software_paths = SoftwarePaths(
                paml_path=paml_root, mcmctree_bin=str(explicit_bin)
            )
            method = MCMCTreeMethod(
                mcmctree_config,
                temp_dir,
                software_paths=software_paths,
                common_config=common_config,
            )
            assert method._get_executable("mcmctree") == str(explicit_bin)

    def test_get_paml_env_prepends_bin_dir(
        self, temp_dir, mcmctree_config, common_config
    ):
        """测试 _get_paml_env 将 paml_path/bin 前置到 PATH"""
        from phylodater.infrastructure.configuration import SoftwarePaths

        with tempfile.TemporaryDirectory() as paml_root:
            bin_dir = Path(paml_root) / "bin"
            bin_dir.mkdir()
            fake_bin = bin_dir / "mcmctree"
            fake_bin.write_text("#!/bin/sh\necho test")
            fake_bin.chmod(0o755)

            software_paths = SoftwarePaths(paml_path=paml_root)
            method = MCMCTreeMethod(
                mcmctree_config,
                temp_dir,
                software_paths=software_paths,
                common_config=common_config,
            )
            env = method._get_paml_env()
            path_parts = env["PATH"].split(os.pathsep)
            assert str(bin_dir.resolve()) == path_parts[0]

    def test_get_paml_env_falls_back_to_exe_parent(
        self, temp_dir, mcmctree_config, common_config
    ):
        """测试未设置 paml_path 但 mcmctree_bin 为完整路径时，_get_paml_env 使用其目录"""
        from phylodater.infrastructure.configuration import SoftwarePaths

        with tempfile.TemporaryDirectory() as exe_dir:
            fake_bin = Path(exe_dir) / "mcmctree"
            fake_bin.write_text("#!/bin/sh\necho test")
            fake_bin.chmod(0o755)

            software_paths = SoftwarePaths(mcmctree_bin=str(fake_bin))
            method = MCMCTreeMethod(
                mcmctree_config,
                temp_dir,
                software_paths=software_paths,
                common_config=common_config,
            )
            env = method._get_paml_env()
            path_parts = env["PATH"].split(os.pathsep)
            assert str(fake_bin.parent.resolve()) == path_parts[0]

    def test_config_parameters(self, temp_dir, mcmctree_config, common_config):
        """测试配置参数"""
        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)

        assert method.config.clock == 2
        assert method.config.burnin == 1000
        assert method.config.nsample == 2000
        assert method.config.num_runs == 1

    def test_has_root_calibration_detection(
        self, temp_dir, mcmctree_config, common_config, mock_calibrations
    ):
        """测试根节点校准检测"""
        MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)

        # 检查校准点中是否有根节点
        has_root = any(cal.is_root_node for cal in mock_calibrations)
        assert has_root is True

    def test_checkpoint_manager_initialization(
        self, temp_dir, mcmctree_config, common_config
    ):
        """测试检查点管理器初始化"""
        method = MCMCTreeMethod(
            mcmctree_config,
            temp_dir,
            enable_checkpoint=True,
            common_config=common_config,
        )
        assert method._enable_checkpoint is True

    def test_work_dir_isolation(self, temp_dir, mcmctree_config, common_config):
        """测试工作目录隔离"""
        method1 = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)
        method2 = MCMCTreeMethod(
            mcmctree_config, temp_dir / "other", common_config=common_config
        )

        assert method1.work_dir != method2.work_dir

    def test_validate_calibrations_rejects_fixed_root(
        self, temp_dir, mcmctree_config, common_config, mock_calibrations
    ):
        """校验固定根节点年龄会被 MCMCTree 拦截"""
        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)
        errors = method.validate_calibrations(mock_calibrations)
        assert len(errors) >= 1
        assert "不支持固定根节点年龄" in errors[0]

    def test_validate_calibrations_accepts_uniform_root(
        self, temp_dir, mcmctree_config, common_config
    ):
        """校验 uniform 根节点约束可通过 MCMCTree 检查"""
        method = MCMCTreeMethod(mcmctree_config, temp_dir, common_config=common_config)

        root_cal = Mock(spec=CalibrationPoint)
        root_cal.name = "Root"
        root_cal.is_root_node = True
        root_cal.age_constraint = UniformAgeConstraint(min_age=499.9, max_age=500.1)

        errors = method.validate_calibrations([root_cal])
        assert errors == []


class TestMCMCTreeConfig:
    """MCMCTree 配置测试"""

    def test_default_values(self):
        """测试默认值"""
        config = MCMCTreeConfig()
        assert config.clock == 2
        assert config.num_runs == 2
        assert config.burnin == 20000
        assert config.nsample == 50000

    def test_custom_params(self):
        """测试自定义参数"""
        config = MCMCTreeConfig(clock=3, num_runs=4, burnin=50000)
        assert config.clock == 3
        assert config.num_runs == 4
        assert config.burnin == 50000

    def test_clock_validation(self):
        """测试 clock 验证"""
        with pytest.raises(ValueError):
            MCMCTreeConfig(clock=4)

    def test_burnin_validation(self):
        """测试 burnin 验证"""
        with pytest.raises(ValueError):
            MCMCTreeConfig(burnin=-1)


class TestMCMCTreeAnnotationParsing:
    """MCMCTree FigTree/NEXUS 节点年龄解析测试"""

    def test_annotation_to_age_paml_48_format(self):
        """测试解析 PAML 4.8 的 [&95%={lower,upper}] 注解"""
        from dendropy import Tree

        from phylodater.adapters.mcmctree_method import _annotation_to_age

        tree_data = (
            "((A:0.1,B:0.1):0.2[&95%={0.1,0.3}],(C:0.1,D:0.1):0.3[&95%={0.2,0.4}]);"
        )
        tree = Tree.get(data=tree_data, schema="newick", preserve_underscores=True)

        ages = []
        for node in tree.preorder_node_iter():
            age = _annotation_to_age(node)
            if age:
                ages.append(age)

        assert len(ages) == 2
        assert ages[0]["mean"] == pytest.approx(0.2)
        assert ages[0]["lower"] == pytest.approx(0.1)
        assert ages[0]["upper"] == pytest.approx(0.3)
        assert ages[1]["mean"] == pytest.approx(0.3)
        assert ages[1]["lower"] == pytest.approx(0.2)
        assert ages[1]["upper"] == pytest.approx(0.4)

    def test_annotation_to_age_paml_4108_format(self):
        """测试解析 PAML 4.10.8+ 的 [&95%HPD={lower,upper}] 注解"""
        from dendropy import Tree

        from phylodater.adapters.mcmctree_method import _annotation_to_age

        tree_data = "((A:0.1,B:0.1):0.2[&95%HPD={0.1,0.3}],(C:0.1,D:0.1):0.3[&95%HPD={0.2,0.4}]);"
        tree = Tree.get(data=tree_data, schema="newick", preserve_underscores=True)

        ages = []
        for node in tree.preorder_node_iter():
            age = _annotation_to_age(node)
            if age:
                ages.append(age)

        assert len(ages) == 2
        assert ages[0]["mean"] == pytest.approx(0.2)
        assert ages[1]["mean"] == pytest.approx(0.3)

    def test_parse_figtree_nexus_paml_48_tree_keyword(self):
        """测试 _parse_figtree_nexus 同时兼容 UTREE 和 TREE 关键字"""
        from phylodater.adapters.mcmctree_method import _parse_figtree_nexus
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        nexus = """#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS A B C D;
END;
BEGIN TREES;
    TREE 1 = ((A:0.1,B:0.1):0.2[&95%={0.1,0.3}],(C:0.1,D:0.1):0.3[&95%={0.2,0.4}]);
END;
"""
        cal = CalibrationPoint(
            name="AB",
            mrca_leaf_pair=("A", "B"),
            age_constraint=FixedAgeConstraint(fixed_age=0.2),
        )
        ages = _parse_figtree_nexus(nexus, [cal], logger=None)

        assert "AB" in ages
        assert ages["AB"]["mean"] == pytest.approx(0.2)
        assert ages["AB"]["lower"] == pytest.approx(0.1)
        assert ages["AB"]["upper"] == pytest.approx(0.3)
        assert "internal_node_1" in ages
        assert ages["internal_node_1"]["mean"] == pytest.approx(0.3)


# --------------------------------------------------------------------------- #
# 多链汇总的区间标注（审阅报告 B-3）
# --------------------------------------------------------------------------- #
class TestMCMCTreeMultiChainSummaryLabelling:
    """多链汇总：ci_type 必须与实际算出的量一致，mean/median 必须同池。"""

    def _write_run(self, method, run_id, figtree_text):
        run_dir = method.work_dir / f"run{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "FigTree.tre").write_text(figtree_text, encoding="utf-8")

    def _method(self, temp_dir, num_runs):
        config = MCMCTreeConfig(num_runs=num_runs, skip_rate_estimation=True)
        method = MCMCTreeMethod(
            config, temp_dir, enable_checkpoint=False, common_config=CommonConfig()
        )
        method._calibrations = [
            CalibrationPoint(
                name="AB",
                mrca_leaf_pair=("A", "B"),
                age_constraint=UniformAgeConstraint(min_age=100.0, max_age=300.0),
            ),
            CalibrationPoint(
                name="CD",
                mrca_leaf_pair=("C", "D"),
                age_constraint=UniformAgeConstraint(min_age=200.0, max_age=400.0),
            ),
        ]
        return method

    def test_single_chain_interval_is_labelled_hpd95(self, temp_dir):
        """只有 1 条链报告该节点时，区间就是 MCMCTree 自己的 95% HPD。"""
        from phylodater.models import CIType

        method = self._method(temp_dir, num_runs=1)
        self._write_run(
            method,
            1,
            "((A:0.1,B:0.1):0.2[&95%HPD={0.1,0.3}],"
            "(C:0.1,D:0.1):0.3[&95%HPD={0.2,0.4}]);",
        )

        result = method.parse_results()

        assert result.node_ages["AB"].ci_type is CIType.HPD95
        assert result.node_ages["AB"].mean_age == pytest.approx(200.0)
        assert result.node_ages["AB"].ci_lower == pytest.approx(100.0)
        assert result.node_ages["AB"].ci_upper == pytest.approx(300.0)
        assert result.metadata["node_age_interval_kinds"]["AB"] == (
            "single-chain 95% HPD as reported by MCMCTree"
        )

    def test_multi_chain_envelope_is_not_labelled_hpd95(self, temp_dir):
        """2 条链的外包络不得再标 HPD95；mean/median 取同一批链。"""
        from phylodater.models import CIType

        method = self._method(temp_dir, num_runs=2)
        self._write_run(
            method,
            1,
            "((A:0.1,B:0.1):0.2[&95%HPD={0.1,0.3}],"
            "(C:0.1,D:0.1):0.3[&95%HPD={0.2,0.4}]);",
        )
        self._write_run(
            method,
            2,
            "((A:0.1,B:0.1):0.2[&95%HPD={0.15,0.35}],"
            "(C:0.1,D:0.1):0.3[&95%HPD={0.25,0.45}]);",
        )

        result = method.parse_results()

        ab = result.node_ages["AB"]
        # 旧实现：ci_type=HPD95（冒称 HPD）、median 只取 run1（0.2）
        assert ab.ci_type is CIType.RANGE
        assert ab.mean_age == pytest.approx(225.0)  # (0.2+0.25)/2 × 1000
        assert ab.median_age == pytest.approx(225.0)  # 与 mean 同一批链
        assert ab.ci_lower == pytest.approx(100.0)  # min(0.1, 0.15)
        assert ab.ci_upper == pytest.approx(350.0)  # max(0.3, 0.35)
        assert "outer envelope" in result.metadata["node_age_interval_kinds"]["AB"]
        assert "NOT a 95% HPD" in result.metadata["node_age_interval_kinds"]["AB"]

    def test_disjoint_chain_intervals_are_reported_as_conflict(self, temp_dir):
        """两条链的 95% 区间完全不相交时必须暴露，而不是被外包络抹平。"""
        method = self._method(temp_dir, num_runs=2)
        recorder = Mock()
        method.logger = recorder
        self._write_run(
            method,
            1,
            "((A:0.1,B:0.1):0.2[&95%HPD={0.1,0.2}],"
            "(C:0.1,D:0.1):0.3[&95%HPD={0.2,0.3}]);",
        )
        self._write_run(
            method,
            2,
            "((A:0.1,B:0.1):0.2[&95%HPD={0.5,0.7}],"
            "(C:0.1,D:0.1):0.3[&95%HPD={0.4,0.6}]);",
        )

        result = method.parse_results()

        assert result.metadata["chain_interval_conflicts"]["AB"] == [[1, 2]]
        assert "AB" not in result.metadata.get("chain_interval_conflicts", {}).get(
            "CD", []
        )
        errors = " ".join(
            " ".join(str(a) for a in c.args) for c in recorder.error.call_args_list
        )
        assert "Multi-chain conflict" in errors or "do not overlap" in errors

    def test_node_reported_by_fewer_runs_is_warned(self, temp_dir):
        """num_runs=2 但只有 run1 有 FigTree.tre 时必须告警（不再静默按 run1 取值）。"""
        from phylodater.models import CIType

        method = self._method(temp_dir, num_runs=2)
        recorder = Mock()
        method.logger = recorder
        self._write_run(
            method,
            1,
            "((A:0.1,B:0.1):0.2[&95%HPD={0.1,0.3}],"
            "(C:0.1,D:0.1):0.3[&95%HPD={0.2,0.4}]);",
        )

        result = method.parse_results()

        assert result.node_ages["AB"].ci_type is CIType.HPD95
        assert result.metadata["runs_with_figtree"] == 1
        warnings = " ".join(
            " ".join(str(a) for a in c.args) for c in recorder.warning.call_args_list
        )
        assert "only 1 of 2 configured runs" in warnings


# --------------------------------------------------------------------------- #
# PAML .dat 模型矩阵（审阅报告 B-5）
# --------------------------------------------------------------------------- #
# 与随附 PAML 4.8 / 4.10.8 的 dat/ 目录逐字一致（两版本内容相同，共 18 个文件）。
UPSTREAM_PAML_DAT_FILES = [
    "MtZoa.dat",
    "cpREV10.dat",
    "cpREV64.dat",
    "dayhoff-dcmut.dat",
    "dayhoff.dat",
    "g1974a.dat",
    "g1974c.dat",
    "g1974p.dat",
    "g1974v.dat",
    "grantham.dat",
    "jones-dcmut.dat",
    "jones.dat",
    "lg.dat",
    "miyata.dat",
    "mtArt.dat",
    "mtREV24.dat",
    "mtmam.dat",
    "wag.dat",
]


class TestMCMCTreePamlDatFiles:
    """.dat 复制必须按上游真实文件名（大小写敏感通路），失败必须可见。"""

    def _install_fake_paml(self, root, dat_files):
        bin_dir = Path(root) / "bin"
        bin_dir.mkdir(parents=True, exist_ok=True)
        exe = bin_dir / "mcmctree"
        exe.write_text("#!/bin/sh\necho 'MCMCTREE in paml version 4.10.8'\n")
        exe.chmod(0o755)
        dat_dir = Path(root) / "dat"
        dat_dir.mkdir(exist_ok=True)
        for name in dat_files:
            (dat_dir / name).write_text("fake matrix\n", encoding="utf-8")
        return dat_dir

    def _method_for(self, temp_dir, root):
        from phylodater.infrastructure.configuration import SoftwarePaths

        method = MCMCTreeMethod(
            MCMCTreeConfig(),
            temp_dir,
            software_paths=SoftwarePaths(paml_path=str(root)),
            common_config=CommonConfig(),
        )
        return method

    def test_expected_list_matches_upstream_file_names(self):
        """清单必须与上游 dat/ 目录**逐字**一致（含大小写）。

        这条断言就是跨平台契约：旧清单里的 mtART.dat / mtZOA.dat 在 Linux 上
        永远复制不到，而 mtMet.dat / cpREV45.dat / aamodel.dat / figual.dat
        上游根本不存在，常用的 lg.dat 反而漏了。
        """
        assert set(MCMCTreeMethod._PAML_DAT_FILES) == set(UPSTREAM_PAML_DAT_FILES)
        assert "mtArt.dat" in MCMCTreeMethod._PAML_DAT_FILES
        assert "MtZoa.dat" in MCMCTreeMethod._PAML_DAT_FILES
        assert "lg.dat" in MCMCTreeMethod._PAML_DAT_FILES
        for phantom in (
            "mtART.dat",
            "mtZOA.dat",
            "mtMet.dat",
            "cpREV45.dat",
            "aamodel.dat",
            "figual.dat",
        ):
            assert phantom not in MCMCTreeMethod._PAML_DAT_FILES

    def test_every_expected_entry_is_matchable_case_insensitively(self, temp_dir):
        """大小写不敏感的匹配检查（在大小写敏感的 Linux 上同样成立）。"""
        self._install_fake_paml(temp_dir / "paml", UPSTREAM_PAML_DAT_FILES)
        method = self._method_for(temp_dir, temp_dir / "paml")

        audit = method._audit_paml_dat_files()

        assert audit["dat_dir"] == (temp_dir / "paml" / "dat").resolve()
        assert audit["missing"] == []
        assert audit["case_fixed"] == {}
        assert audit["required_missing"] == []

    def test_copy_all_dat_files_with_real_names(self, temp_dir):
        """复制走整目录枚举，文件名保留上游大小写。"""
        self._install_fake_paml(temp_dir / "paml", UPSTREAM_PAML_DAT_FILES)
        method = self._method_for(temp_dir, temp_dir / "paml")

        method._copy_paml_dat_files()

        copied = sorted(p.name for p in method.work_dir.glob("*.dat"))
        assert copied == sorted(UPSTREAM_PAML_DAT_FILES)
        # 按源文件真实大小写落盘（不能用 work_dir/"mtART.dat".exists() 断言：
        # APFS 默认大小写不敏感，两种拼写指向同一文件，只有枚举真实名字才有效）
        assert "mtArt.dat" in copied
        assert "MtZoa.dat" in copied
        assert method._missing_dat_files == []
        assert method._required_dat_files_missing == []

    def test_missing_required_dat_is_error_and_fails_environment(self, temp_dir):
        """控制文件依赖的 wag.dat 缺失 → ERROR + 环境检查失败（不再只记 debug）。"""
        files = [f for f in UPSTREAM_PAML_DAT_FILES if f != "wag.dat"]
        self._install_fake_paml(temp_dir / "paml", files)
        method = self._method_for(temp_dir, temp_dir / "paml")
        recorder = Mock()
        method.logger = recorder

        method._copy_paml_dat_files()

        assert method._required_dat_files_missing == ["wag.dat"]
        errors = " ".join(
            " ".join(str(a) for a in c.args) for c in recorder.error.call_args_list
        )
        assert "wag.dat" in errors

        with patch("phylodater.adapters.mcmctree_method.ProcessRunner") as runner_cls:
            runner = Mock()
            runner.check_executable.return_value = True
            runner.get_version.return_value = (
                "MCMCTREE in paml version 4.10.8, 15 May 2025"
            )
            runner_cls.return_value = runner
            assert method.validate_environment() is False

    def test_missing_dat_directory_logs_warning_not_debug(self, temp_dir):
        """dat 目录找不到时至少 WARNING（旧实现只记 debug，用户与 CI 都看不到）。"""
        from phylodater.infrastructure.configuration import SoftwarePaths

        empty_root = temp_dir / "nodir"
        (empty_root / "bin").mkdir(parents=True)
        exe = empty_root / "bin" / "mcmctree"
        exe.write_text("#!/bin/sh\nexit 0\n")
        exe.chmod(0o755)
        method = MCMCTreeMethod(
            MCMCTreeConfig(),
            temp_dir,
            software_paths=SoftwarePaths(paml_path=str(empty_root)),
            common_config=CommonConfig(),
        )
        recorder = Mock()
        method.logger = recorder

        method._copy_paml_dat_files()

        assert method._dat_dir is None
        assert method._missing_dat_files == []
        warnings = " ".join(
            " ".join(str(a) for a in c.args) for c in recorder.warning.call_args_list
        )
        assert "dat" in warnings
        assert recorder.error.called is False


# --------------------------------------------------------------------------- #
# 校准写树覆盖（审阅报告 B-6）与根校准 mrca_leaf_pair=None（B-20）
# --------------------------------------------------------------------------- #
class _FakeTree:
    """PhylogeneticTree 的最小替身，只提供 MCMCTree 写树/审计用到的接口。

    不 import ete3（校准标注的真实实现走 ete3，本机解释器上不可用），
    而是按 ``models/tree.py`` 的行为复现：非根校准若缺少 mrca_leaf_pair 就被
    静默 continue，同一 MRCA 节点上的多条校准后者覆盖前者。
    """

    def __init__(self, newick="", tip_names=None, mrca_map=None, drop=(), collide=()):
        self.newick = newick
        self.tip_names = list(tip_names or [])
        self.num_tips = len(self.tip_names)
        self._mrca_map = dict(mrca_map or {})
        self._drop = set(drop)
        self._collide = set(collide)
        self.written_path = None

    # -- 写树路径 ----------------------------------------------------------
    def with_calibration_annotations(self, calibrations):
        parts = []
        seen_mrca = {}
        for cal in calibrations:
            if cal.is_root_node:
                continue
            if not cal.mrca_leaf_pair or cal.name in self._drop:
                continue  # 与 models/tree.py 一样：静默丢弃
            key = frozenset(
                self._mrca_map.get(tuple(cal.mrca_leaf_pair), cal.mrca_leaf_pair)
            )
            if key in seen_mrca:
                parts.remove(seen_mrca[key])  # 与 tree.py 一样：后者覆盖前者
            rendered = cal.age_constraint.to_mcmctree_calib_string()
            entry = f"'{rendered}'" if rendered else ""
            seen_mrca[key] = entry
            parts.append(entry)
        leaves = ",".join(self.tip_names[:2] or ["A", "B"])
        rest = ",".join(self.tip_names[2:] or ["C", "D"])
        body = f"(({leaves}){parts[0] if parts else ''},({rest}){parts[1] if len(parts) > 1 else ''});"
        return _FakeTree(body, tip_names=self.tip_names, mrca_map=self._mrca_map)

    def with_paml48_calibration_annotations(self, calibrations):
        inner = self.with_calibration_annotations(calibrations)
        return _FakeTree(f"{self.num_tips} 1\n{inner.newick}", tip_names=self.tip_names)

    def with_renamed_leaves(self, mapping):
        nk = self.newick
        for old, new in mapping.items():
            nk = nk.replace(old, new)
        return _FakeTree(nk, tip_names=self.tip_names, mrca_map=self._mrca_map)

    def get_mrca(self, tip_names):
        return self._mrca_map.get(tuple(tip_names))

    def write(self, path):
        self.written_path = Path(path)
        self.written_path.write_text(self.newick, encoding="utf-8")


def _cal(name, pair, constraint=None, is_root=False):
    return CalibrationPoint(
        name=name,
        age_constraint=constraint or UniformAgeConstraint(min_age=100.0, max_age=200.0),
        mrca_leaf_pair=pair,
        is_root_node=is_root,
    )


class TestMCMCTreeCalibrationCoverage:
    """B-6：校准被静默丢弃 / 同节点互相覆盖时必须报错，不再静默少跑校准。"""

    def _method(self, temp_dir):
        return MCMCTreeMethod(MCMCTreeConfig(), temp_dir, common_config=CommonConfig())

    def test_missing_leaf_pair_is_a_validation_error(self, temp_dir):
        method = self._method(temp_dir)
        calibrations = [_cal("AB", None), _cal("CD", ("C", "D"))]

        errors = method.validate_calibrations(calibrations)

        assert any("mrca_leaf_pair" in e for e in errors)
        assert any("'AB'" in e for e in errors)

    def test_root_calibration_without_leaf_pair_is_legal(self, temp_dir):
        """B-20：根校准的 mrca_leaf_pair=None 是设计允许的，不得崩溃/报错。"""
        method = self._method(temp_dir)
        root = _cal(
            "LUCA",
            None,
            constraint=UniformAgeConstraint(min_age=3000.0, max_age=4000.0),
            is_root=True,
        )

        errors = method.validate_calibrations([root])

        assert errors == []
        tree = _FakeTree(tip_names=["A", "B", "C", "D"])
        method._audit_calibration_coverage(tree, [root])  # 不抛异常即通过

    def test_identical_leaf_pairs_are_flagged(self, temp_dir):
        method = self._method(temp_dir)
        calibrations = [
            _cal("AB_1", ("A", "B")),
            _cal("AB_2", ("A", "B")),
        ]

        errors = method.validate_calibrations(calibrations)

        assert any("同一 MRCA 叶对" in e for e in errors)

    def test_same_mrca_node_different_pairs_detected_from_topology(self, temp_dir):
        """嵌套分类群：叶对不同但 MRCA 相同 → 同样会被覆盖，必须检出。"""
        method = self._method(temp_dir)
        tree = _FakeTree(
            tip_names=["A", "B", "C", "D"],
            mrca_map={("A", "B"): ["A", "B", "C"], ("A", "C"): ["A", "B", "C"]},
        )

        with pytest.raises(CalibrationError):
            method._audit_calibration_coverage(
                tree, [_cal("AB", ("A", "B")), _cal("AC", ("A", "C"))]
            )

    def test_audit_rejects_calibration_without_leaf_pair(self, temp_dir):
        method = self._method(temp_dir)
        tree = _FakeTree(tip_names=["A", "B", "C", "D"])

        with pytest.raises(CalibrationError):
            method._audit_calibration_coverage(tree, [_cal("AB", None)])

    def test_verify_written_detects_dropped_annotation(self, temp_dir):
        method = self._method(temp_dir)
        written = _FakeTree("((A,B)'B(0.10, 0.20, 0.0)',(C,D));")

        with pytest.raises(CalibrationError):
            method._verify_calibrations_written(
                written, [_cal("AB", ("A", "B")), _cal("CD", ("C", "D"))]
            )

    def test_verify_written_detects_overwritten_duplicate(self, temp_dir):
        """两条校准渲染出同一约束串、但树里只出现一次时也要检出。"""
        method = self._method(temp_dir)
        only_once = _cal("AB", ("A", "B")).age_constraint.to_mcmctree_calib_string()
        written = _FakeTree(f"((A,B)'{only_once}',(C,D));")

        with pytest.raises(CalibrationError):
            method._verify_calibrations_written(
                written,
                [
                    _cal("AB", ("A", "B")),
                    _cal(
                        "AB2",
                        (
                            "A",
                            "B",
                        ),
                    ),
                ],
            )

    def test_verify_written_passes_when_all_present(self, temp_dir):
        method = self._method(temp_dir)
        calibrations = [_cal("AB", ("A", "B")), _cal("CD", ("C", "D"))]
        ab = calibrations[0].age_constraint.to_mcmctree_calib_string()
        cd = UniformAgeConstraint(
            min_age=300.0, max_age=400.0
        ).to_mcmctree_calib_string()
        calibrations[1].age_constraint = UniformAgeConstraint(
            min_age=300.0, max_age=400.0
        )
        written = _FakeTree(f"((A,B)'{ab}',(C,D)'{cd}');")

        method._verify_calibrations_written(written, calibrations)  # 不抛异常

    def test_prepare_inputs_fails_loudly_on_dropped_calibration(self, temp_dir):
        """端到端：写树阶段丢掉一条校准时 prepare_inputs 必须 raise。"""
        method = self._method(temp_dir)
        tree = _FakeTree(
            tip_names=["A", "B", "C", "D"],
            mrca_map={("A", "B"): ["A", "B"], ("C", "D"): ["C", "D"]},
            drop={"CD"},  # 模拟 models/tree.py 静默丢弃
        )

        with pytest.raises(CalibrationError):
            method.prepare_inputs(
                tree, [_cal("AB", ("A", "B")), _cal("CD", ("C", "D"))]
            )

    def test_prepare_inputs_accepts_complete_calibration_set(self, temp_dir):
        """正常路径（含一条根校准，mrca_leaf_pair=None）必须跑通。"""
        method = self._method(temp_dir)
        tree = _FakeTree(
            tip_names=["A", "B", "C", "D"],
            mrca_map={("A", "B"): ["A", "B"], ("C", "D"): ["C", "D"]},
        )
        calibrations = [
            _cal("AB", ("A", "B")),
            _cal("CD", ("C", "D")),
            _cal("LUCA", None, is_root=True),
        ]

        inputs = method.prepare_inputs(tree, calibrations)

        assert Path(inputs["tree"]).exists()
        assert method._has_root_calibration is True
