"""
MDCatMethod 适配器单元测试

测试 MD-Cat 定年方法的各项功能
"""

import shutil
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.adapters.mdcat_method import MDCatMethod
from phylodater.core.exceptions import CalibrationError
from phylodater.infrastructure.configuration import CommonConfig, MDCatConfig
from phylodater.models import (
    CalibrationPoint,
    FixedAgeConstraint,
    PhylogeneticTree,
    UniformAgeConstraint,
)


class TestMDCatMethod:
    """MDCatMethod 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="mdcat_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def mdcat_config(self):
        """创建 MD-Cat 配置"""
        return MDCatConfig(ncat=50, nrep=10, max_iter=50)

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

    def test_init_creates_work_directory(self, temp_dir, mdcat_config, common_config):
        """测试初始化创建工作目录"""
        method = MDCatMethod(mdcat_config, temp_dir, common_config=common_config)

        assert method.work_dir.exists()
        assert method.method_name == "mdcat"

    def test_method_name_property(self, temp_dir, mdcat_config, common_config):
        """测试方法名称属性"""
        method = MDCatMethod(mdcat_config, temp_dir, common_config=common_config)
        assert method.method_name == "mdcat"

    def test_common_config_usage(self, temp_dir, mdcat_config, common_config):
        """测试 common_config 使用"""
        method = MDCatMethod(mdcat_config, temp_dir, common_config=common_config)
        assert method.common_config.seed == 42
        assert method.common_config.verbose is False

    def test_config_parameters(self, temp_dir, mdcat_config, common_config):
        """测试配置参数"""
        method = MDCatMethod(mdcat_config, temp_dir, common_config=common_config)

        assert method.config.ncat == 50
        assert method.config.nrep == 10
        assert method.config.max_iter == 50

    def test_backward_time_default(self, temp_dir, mdcat_config, common_config):
        """测试向后时间默认值"""
        method = MDCatMethod(mdcat_config, temp_dir, common_config=common_config)
        assert method.config.backward_time is True

    def test_ci_parameters(self, temp_dir, mdcat_config, common_config):
        """测试 CI 参数：bootstrap 默认关闭（与上游 md_cat.py 一致）"""
        method = MDCatMethod(mdcat_config, temp_dir, common_config=common_config)
        assert method.config.ci_nboots == 0
        assert method.config.ci_plower == 0.025
        assert method.config.ci_pupper == 0.975

    def test_work_dir_isolation(self, temp_dir, mdcat_config, common_config):
        """测试工作目录隔离"""
        method1 = MDCatMethod(mdcat_config, temp_dir, common_config=common_config)
        method2 = MDCatMethod(
            mdcat_config, temp_dir / "other", common_config=common_config
        )

        assert method1.work_dir != method2.work_dir


class TestMDCatConfig:
    """MD-Cat 配置测试"""

    def test_default_values(self):
        """测试默认值"""
        config = MDCatConfig()
        assert config.ncat == 50
        assert config.nrep == 100
        assert config.max_iter == 100

    def test_custom_values(self):
        """测试自定义值"""
        config = MDCatConfig(ncat=20, nrep=50, max_iter=200)
        assert config.ncat == 20
        assert config.nrep == 50
        assert config.max_iter == 200

    def test_ncat_validation(self):
        """测试 ncat 验证"""
        with pytest.raises(ValueError):
            MDCatConfig(ncat=0)

    def test_ci_validation(self):
        """测试 CI 参数验证"""
        with pytest.raises(ValueError):
            MDCatConfig(ci_plower=0.975, ci_pupper=0.025)


class TestMDCatCliContract:
    """md-cat 命令行必须匹配上游 md_cat.py 的 argparse（--CI 为单个空白分隔串）。

    上游代码：
        tokens = args["CI"].split()          # 单一 argparse 值，内部按空白切分
        nboots, p_lower, p_upper = tokens[0], tokens[1], tokens[2]
    因此三个裸 token 会被 argparse 拒（unrecognized arguments），逗号串会让
    ``int()`` 报 ValueError；唯一合法形式是 ``--CI "100 0.025 0.975"``。
    """

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="mdcat_cli_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def _prep(self, temp_dir):
        m = MDCatMethod(
            MDCatConfig(ncat=50, ci_nboots=100),
            temp_dir,
            common_config=CommonConfig(nthreads=2),
        )
        m._use_direct_import = False
        wd = m.work_dir
        tree_f = wd / "t.nwk"
        tree_f.write_text("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")
        cons = wd / "c.txt"
        cons.write_text("Primates=human+chimp\t7.0\n")
        m._input_files = {"tree": tree_f, "constraints": cons, "seq_len": 1000}
        return m

    def _run_capture(self, m, monkeypatch_rec):
        def fake_run(cmd):
            monkeypatch_rec["cmd"] = list(cmd)
            outp = Path(cmd[cmd.index("-o") + 1])
            outp.write_text("((human:1.0,chimp:1.0)X:1.0);")
            r = Mock()
            r.returncode = 0
            r.stderr = ""
            r.stdout = ""
            return r

        with patch("phylodater.adapters.mdcat_method.ProcessRunner") as PR:
            PR.return_value.run.side_effect = fake_run
            assert m._execute_subprocess() is True

    def test_ci_is_one_whitespace_separated_token(self, temp_dir):
        rec = {}
        m = self._prep(temp_dir)
        self._run_capture(m, rec)
        cmd = rec["cmd"]
        assert "--CI" in cmd
        i = cmd.index("--CI")
        ci_val = cmd[i + 1]
        # 上游只吃一个 token：形如 "n l u"，两个空白，且后面紧跟的是开关而非裸值
        assert ci_val.count(",") == 0
        assert ci_val.split() == ["100", "0.025", "0.975"]
        assert cmd[i + 2].startswith("-")

    def test_flags_match_upstream_argparse(self, temp_dir):
        rec = {}
        m = self._prep(temp_dir)
        self._run_capture(m, rec)
        cmd = rec["cmd"]
        for flag in [
            "-i",
            "-t",
            "-l",
            "-k",
            "-p",
            "--maxIter",
            "-o",
            "-b",
            "--rootTime",
            "--leafTime",
            "--annotate",
            "--CI",
        ]:
            assert flag in cmd, f"missing flag {flag}"
        # 除 -b 外，每个带值选项后必须紧跟一个值（不能有游离位置参数把 --CI 撑爆）
        assert "-b" in cmd  # backward_time 默认 True


class TestMDCatDegradationDisclosure:
    """B-27：拓扑匹配失败时的降级披露不得因把异常对象喂给 logger 而崩溃。"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="mdcat_disc_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_label_fallback_logs_string_not_object(self, temp_dir):
        from phylodater.models import CalibrationPoint

        m = MDCatMethod(MDCatConfig(), temp_dir, common_config=CommonConfig(nthreads=2))
        cal = CalibrationPoint(
            name="Primates",
            age_constraint=FixedAgeConstraint(fixed_age=7.0),
            resolved_taxa=["human", "chimp"],
            mrca_leaf_pair=("human", "chimp"),
            is_root_node=False,
        )
        m._calibrations = [cal]
        # 带 [t=...] 注解、且 ete3 在本环境不可导入 -> 必然进入标签匹配披露分支
        annot = (
            "((human:0.1,chimp:0.1)Primates[t=7.5,mu=0.1]:0.2,"
            "mouse:0.3,rat:0.3)Root[t=1.0,mu=0.05]:0.0;"
        )
        with patch.object(m.logger, "warning") as wmock:
            ages, warns = m._parse_tree_ages(annot)  # 不得抛 TypeError
        for call in wmock.call_args_list:
            assert isinstance(call.args[0], str), "logger.warning 必须收到字符串"
        assert any(isinstance(w, Exception) for w in warns)


class TestMDCatCiSemantics:
    """C-52：区间语义由配置分位数推导，而不是无条件写死 CI95。"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="mdcat_ci_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_default_quantiles_are_ci95(self, temp_dir):
        from phylodater.models import CIType

        m = MDCatMethod(MDCatConfig(), temp_dir, common_config=CommonConfig(nthreads=2))
        assert m._resolved_ci_type() == CIType.CI95

    def test_custom_quantiles_are_range(self, temp_dir):
        from phylodater.models import CIType

        m = MDCatMethod(
            MDCatConfig(ci_plower=0.05, ci_pupper=0.95),
            temp_dir,
            common_config=CommonConfig(nthreads=2),
        )
        assert m._resolved_ci_type() == CIType.RANGE

    def test_ci_numeric_regex_accepts_scientific(self, temp_dir):

        m = MDCatMethod(MDCatConfig(), temp_dir, common_config=CommonConfig(nthreads=2))
        m._calibrations = []
        annot = (
            "((human:0.1,chimp:0.1)Primates[t=4e-02,mu=1e-03,"
            "CI={3e-02,5e-02}]:0.2,mouse:0.3,rat:0.3)Root:0.0;"
        )
        # 不应因旧 [\\d.]+ 漏掉科学计数法而完全解析不到（这里只断言不崩溃、且无 CI 崩溃）
        ages, warns = m._parse_tree_ages(annot)
        assert isinstance(ages, dict)


class TestMDCatCalibrationGuards:
    """C-22 / C-58：约束文件的兜底检查与直接导入回退的收窄。"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="mdcat_guard_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_prepare_inputs_raises_when_nothing_written(self, temp_dir):
        from phylodater.models import (
            CalibrationPoint,
            FixedAgeConstraint,
            PhylogeneticTree,
        )

        m = MDCatMethod(MDCatConfig(), temp_dir, common_config=CommonConfig(nthreads=2))
        tree = Mock(spec=PhylogeneticTree)
        tree.newick = "((human:0.1,chimp:0.1):0.2);"
        tree.write = Mock()
        bad = CalibrationPoint(
            name="NoTaxa",
            age_constraint=FixedAgeConstraint(fixed_age=7.0),
            resolved_taxa=[],
            mrca_leaf_pair=None,
            is_root_node=False,
        )
        with pytest.raises(CalibrationError):
            m.prepare_inputs(tree, [bad], alignment_path=None)

    def test_direct_import_re_raises_signature_errors(self, temp_dir, monkeypatch):
        import sys

        m = MDCatMethod(MDCatConfig(), temp_dir, common_config=CommonConfig(nthreads=2))
        # 让 _execute_direct 内的 `import treeswift` 成功，从而真正触及 MDCat 调用
        monkeypatch.setitem(sys.modules, "treeswift", Mock())
        m._use_direct_import = True
        m._mdcat_func = Mock(side_effect=TypeError("emd API signature changed"))
        m._input_files = {
            "tree": temp_dir / "t.nwk",
            "constraints": temp_dir / "c.txt",
            "seq_len": 1000,
        }
        m._execute_subprocess = Mock(return_value=True)
        # C-58：签名漂移(TypeError)应上抛，而不是被静默吞掉后改跑 CLI
        with pytest.raises(TypeError):
            m._execute_direct()
        m._execute_subprocess.assert_not_called()
