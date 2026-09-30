"""
MCMCTree 阶段 0 速率估算单元测试

覆盖：
- skip_rate_estimation 配置项（默认 False=不跳过）
- "Substitution rate is per time unit" 提取（4.8 紧凑格式与 4.10.8 空行格式）
- 根年龄解析（根校准约束 / root_age 配置 / 均无）
- 速率估算树生成（去校准标注、根部追加 @age）
- 完整阶段流程（用假 baseml 可执行文件验证 rgene_gamma 由估算速率推导）
- 估算失败时的默认先验回退
"""

import shutil
import tempfile
from pathlib import Path
from unittest.mock import Mock

import pytest

from phylodater.adapters.mcmctree_method import MCMCTreeMethod
from phylodater.infrastructure.configuration import CommonConfig, MCMCTreeConfig
from phylodater.models import (
    CalibrationPoint,
    FixedAgeConstraint,
    MaximumAgeConstraint,
    SoftBoundsConstraint,
    UniformAgeConstraint,
)


@pytest.fixture
def temp_dir():
    temp_path = Path(tempfile.mkdtemp(prefix="mcmctree_rate_test_"))
    yield temp_path
    shutil.rmtree(temp_path, ignore_errors=True)


@pytest.fixture
def config():
    return MCMCTreeConfig(clock=2, burnin=1000, nsample=2000, num_runs=1)


@pytest.fixture
def method(temp_dir, config):
    m = MCMCTreeMethod(
        config, output_dir=temp_dir, common_config=CommonConfig(nthreads=4)
    )
    m._calibrations = []
    return m


def _make_root_cal(constraint):
    cal = Mock(spec=CalibrationPoint)
    cal.name = "Root"
    cal.is_root_node = True
    cal.age_constraint = constraint
    return cal


class TestSkipRateEstimationConfig:
    def test_default_is_not_skipped(self):
        assert MCMCTreeConfig().skip_rate_estimation is False

    def test_should_skip_respects_config(self, temp_dir):
        cfg = MCMCTreeConfig(skip_rate_estimation=True)
        m = MCMCTreeMethod(cfg, output_dir=temp_dir)
        assert m._should_skip_rate_estimation() is True

        cfg2 = MCMCTreeConfig(skip_rate_estimation=False)
        m2 = MCMCTreeMethod(cfg2, output_dir=temp_dir)
        assert m2._should_skip_rate_estimation() is False

    def test_config_coerces_cli_string(self):
        # --method-args "skip_rate_estimation=true" 传入字符串应转为 bool
        cfg = MCMCTreeConfig(skip_rate_estimation="true")
        assert cfg.skip_rate_estimation is True
        cfg = MCMCTreeConfig(skip_rate_estimation="false")
        assert cfg.skip_rate_estimation is False


class TestExtractSubstitutionRate:
    def test_paml48_format(self, method, temp_dir):
        # PAML 4.8: 数值紧跟标题行
        f = temp_dir / "mlb"
        f.write_text("Substitution rate is per time unit\n    4.695519\n")
        assert method._extract_substitution_rate(f) == pytest.approx(4.695519)

    def test_paml4108_format_blank_line(self, method, temp_dir):
        # PAML 4.10.8: 标题与数值之间有一个空行
        f = temp_dir / "mlb"
        f.write_text("Substitution rate is per time unit\n\n    4.695467\n")
        assert method._extract_substitution_rate(f) == pytest.approx(4.695467)

    def test_scientific_notation(self, method, temp_dir):
        f = temp_dir / "mlb"
        f.write_text("Substitution rate is per time unit\n  1.234e-3\n")
        assert method._extract_substitution_rate(f) == pytest.approx(1.234e-3)

    def test_missing_line_returns_none(self, method, temp_dir):
        f = temp_dir / "mlb"
        f.write_text("no rate info here\n")
        assert method._extract_substitution_rate(f) is None

    def test_missing_file_returns_none(self, method, temp_dir):
        assert method._extract_substitution_rate(temp_dir / "nonexistent") is None


class TestResolveRootAge:
    def test_from_fixed_root_calibration(self, method):
        method._calibrations = [_make_root_cal(FixedAgeConstraint(fixed_age=100.0))]
        assert method._resolve_root_age_ga_for_rate_estimation() == pytest.approx(0.1)

    def test_from_uniform_root_calibration(self, method):
        method._calibrations = [
            _make_root_cal(UniformAgeConstraint(min_age=600.0, max_age=800.0))
        ]
        # 中点 700 Ma -> 0.7 Ga
        assert method._resolve_root_age_ga_for_rate_estimation() == pytest.approx(0.7)

    def test_from_soft_bounds_root_calibration(self, method):
        method._calibrations = [
            _make_root_cal(SoftBoundsConstraint(min_age=400.0, max_age=500.0))
        ]
        assert method._resolve_root_age_ga_for_rate_estimation() == pytest.approx(0.45)

    def test_from_maximum_root_calibration(self, method):
        method._calibrations = [_make_root_cal(MaximumAgeConstraint(max_age=900.0))]
        assert method._resolve_root_age_ga_for_rate_estimation() == pytest.approx(0.9)

    def test_fallback_to_config_root_age(self, method, config):
        config.root_age = 50.0  # Ma
        assert method._resolve_root_age_ga_for_rate_estimation() == pytest.approx(0.05)

    def test_root_calibration_wins_over_config(self, method, config):
        config.root_age = 50.0
        method._calibrations = [_make_root_cal(FixedAgeConstraint(fixed_age=100.0))]
        assert method._resolve_root_age_ga_for_rate_estimation() == pytest.approx(0.1)

    def test_none_when_no_source(self, method):
        assert method._resolve_root_age_ga_for_rate_estimation() is None


class TestPrepareRateEstimationTree:
    def test_strips_calibrations_appends_root_age_keeps_header(self, method, temp_dir):
        src = temp_dir / "in.tree"
        src.write_text(
            "4 1\n(((human,chimp)'B(0.006, 0.008, 0.0)'),(mouse,rat)'U(0.012)')'<0.1';\n"
        )
        dest = temp_dir / "rate_est.tree"
        assert method._prepare_rate_estimation_tree(src, dest, 0.1) is True
        content = dest.read_text()
        lines = content.strip().split("\n")
        # baseml/codeml 要求 "n_species n_trees" 头部，必须保留
        assert lines[0].strip() == "4 1"
        body = "\n".join(lines[1:])
        assert body.startswith("(")
        assert "B(" not in body
        assert "U(" not in body
        assert "'" not in body
        assert body.strip().endswith("@0.1;")

    def test_rejects_unparseable_tree(self, method, temp_dir):
        src = temp_dir / "in.tree"
        src.write_text("not a tree at all")
        dest = temp_dir / "rate_est.tree"
        assert method._prepare_rate_estimation_tree(src, dest, 0.1) is False


class TestGenerateControlFileUsesEstimatedBeta:
    def test_4108_uses_calculated_beta(self, method, temp_dir, config):
        config.paml_version = "4.10.8"
        method._calculated_beta = 3.14
        ctl = temp_dir / "mcmctree.ctl"
        method._generate_control_file(ctl, usedata=2)
        assert "rgene_gamma  = 2.0000 3.1400" in ctl.read_text()

    def test_48_uses_calculated_beta(self, method, temp_dir, config):
        config.paml_version = "4.8"
        method._calculated_beta = 3.14
        ctl = temp_dir / "mcmctree.ctl"
        method._generate_control_file(ctl, usedata=2)
        assert "rgene_gamma  = 2.0000 3.1400" in ctl.read_text()

    def test_default_beta_when_no_estimation(self, method, temp_dir, config):
        # 跳过速率估算时（_calculated_beta 为 None），4.10.8 回退默认
        # beta = rate_alpha / base_rate(0.1) = 20
        config.paml_version = "4.10.8"
        method._calculated_beta = None
        ctl = temp_dir / "mcmctree.ctl"
        method._generate_control_file(ctl, usedata=2)
        assert "rgene_gamma  = 2.0000 20.0000" in ctl.read_text()


class TestRateEstimationStage:
    def _setup_inputs(self, method, temp_dir):
        seq = temp_dir / "mcmctree.phy"
        seq.write_text("4 100\nhuman  ACGT...\n")
        tree = temp_dir / "mcmctree.tree"
        tree.write_text("4 1\n(((human,chimp)'B(0.006, 0.008)'),(mouse,rat));\n")
        method._input_files = {"tree": tree, "seq": seq}
        method._seqtype = 0  # DNA

    def test_full_stage_with_fake_baseml(self, method, temp_dir, config):
        """假 baseml 返回速率行 → rgene_gamma 由估算速率推导"""
        self._setup_inputs(method, temp_dir)
        method._calibrations = [_make_root_cal(FixedAgeConstraint(fixed_age=100.0))]

        fake_bin = temp_dir / "fakebin"
        fake_bin.mkdir(exist_ok=True)
        rate_value = "0.084825"
        script = fake_bin / "baseml"
        script.write_text(
            "#!/bin/sh\n"
            'printf "Substitution rate is per time unit\\n  '
            f"{rate_value}\\n"
            '" > mlb\nexit 0\n'
        )
        script.chmod(0o755)

        orig = method._get_executable
        method._get_executable = lambda name, _orig=orig: (
            str(script) if name == "baseml" else _orig(name)
        )

        method._run_rate_estimation_stage()

        # beta = rate_alpha / (rate * rate_scale) = 2 / 0.084825
        assert method._estimated_rate == pytest.approx(0.084825)
        assert method._calculated_beta == pytest.approx(2.0 / 0.084825)

        rate_dir = method.work_dir / "rate_estimation"
        assert (rate_dir / "baseml.ctl").exists()
        assert "clock = 1" in (rate_dir / "baseml.ctl").read_text()
        # 速率估算树: 去校准标注 + 根部 @age (100 Ma -> 0.1 Ga)
        rate_tree = (rate_dir / "rate_est.tree").read_text()
        assert "B(" not in rate_tree
        assert rate_tree.strip().endswith("@0.1;")
        assert (rate_dir / "estimated_rate.txt").exists()

    def test_fallback_when_exe_fails(self, method, temp_dir, config):
        """baseml 运行失败 → 回退默认先验 beta = rate_alpha/20"""
        self._setup_inputs(method, temp_dir)
        method._calibrations = [_make_root_cal(FixedAgeConstraint(fixed_age=100.0))]

        script = temp_dir / "fake_baseml_fail"
        script.write_text("#!/bin/sh\nexit 1\n")
        script.chmod(0o755)

        orig = method._get_executable
        method._get_executable = lambda name, _orig=orig: (
            str(script) if name == "baseml" else _orig(name)
        )

        method._run_rate_estimation_stage()
        assert method._estimated_rate is None
        assert method._calculated_beta == pytest.approx(config.rate_alpha / 20.0)

    def test_fallback_without_root_age(self, method, temp_dir):
        """无根校准且未配置 root_age → 回退默认先验，不运行 baseml"""
        self._setup_inputs(method, temp_dir)
        method._calibrations = []

        def _unexpected(name):
            raise AssertionError(f"baseml should not be invoked, got {name}")

        orig = method._get_executable
        method._get_executable = lambda name, _orig=orig: (
            _unexpected(name) if name == "baseml" else _orig(name)
        )

        method._run_rate_estimation_stage()
        assert method._estimated_rate is None
        assert method._calculated_beta == pytest.approx(method.config.rate_alpha / 20.0)

    def test_codeml_used_for_protein(self, method, temp_dir, config):
        """蛋白质序列应生成 codeml 控制文件并解析 mlc"""
        self._setup_inputs(method, temp_dir)
        method._seqtype = 2  # 蛋白
        method._calibrations = [_make_root_cal(FixedAgeConstraint(fixed_age=100.0))]

        script = temp_dir / "fake_codeml"
        script.write_text(
            "#!/bin/sh\n"
            'printf "Substitution rate is per time unit\\n\\n  1.5\\n" > mlc\n'
            "exit 0\n"
        )
        script.chmod(0o755)

        calls = {}

        def _fake_exec(name):
            calls["name"] = name
            return str(script)

        method._get_executable = _fake_exec

        method._run_rate_estimation_stage()
        assert calls["name"] == "codeml"
        rate_dir = method.work_dir / "rate_estimation"
        assert (rate_dir / "codeml.ctl").exists()
        assert "seqtype = 2" in (rate_dir / "codeml.ctl").read_text()
        assert method._estimated_rate == pytest.approx(1.5)
        assert method._calculated_beta == pytest.approx(2.0 / 1.5)


class TestRootCalibrationWithoutLeafPair:
    """B-20：根校准的 mrca_leaf_pair=None 是合法情形，速率估算路径不得崩溃。"""

    def test_real_root_calibration_object_with_none_pair(self, method, config):
        root_cal = CalibrationPoint(
            name="LUCA",
            age_constraint=UniformAgeConstraint(min_age=3800.0, max_age=4200.0),
            mrca_leaf_pair=None,
            is_root_node=True,
        )
        method._calibrations = [root_cal]

        # 中点 4000 Ma -> 4.0 Ga；不得因 pair[0] 下标而抛 TypeError
        assert method._resolve_root_age_ga_for_rate_estimation() == pytest.approx(4.0)

    def test_non_root_calibrations_are_ignored_for_root_age(self, method, config):
        method._calibrations = [
            CalibrationPoint(
                name="AB",
                age_constraint=UniformAgeConstraint(min_age=100.0, max_age=200.0),
                mrca_leaf_pair=("A", "B"),
                is_root_node=False,
            )
        ]
        config.root_age = None

        assert method._resolve_root_age_ga_for_rate_estimation() is None

    def test_rate_stage_falls_back_loudly_without_root_age(self, method, temp_dir):
        """无根年龄可确定时回退默认先验，并留下 WARNING（不静默换先验）。"""
        seq = temp_dir / "mcmctree.phy"
        seq.write_text("4 100\nhuman  ACGT...\n")
        tree = temp_dir / "mcmctree.tree"
        tree.write_text("4 1\n(((human,chimp)'B(0.006, 0.008)'),(mouse,rat));\n")
        method._input_files = {"tree": tree, "seq": seq}
        method._seqtype = 0
        method._calibrations = []
        method.config.root_age = None

        recorder = Mock()
        method.logger = recorder
        method._run_rate_estimation_stage()

        assert method._estimated_rate is None
        assert method._calculated_beta == pytest.approx(method.config.rate_alpha / 20.0)
        warnings = " ".join(
            " ".join(str(a) for a in c.args) for c in recorder.warning.call_args_list
        )
        assert "No absolute root age available" in warnings
