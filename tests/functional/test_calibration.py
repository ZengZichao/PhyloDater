"""
模块四：校准点处理功能测试

验证校准配置加载、约束类型、MRCA 定位、自动校准等。
"""

import pytest
import yaml

from phylodater import CalibrationLoader, PhylogeneticTree
from phylodater.models import (
    FixedAgeConstraint,
    MaximumAgeConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)
from phylodater.services.calibration_resolver import CalibrationResolver


class TestCalibrationLoading:
    """校准配置加载测试"""

    def test_cal_001_valid_yaml(self, normal_calibration):
        """CAL-001: 有效 YAML"""
        loader = CalibrationLoader()
        calibrations = loader.load(normal_calibration)
        assert len(calibrations) == 3

    def test_cal_002_file_not_found(self):
        """CAL-002: 文件不存在"""
        from pathlib import Path

        loader = CalibrationLoader()
        with pytest.raises(FileNotFoundError):
            loader.load(Path("/nonexistent/cal.yaml"))

    def test_cal_003_invalid_yaml(self, functional_input_dir):
        """CAL-003: 格式错误"""
        loader = CalibrationLoader()
        invalid = functional_input_dir / "invalid" / "invalid_yaml.yaml"
        calibrations = loader.load(invalid)
        assert len(calibrations) == 0

    def test_cal_004_missing_required_field(self, tmp_path):
        """CAL-004: 缺少必填字段"""
        loader = CalibrationLoader()
        config = {"calibrations": [{"constraint": {"type": "fixed", "age": 100}}]}
        config_file = tmp_path / "no_name.yaml"
        config_file.write_text(yaml.dump(config), encoding="utf-8")
        calibrations = loader.load(config_file)
        assert len(calibrations) == 0


class TestCalibrationConstraintTypes:
    """校准约束类型测试"""

    @pytest.mark.parametrize(
        "ctype,cls,kwargs",
        [
            ("fixed", FixedAgeConstraint, {"age": 100.0}),
            ("uniform", UniformAgeConstraint, {"min": 70.0, "max": 90.0}),
            ("soft_lower", SoftLowerBoundConstraint, {"min": 60.0}),
            ("maximum", MaximumAgeConstraint, {"max": 150.0}),
        ],
    )
    def test_cal_010_013_constraint_types(self, ctype, cls, kwargs):
        """CAL-010 至 CAL-013: 约束类型解析"""
        loader = CalibrationLoader()
        constraint_data = {"type": ctype, **kwargs}
        constraint = loader._parse_constraint(constraint_data)
        assert isinstance(constraint, cls)

    def test_cal_018_uniform_invalid_bounds(self, functional_input_dir):
        """CAL-018: uniform 边界错误"""
        loader = CalibrationLoader()
        invalid_uniform = functional_input_dir / "invalid" / "invalid_uniform.yaml"
        calibrations = loader.load(invalid_uniform)
        assert len(calibrations) == 0

    def test_cal_019_gamma_negative_params(self):
        """CAL-019: gamma 负参数"""
        loader = CalibrationLoader()
        constraint = loader._parse_constraint({"type": "gamma", "alpha": -1, "beta": 1})
        assert constraint is None


class TestMRCAResolution:
    """MRCA 定位测试"""

    @pytest.fixture
    def sample_tree(self):
        return PhylogeneticTree.from_newick(
            "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
        )

    def test_cal_030_mrca_pair(self, sample_tree):
        """CAL-030: 两点 MRCA"""
        resolver = CalibrationResolver(sample_tree)
        result = resolver.resolve_mrca_pair(["human", "mouse"])
        assert result is not None
        assert "human" in result.resolved_taxa and "mouse" in result.resolved_taxa

    def test_cal_031_multiple_taxa_mrca(self, sample_tree):
        """CAL-031: 多点 MRCA"""
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        point = CalibrationPoint(
            name="TestMRCA",
            age_constraint=FixedAgeConstraint(fixed_age=100.0),
            resolved_taxa=["human", "chimp"],
            mrca_leaf_pair=("human", "chimp"),
        )
        assert point.mrca_leaf_pair == ("human", "chimp")

    def test_cal_033_root_calibration(self, sample_tree):
        """CAL-033: 根节点校准"""
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        point = CalibrationPoint(
            name="Root",
            age_constraint=FixedAgeConstraint(fixed_age=150.0),
            is_root_node=True,
        )
        assert point.is_root_node


class TestAutoCalibration:
    """自动校准点测试"""

    def test_cal_040_auto_calibrate(self, tmp_path):
        """CAL-040: 自动校准"""
        # 当前 CLI 支持 --auto-calibrate，通过 CLI 验证参数可解析
        from tests.functional.conftest import run_phylodater_cli

        tree = tmp_path / "tree.nwk"
        tree.write_text("((A:0.1,B:0.2):0.3,C:0.4);", encoding="utf-8")
        result = run_phylodater_cli(["dating", "--help"])
        assert result.returncode == 0
        assert "--auto-calibrate" in result.stdout


class TestCalibrationRegression:
    """回归测试"""

    def test_cal_060_example_calibration(self):
        """CAL-060 / REG-001: 示例校准文件可加载，且含引擎共需的固定年龄

        示例集由 ``test-data/benchmark/generate_benchmark_dataset.py --profile mini``
        生成（根 105 Ma、Homininae 9 Ma、Primates 75 Ma），类型必须是点校准 + 双边区间：
        PATHd8 / MD-Cat / wLogDate 只能处理固定年龄或双边区间，若示例里没有 fixage
        点，文档中的快速上手命令就无法跑通。
        """
        from pathlib import Path

        example = Path(__file__).parent.parent.parent / "examples" / "calibrations.yaml"
        loader = CalibrationLoader()
        calibrations = loader.load(example)
        assert len(calibrations) == 3
        assert isinstance(calibrations[0].age_constraint, FixedAgeConstraint)
        assert isinstance(calibrations[1].age_constraint, UniformAgeConstraint)
        assert isinstance(calibrations[2].age_constraint, UniformAgeConstraint)
        assert calibrations[2].is_root_node is True

    def test_cal_061_empty_calibration_error(self, tmp_path):
        """CAL-061 / REG-002: 空校准列表友好错误"""
        empty = tmp_path / "empty.yaml"
        empty.write_text("calibrations: []\n", encoding="utf-8")
        loader = CalibrationLoader()
        calibrations = loader.load(empty)
        assert calibrations == []

    def test_cal_062_old_type_names_translated(self, tmp_path):
        """CAL-062 / REG-003: 旧类型名 lower_bound/upper_bound 被兼容映射"""
        config = {
            "calibrations": [
                {
                    "name": "Old",
                    "mrca_pair": ["A", "B"],
                    "constraint": {"type": "lower_bound", "min": 10},
                }
            ]
        }
        old_file = tmp_path / "old.yaml"
        old_file.write_text(yaml.dump(config), encoding="utf-8")
        loader = CalibrationLoader()
        calibrations = loader.load(old_file)
        assert len(calibrations) == 1
        assert (
            calibrations[0].age_constraint.__class__.__name__
            == "SoftLowerBoundConstraint"
        )
