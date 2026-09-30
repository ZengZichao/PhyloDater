"""Unit tests for CalibrationLoader."""

import pytest

from phylodater.core.exceptions import CalibrationError
from phylodater.services import CalibrationLoader


class TestCalibrationLoader:
    def test_load_valid_uniform(self, valid_calibration_path):
        loader = CalibrationLoader()
        cals = loader.load(valid_calibration_path)
        assert len(cals) == 1
        assert cals[0].name == "Primates"

    def test_load_valid_fixed(self, fixed_calibration_path):
        loader = CalibrationLoader()
        cals = loader.load(fixed_calibration_path)
        assert len(cals) == 1
        assert cals[0].name == "Rodents"

    def test_load_missing_type(self, tmp_test_dir):
        from tests.validation.helpers import make_calibration_file

        path = tmp_test_dir / "missing_type.yaml"
        make_calibration_file(
            path,
            [{"name": "Bad", "mrca_pair": ["A", "B"], "constraint": {"age": 10.0}}],
        )
        loader = CalibrationLoader()
        cals = loader.load(path)
        assert len(cals) == 0

    def test_load_nan_age(self, tmp_test_dir):
        from tests.validation.helpers import make_calibration_file

        path = tmp_test_dir / "nan.yaml"
        make_calibration_file(
            path,
            [
                {
                    "name": "Bad",
                    "mrca_pair": ["A", "B"],
                    "constraint": {"type": "fixed", "age": float("nan")},
                }
            ],
        )
        loader = CalibrationLoader()
        cals = loader.load(path)
        assert len(cals) == 0

    def test_load_empty_calibrations(self, empty_calibration_path):
        loader = CalibrationLoader()
        cals = loader.load(empty_calibration_path)
        assert len(cals) == 0

    def test_strict_mode_raises(self, tmp_test_dir):
        from tests.validation.helpers import make_calibration_file

        path = tmp_test_dir / "bad.yaml"
        make_calibration_file(
            path,
            [{"name": "Bad", "mrca_pair": ["A", "B"], "constraint": {"age": 10.0}}],
        )
        loader = CalibrationLoader()
        with pytest.raises(CalibrationError):
            loader.load(path, strict=True)

    def test_load_all_types(self, fixtures_dir):
        loader = CalibrationLoader()
        cals = loader.load(fixtures_dir / "calibrations" / "all_types.yaml")
        assert len(cals) == 5
        names = {c.name for c in cals}
        assert names == {"Fixed", "Uniform", "SoftLower", "Maximum", "SoftBounds"}
