"""Input validation and injection security tests."""

import pytest

from tests.validation.helpers import run_cli


@pytest.mark.security
class TestInputSafety:
    def test_malicious_unicode_tree_rejected(
        self, fixtures_dir, small_alignment_path, valid_calibration_path, tmp_test_dir
    ):
        """Trees containing bi-directional override characters should not be silently accepted."""
        malicious_tree = fixtures_dir / "trees" / "malicious_tree.nwk"
        result = run_cli(
            [
                "dating",
                "-t",
                str(malicious_tree),
                "-s",
                str(small_alignment_path),
                "-c",
                str(valid_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--dry-run",
            ]
        )
        # The CLI may either reject it outright or fail validation; either is safe.
        assert result.returncode != 0, result.stderr

    def test_command_injection_in_method_args(
        self,
        small_tree_path,
        small_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        """Method arguments must not be interpreted as shell metacharacters."""
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_alignment_path),
                "-c",
                str(small_tree_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--method-args",
                "burnin=100; nsample=200 | rm -rf /",
                "--dry-run",
            ]
        )
        # Should fail because the argument is malformed, but must not execute commands.
        assert result.returncode != 0, result.stderr

    def test_negative_branch_length_rejected(
        self, fixtures_dir, small_alignment_path, valid_calibration_path, tmp_test_dir
    ):
        """Trees with negative branch lengths should be rejected."""
        bad_tree = fixtures_dir / "trees" / "negative_branch.nwk"
        result = run_cli(
            [
                "dating",
                "-t",
                str(bad_tree),
                "-s",
                str(small_alignment_path),
                "-c",
                str(valid_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--dry-run",
            ]
        )
        assert result.returncode != 0, result.stderr

    def test_calibration_with_nan_rejected(
        self, fixtures_dir, small_tree_path, small_alignment_path, tmp_test_dir
    ):
        """Calibrations containing NaN ages must be rejected."""
        nan_cal = fixtures_dir / "calibrations" / "nan_age.yaml"
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_alignment_path),
                "-c",
                str(nan_cal),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--dry-run",
            ]
        )
        assert result.returncode != 0, result.stderr

    def test_missing_calibration_type_rejected(
        self, fixtures_dir, small_tree_path, small_alignment_path, tmp_test_dir
    ):
        """Calibrations with missing constraint type must be rejected."""
        missing_type = fixtures_dir / "calibrations" / "missing_type.yaml"
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_alignment_path),
                "-c",
                str(missing_type),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--dry-run",
            ]
        )
        assert result.returncode != 0, result.stderr
