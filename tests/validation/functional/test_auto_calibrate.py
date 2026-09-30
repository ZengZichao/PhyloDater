"""Functional tests for --auto-calibrate."""

from tests.validation.helpers import run_cli


class TestAutoCalibrate:
    def test_auto_calibrate_cli(
        self,
        taxonomy_tree_path,
        taxonomy_tree_alignment_path,
        taxonomy_tree_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(taxonomy_tree_path),
                "-s",
                str(taxonomy_tree_alignment_path),
                "-c",
                str(taxonomy_tree_calibration_path),
                "-o",
                str(tmp_test_dir / "out"),
                "--method",
                "mcmctree",
                "--dry-run",
                "--auto-calibrate",
                "Bact:50",
            ]
        )
        assert result.returncode == 0, result.stderr
