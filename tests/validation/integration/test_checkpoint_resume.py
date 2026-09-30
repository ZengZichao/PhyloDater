"""Integration tests for checkpoint resume behavior."""

import pytest

from tests.validation.helpers import run_cli


@pytest.mark.slow
@pytest.mark.requires_mcmctree
class TestCheckpointResume:
    def test_resume_after_success(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        # First run
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_tree_alignment_path),
                "-c",
                str(small_tree_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--method-args",
                "burnin=100,nsample=200,paml_version=4.10.10",
                "--threads",
                "1",
            ]
        )
        assert result.returncode == 0, result.stderr

        # Resume with identical input should succeed quickly
        result2 = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_tree_alignment_path),
                "-c",
                str(small_tree_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--method-args",
                "burnin=100,nsample=200,paml_version=4.10.10",
                "--threads",
                "1",
                "--resume",
                "--no-clobber",
            ]
        )
        assert result2.returncode == 0, result2.stderr
