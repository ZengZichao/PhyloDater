"""End-to-end integration test for MCMCTree."""

import pytest

from tests.validation.helpers import run_cli


@pytest.mark.slow
@pytest.mark.requires_mcmctree
class TestEndToEndMCMCTree:
    def test_mcmctree_dry_run(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
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
                "--dry-run",
            ]
        )
        assert result.returncode == 0, result.stderr

    def test_mcmctree_real_run(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
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
        assert (
            (tmp_test_dir / "mcmctree_dated_tree.nhx").exists()
            or (tmp_test_dir / "mcmctree_dated_tree.nwk").exists()
            or (tmp_test_dir / "mcmctree_dated_tree.nexus").exists()
        )
