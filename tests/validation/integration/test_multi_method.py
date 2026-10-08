"""Integration test running multiple methods."""

import pytest

from tests.validation.helpers import run_cli


@pytest.mark.slow
@pytest.mark.requires_mcmctree
class TestMultiMethod:
    """Runs the real MCMCTree engine end to end through the CLI.

    ``requires_mcmctree`` is required in addition to ``slow``: this class shells
    out to ``phylodater dating --method mcmctree``, so without PAML on PATH the
    CLI legitimately exits 1 and the test reports a *failure* rather than an
    unavailable engine. The sibling ``test_end_to_end_mcmctree.py`` already
    carries both markers; ``pytest_runtest_setup`` in
    ``tests/validation/conftest.py`` turns the marker into a skip.
    """
    def test_mcmctree_only(
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
