"""Functional tests for the 'phylodater fold' CLI."""

from tests.validation.helpers import run_cli


class TestCliFold:
    def test_fold_not_implemented(self, small_tree_path, tmp_test_dir):
        result = run_cli(
            [
                "fold",
                "-t",
                str(small_tree_path),
                "-o",
                str(tmp_test_dir),
                "--fold-rank",
                "family",
            ]
        )
        # fold is experimental and not yet implemented -> runtime error
        assert result.returncode == 1, result.stderr
