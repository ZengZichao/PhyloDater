"""Functional tests for CLI error handling."""

from tests.validation.helpers import run_cli


class TestCliErrors:
    def test_missing_subcommand(self):
        result = run_cli([])
        assert result.returncode == 2, result.stderr

    def test_dating_nonexistent_tree(self):
        # A non-existent tree path is validated before argparse missing-arg checks.
        result = run_cli(
            [
                "dating",
                "-t",
                "no_such_tree.nwk",
                "-s",
                "seq.fa",
                "-c",
                "cal.yaml",
                "-o",
                "out",
                "--method",
                "mcmctree",
            ]
        )
        assert result.returncode == 3, result.stderr

    def test_dating_nonexistent_sequence(self, small_tree_path):
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                "no_such_seq.fa",
                "-c",
                "cal.yaml",
                "-o",
                "out",
                "--method",
                "mcmctree",
            ]
        )
        assert result.returncode == 3, result.stderr

    def test_dating_nonexistent_calibration(
        self, small_tree_path, small_alignment_path, tmp_test_dir
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_alignment_path),
                "-c",
                "nonexistent.yaml",
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
            ]
        )
        assert result.returncode == 3, result.stderr

    def test_unknown_method(
        self,
        small_tree_path,
        small_alignment_path,
        valid_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_alignment_path),
                "-c",
                str(valid_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "notamethod",
            ]
        )
        assert result.returncode == 2, result.stderr
