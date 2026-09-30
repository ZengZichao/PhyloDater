"""YAML parsing and calibration-file safety tests."""

import pytest

from tests.validation.helpers import run_cli, write_text


@pytest.mark.security
class TestYAMLSafety:
    def test_calibration_file_with_anchors_rejected(
        self, small_tree_path, small_alignment_path, tmp_test_dir
    ):
        """YAML aliases/anchors should not cause unexpected object expansion."""
        cal = tmp_test_dir / "anchors.yaml"
        write_text(
            cal,
            "calibrations:\n"
            "  - &cal\n"
            "    name: A\n"
            "    mrca_pair: [A, B]\n"
            "    constraint: &c {type: uniform, min: 10, max: 20}\n"
            "  - *cal\n",
        )
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_alignment_path),
                "-c",
                str(cal),
                "-o",
                str(tmp_test_dir / "out"),
                "--method",
                "mcmctree",
                "--dry-run",
            ]
        )
        # Duplicate-name calibrations should be rejected or fail validation.
        assert result.returncode != 0, result.stderr

    def test_calibration_file_with_exe_tag_rejected(
        self, small_tree_path, small_alignment_path, tmp_test_dir
    ):
        """Custom YAML tags that could construct arbitrary objects must be rejected."""
        cal = tmp_test_dir / "tag.yaml"
        write_text(
            cal,
            "calibrations:\n"
            "  - name: A\n"
            "    mrca_pair: [A, B]\n"
            "    constraint: !python/object/apply:os.system ['echo pwned']\n",
        )
        result = run_cli(
            [
                "dating",
                "-t",
                str(small_tree_path),
                "-s",
                str(small_alignment_path),
                "-c",
                str(cal),
                "-o",
                str(tmp_test_dir / "out"),
                "--method",
                "mcmctree",
                "--dry-run",
            ]
        )
        assert result.returncode != 0, result.stderr
