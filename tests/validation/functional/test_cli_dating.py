"""Functional tests for the 'phylodater dating' CLI."""

from tests.validation.helpers import run_cli


class TestCliDating:
    def test_dating_dry_run(
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

    def test_dating_generates_config(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        config_path = tmp_test_dir / "generated.yaml"
        # --generate-config is checked after argument parsing, so provide all required args.
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
                str(tmp_test_dir / "out"),
                "--method",
                "mcmctree",
                "--generate-config",
                str(config_path),
            ]
        )
        assert result.returncode == 0, result.stderr
        assert config_path.exists()

    def test_dating_force_overwrite(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        marker = tmp_test_dir / ".phylodater_run"
        marker.write_text("completed", encoding="utf-8")
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
                "--force",
            ]
        )
        assert result.returncode == 0, result.stderr

    def test_dating_no_clobber(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        marker = tmp_test_dir / ".phylodater_run"
        marker.write_text("completed", encoding="utf-8")
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
                "--no-clobber",
            ]
        )
        # --no-clobber skips existing outputs and exits successfully.
        assert result.returncode == 0, result.stderr

    def test_dating_force_and_no_clobber_mutually_exclusive(
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
                "--force",
                "--no-clobber",
            ]
        )
        assert result.returncode == 2, result.stderr

    def test_dating_low_memory(
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
                "--low-memory",
            ]
        )
        assert result.returncode == 0, result.stderr

    def test_dating_strip_annotations(
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
                "--strip-annotations",
            ]
        )
        assert result.returncode == 0, result.stderr

    def test_dating_auto_calibrate(
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
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--dry-run",
                "--auto-calibrate",
                "Bact:50",
            ]
        )
        assert result.returncode == 0, result.stderr

    def test_dating_with_taxonomy(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        taxonomy_table_path,
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
                "--taxonomy-file",
                str(taxonomy_table_path),
            ]
        )
        # The small tree labels do not match the taxonomy table; this is a data error.
        assert result.returncode in (0, 3), result.stderr
