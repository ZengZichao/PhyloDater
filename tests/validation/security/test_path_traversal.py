"""Path traversal and output isolation security tests."""

import pytest

from tests.validation.helpers import run_cli


@pytest.mark.security
class TestPathTraversal:
    def test_output_confined_to_specified_directory(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        """A dry run may create logging artifacts, but all must stay inside the output directory."""
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
        # Nothing should escape tmp_test_dir
        for path in tmp_test_dir.rglob("*"):
            assert str(path.resolve()).startswith(str(tmp_test_dir.resolve()))

    @pytest.mark.requires_mcmctree
    def test_output_directory_isolation(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        """A real run must create its artifacts only under the specified output directory.

        这条与下面 ``--force`` 那条都真的跑一次 MCMCTree（断言 returncode == 0
        且产物存在），所以必须标 ``requires_mcmctree``：本套件的既定约定是
        "依赖外部引擎的测试在没有该引擎时带原因跳过，而不是失败"（conftest 的
        ``pytest_runtest_setup`` 就是干这个的）。不标的话，在未装 PAML 的机器上
        （例如从 sdist 干装后跑 ``pytest tests/validation -m "not slow"``）它们会
        报失败。装了 PAML 的环境里两者照旧执行，断言未变。
        """
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
                "burnin=100,nsample=50,paml_version=4.10.10",
                "--threads",
                "1",
            ]
        )
        assert result.returncode == 0, result.stderr
        for path in tmp_test_dir.rglob("*"):
            assert str(path.resolve()).startswith(str(tmp_test_dir.resolve()))

    def test_no_clobber_respects_existing_output(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        """Without --force, CLI must refuse to overwrite an existing PhyloDater output directory."""
        marker = tmp_test_dir / ".phylodater_run"
        marker.write_text("previous run", encoding="utf-8")
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
                "burnin=100,nsample=50,paml_version=4.10.10",
                "--threads",
                "1",
            ]
        )
        assert result.returncode != 0, result.stderr
        assert marker.exists()
        assert marker.read_text(encoding="utf-8") == "previous run"

    @pytest.mark.requires_mcmctree
    def test_force_overwrite_existing_output(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        """With --force, CLI may overwrite an existing output directory and produce results.

        同上：这条会真跑一次 MCMCTree 并断言产物存在，因此需要引擎可用；
        没有引擎时按套件约定跳过而不是失败。
        """
        marker = tmp_test_dir / ".phylodater_run"
        marker.write_text("previous run", encoding="utf-8")
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
                "burnin=100,nsample=50,paml_version=4.10.10",
                "--threads",
                "1",
                "--force",
            ]
        )
        assert result.returncode == 0, result.stderr
        assert (
            (tmp_test_dir / "mcmctree_dated_tree.nexus").exists()
            or (tmp_test_dir / "mcmctree_dated_tree.nwk").exists()
            or (tmp_test_dir / "mcmctree_dated_tree.nhx").exists()
        )
