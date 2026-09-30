"""End-to-end integration tests covering ALL dating tools supported by PhyloDater.

The tests run each supported method on a realistic benchmark dataset
(7-taxon primate tree + 500 bp alignment, `benchmark/` fixtures) and assert:

1. the PhyloDater CLI exits successfully (returncode == 0);
2. a dated tree file is produced per method;
3. the comparison table exists and reports a *finite, positive* root age
   (scientific sanity check: a dating run that returns zero or negative root
   ages is treated as a failure).

These tests require the corresponding external software to be installed; they
are skipped automatically when a binary/Python backend is missing
(see `tests/validation/helpers.py::external_available`).
"""

import csv
from pathlib import Path

import pytest

from tests.validation.helpers import run_cli


def _find_dated_tree(output_dir: Path) -> Path:
    """Return the first dated-tree file produced by a method run."""
    candidates = sorted(output_dir.glob("*dated_tree*"))
    assert candidates, f"no dated tree produced in {output_dir}"
    return candidates[0]


def _assert_comparison_has_numeric_age(output_dir: Path, method_column: str) -> float:
    """Assert the comparison table carries a finite positive *mean* age for the
    Root node of the given method, and return that value."""
    tsv = output_dir / "comparison_table.tsv"
    assert tsv.exists(), f"comparison_table.tsv missing in {output_dir}"
    with tsv.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    root_rows = [r for r in rows if r.get("Node", "").lower() == "root"]
    assert root_rows, f"no Root row in comparison table for {method_column}"
    mean_col = f"{method_column}_mean"
    root = root_rows[0]
    mean = root.get(mean_col)
    assert mean not in (
        None,
        "",
        "NA",
    ), f"Root mean age for {method_column} is {mean!r} (expected a number)"
    value = float(mean)
    assert value > 0, f"Root age for {method_column} must be > 0, got {value}"
    return value


@pytest.mark.slow
@pytest.mark.requires_mcmctree
class TestMCMCTreeBenchmark:
    def test_mcmctree_end_to_end(
        self,
        benchmark_tree_path,
        benchmark_alignment_path,
        benchmark_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(benchmark_tree_path),
                "-s",
                str(benchmark_alignment_path),
                "-c",
                str(benchmark_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mcmctree",
                "--method-args",
                "burnin=200,nsample=500,sampfreq=2,skip_rate_estimation=true,paml_version=4.10.10",
                "--threads",
                "1",
            ]
        )
        assert result.returncode == 0, result.stderr
        _find_dated_tree(tmp_test_dir)
        _assert_comparison_has_numeric_age(tmp_test_dir, "mcmctree")


@pytest.mark.slow
@pytest.mark.requires_lsd2
class TestLSD2Benchmark:
    def test_lsd2_end_to_end(
        self,
        benchmark_tree_path,
        benchmark_alignment_path,
        benchmark_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(benchmark_tree_path),
                "-s",
                str(benchmark_alignment_path),
                "-c",
                str(benchmark_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "lsd2",
                "--method-args",
                "model=JC69",
                "--threads",
                "1",
                "--seed",
                "42",
            ]
        )
        assert result.returncode == 0, result.stderr
        _find_dated_tree(tmp_test_dir)
        _assert_comparison_has_numeric_age(tmp_test_dir, "lsd2")


@pytest.mark.slow
@pytest.mark.requires_r8s
class TestR8SBenchmark:
    def test_r8s_end_to_end(
        self,
        benchmark_tree_path,
        benchmark_alignment_path,
        benchmark_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(benchmark_tree_path),
                "-s",
                str(benchmark_alignment_path),
                "-c",
                str(benchmark_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "r8s",
                "--method-args",
                "method=NPRS",
                "--threads",
                "1",
                "--seed",
                "42",
            ]
        )
        assert result.returncode == 0, result.stderr
        _find_dated_tree(tmp_test_dir)
        _assert_comparison_has_numeric_age(tmp_test_dir, "r8s")


@pytest.mark.slow
@pytest.mark.requires_treepl
class TestTreePLBenchmark:
    def test_treepl_end_to_end(
        self,
        benchmark_tree_path,
        benchmark_alignment_path,
        benchmark_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(benchmark_tree_path),
                "-s",
                str(benchmark_alignment_path),
                "-c",
                str(benchmark_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "treepl",
                "--method-args",
                "cviter=5",
                "--threads",
                "1",
                "--seed",
                "42",
            ]
        )
        assert result.returncode == 0, result.stderr
        _find_dated_tree(tmp_test_dir)
        _assert_comparison_has_numeric_age(tmp_test_dir, "treepl")


@pytest.mark.slow
@pytest.mark.requires_pathd8
class TestPathD8Benchmark:
    def test_pathd8_end_to_end(
        self,
        benchmark_tree_path,
        benchmark_alignment_path,
        benchmark_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(benchmark_tree_path),
                "-s",
                str(benchmark_alignment_path),
                "-c",
                str(benchmark_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "pathd8",
                "--threads",
                "1",
            ]
        )
        assert result.returncode == 0, result.stderr
        _find_dated_tree(tmp_test_dir)
        _assert_comparison_has_numeric_age(tmp_test_dir, "pathd8")


@pytest.mark.slow
@pytest.mark.requires_wlogdate
class TestWLogDateBenchmark:
    def test_wlogdate_end_to_end(
        self,
        benchmark_tree_path,
        benchmark_alignment_path,
        benchmark_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(benchmark_tree_path),
                "-s",
                str(benchmark_alignment_path),
                "-c",
                str(benchmark_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "wlogdate",
                "--method-args",
                "max_iter=2000",
                "--threads",
                "1",
                "--seed",
                "42",
            ]
        )
        assert result.returncode == 0, result.stderr
        _find_dated_tree(tmp_test_dir)
        _assert_comparison_has_numeric_age(tmp_test_dir, "wlogdate")


@pytest.mark.slow
@pytest.mark.requires_mdcat
class TestMDCatBenchmark:
    def test_mdcat_end_to_end(
        self,
        benchmark_tree_path,
        benchmark_alignment_path,
        benchmark_calibration_path,
        tmp_test_dir,
    ):
        result = run_cli(
            [
                "dating",
                "-t",
                str(benchmark_tree_path),
                "-s",
                str(benchmark_alignment_path),
                "-c",
                str(benchmark_calibration_path),
                "-o",
                str(tmp_test_dir),
                "--method",
                "mdcat",
                "--method-args",
                "ncat=2,nrep=3,max_iter=5,ci_nboots=0,use_direct_import=false",
                "--threads",
                "1",
                "--seed",
                "42",
            ]
        )
        assert result.returncode == 0, result.stderr
        _find_dated_tree(tmp_test_dir)
        _assert_comparison_has_numeric_age(tmp_test_dir, "mdcat")
