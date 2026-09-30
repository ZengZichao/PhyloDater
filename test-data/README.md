> **Language:** English (canonical) · [中文](README_CN.md)

# test-data

Test data shipped **with** the package, so that the validation recorded in
[docs/TEST_RESULTS.md](../docs/TEST_RESULTS.md) can be re-run (or re-checked)
without reconstructing anything.

Two kinds of data live here, and they serve different purposes:

| path | role | consumed by |
|------|------|-------------|
| `benchmark/` | **quantitative**: a simulated dataset with known divergence times, used to score how close each method's estimates are to the truth | `scripts/run_reference_benchmark.py`, `tests/validation/integration/` |
| `input/` | **qualitative**: valid / invalid / boundary scenario inputs for manual end-to-end checks of the workflow and its error handling | QA walkthroughs, `tests/functional/` fixtures |

There are no automated tests in this directory; the tests are in [`tests/`](../tests/).

## `benchmark/` — the ground-truth accuracy dataset

| file | content |
|------|---------|
| `tree.nwk` | rooted binary tree, 13 taxa, branch lengths = expected substitutions/site under a strict clock |
| `tree_ml.nwk` | the same topology with **maximum-likelihood** branch lengths estimated from the alignment — the realistic input a user would supply |
| `alignment.fasta` | 13 × 3 000 nt alignment simulated under JC69 on the tree above |
| `calibrations.yaml` | four calibrations on four distinct nodes (one fixed, three uniform ranges) |
| `truth.yaml` | generation parameters (seed, rate, model, length) plus the true age of all 12 internal nodes, keyed by MRCA leaf pair |
| `generate_benchmark_dataset.py` | deterministic generator; `--include-ml-tree` regenerates `tree_ml.nwk` with IQ-TREE 2 |
| `results/` | archived outcome of the release run: `scorecard.tsv`, `node_ages.tsv`, `benchmark_results.json`, plus `ml_tree/` for the ML variant |
| `results/PYTHON_SUPPORT_MATRIX.txt` | evidence behind the declared Python support matrix: per-version dependency resolution and the `tests/unit` result measured on 3.8–3.14 |
| `results/validation_log.txt` | verbatim output of the final verification command set (pytest, `phylodater check`, formatters, `check_docs.py`, mypy, both benchmark variants) |
| `checksums.sha256` | SHA-256 of every file above (per-run scratch under `results/runs/` is excluded) |

Why a simulated dataset: real fossil calibrations cannot serve as ground truth,
because each method would be scored against the very numbers it was given. A
strict-clock JC69 simulation gives each estimator an unambiguous reference. The
topology and ages are synthetic but plausible; they are a numerical reference, not
a claim about evolutionary history. Full design rationale and the acceptance gates:
[docs/TESTING.md](../docs/TESTING.md) §3–§4.

Regenerate and re-score:

```bash
python test-data/benchmark/generate_benchmark_dataset.py --force --include-ml-tree
sha256sum -c test-data/benchmark/checksums.sha256
python scripts/run_reference_benchmark.py --output test-data/benchmark/results
python scripts/run_reference_benchmark.py --tree-file tree_ml.nwk \
    --output test-data/benchmark/results/ml_tree
```

`results/runs/` (per-method dated trees, comparison tables, logs) is regenerated on
demand and git-ignored; the summary files plus the two evidence files above are what
ships with the release.

`results/validation_log.txt` is not produced by the benchmark harness — it is the
archived terminal output of the release-gate commands listed in
[docs/TESTING.md](../docs/TESTING.md) §5, captured in one pass on the validation
host and committed unchanged (its header records the `git HEAD` it was run against;
the log file and its checksum entry are necessarily added by the following commit).
Re-capture it the same way before any release.

## `input/` — scenario inputs

| path | purpose |
|------|---------|
| `input/valid/` | standard workflow input (tree, alignment, calibration config) |
| `input/invalid/` | inputs that must trigger error handling (illegal characters, malformed formats, duplicate labels, …) |
| `input/boundary/` | boundary conditions (empty tree, multiple trees, BOM, GBK encoding, …) |

Usage: pick files from the matching subdirectory, run `phylodater dating`, and check
the reported behaviour (exit code, log, disclosed warnings) against expectations.

## Relationship to automated testing

- `tests/` — pytest-driven automated tests (unit / functional / integration / validation), run in CI
- `test-data/` (here) — the data those tests and the benchmark harness read

The two are complementary: the data is reviewed for scientific content, the tests
for behaviour.
