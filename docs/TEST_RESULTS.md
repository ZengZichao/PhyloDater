> **Language:** English (canonical) · [中文](TEST_RESULTS_CN.md)

# Test Results

Results of the validation campaign described in [TESTING.md](TESTING.md). Every
number below was produced by the commands in §8 on the machine and
environment recorded in §1, and is reproducible from the shipped artefacts in
[`test-data/benchmark/results/`](../test-data/benchmark/results/).

**Release:** PhyloDater 0.1.0 | **Date:** 2026-09-29 | **Verdict: PASS** —
all seven supported dating methods run end-to-end against their real engines and
satisfy the acceptance gates.

---

## 1. Environment

| item | value |
|------|-------|
| OS / kernel | Ubuntu 22.04.4 LTS, Linux 6.8.0-111-generic, x86_64 |
| Hardware | 48 physical / 96 logical cores, 772 GB RAM |
| Python | 3.12.14 (conda env `phylodater`) |
| Package managers | micromamba 2.9.0 (inside the env), pip 24.2 |
| Core libraries | numpy 2.5.3, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2, DendroPy 5.1.0, ete3 3.1.3, PyYAML 6.0.3, biopython 1.88 |
| Test framework | pytest 9.1.1, pytest-cov 7.1.0 |

### 1.1 Dating engines actually executed

| method | engine and version recorded in `runtime_metadata.json` | install route |
|--------|--------------------------------------------------------|---------------|
| `mcmctree` | `MCMCTREE in paml version 4.10.10, 29 Jan 2026` | bioconda `paml` |
| `lsd2` | `IQ-TREE version 3.1.3 for Linux x86 64-bit built Jul 26 2026` | bioconda `iqtree` |
| `r8s` | `pyr8s 0.3.1 (NPRS only)` | pip from source (`iTaxoTools/pyr8s`) |
| `treepl` | `treePL version 1.0` (2.6.3 build) | historic conda build + `patchelf`, or source build |
| `pathd8` | `PATHd8 (no version banner): <prefix>/bin/PATHd8` | built from official source |
| `wlogdate` | `wLogDate 1.0.2` | PyPI `wlogdate` |
| `mdcat` | `1.0.1` | pip from source (`uym2/MD-Cat`) |

Provenance, checksums and the reason each route was chosen:
[EXTERNAL_TOOLS.md](EXTERNAL_TOOLS.md).

## 2. Reference benchmark — headline result

Dataset: `test-data/benchmark` (13 taxa, 3 000 sites, JC69 strict clock,
root = 105 Ma, 12 internal nodes, 4 calibrations; see TESTING.md §3).
Gate: MAE ≤ 10 Ma, worst-node relative error ≤ 0.5, every requested calibration
either honoured or explicitly disclosed as not honoured.

### 2.1 Variant A — input tree with the *expected* (noise-free) clock lengths

| method | status | wall (s) | MAE (Ma) | RMSE (Ma) | max abs err (Ma) | max rel err |
|--------|--------|-----------|----------|-----------|------------------|-------------|
| `mcmctree` | PASS | 14.2 | 0.721 | 0.863 | 1.641 | 0.284 |
| `lsd2` | PASS | 4.2 | 1.881 | 2.373 | 4.912 | 0.303 |
| `r8s` (pyr8s/NPRS) | PASS | 4.4 | 0.183 | 0.256 | 0.528 | 0.013 |
| `treepl` | PASS | 64.4 | 1.099 | 1.527 | 3.044 | 0.029 |
| `pathd8` | PASS | 4.4 | 0.509 | 0.794 | 1.659 | 0.019 |
| `wlogdate` | PASS | 6.1 | 0.000 | 0.000 | 0.000 | 0.000 |
| `mdcat` | PASS | 10.5 | 0.081 | 0.106 | 0.163 | 0.069 |

All seven: tips complete, topology identical, ages non-negative, chronology
monotone, calibrations honoured, all 12 truth nodes recovered.

`wlogdate` returns exactly the truth here, and that is *not* a sign of superior
accuracy: this variant feeds the method a branch-length vector that already is
proportional to time, so recovering the ages is an arithmetic rescaling. It is
kept for exactly that reason — it is the control that proves the scoring pipeline
cannot do better than the input it is given.

### 2.2 Variant B — input tree with ML-estimated branch lengths (the realistic case)

Same alignment, but the input tree's branch lengths were estimated by
maximum likelihood (GTR+G, topology fixed) with IQ-TREE 2 and re-rooted as
described in TESTING.md §3. Total tree length drops from 1.0365 to 0.8793
subs/site, so the input itself is now biased and noisy.

| method | status | wall (s) | MAE (Ma) | RMSE (Ma) | max abs err (Ma) | max rel err | disclosed limitations |
|--------|--------|-----------|----------|-----------|------------------|-------------|------------------------|
| `mcmctree` | PASS | 14.0 | 0.730 | 0.878 | 1.650 | 0.289 | — |
| `lsd2` | PASS | 4.1 | 1.881 | 2.373 | 4.912 | 0.303 | — |
| `r8s` (pyr8s/NPRS) | PASS | 4.8 | 7.862 | 11.219 | 25.023 | 0.426 | `Primates` 54.64 Ma outside [68, 82]; `Glires` 53.75 Ma outside [60, 80] |
| `treepl` | PASS | 63.0 | 4.041 | 5.762 | 14.866 | 0.224 | — |
| `pathd8` | PASS | 4.1 | 4.492 | 6.262 | 15.652 | 0.240 | — |
| `wlogdate` | PASS | 6.0 | 1.850 | 3.199 | 9.666 | 0.180 | — |
| `mdcat` | PASS | 18.5 | 0.298 | 0.472 | 1.116 | 0.413 | — |

Reading the table honestly:

* **Variant B is the result to quote.** Every method degrades except `mcmctree`
  and `mdcat`, which re-estimate rates instead of trusting the supplied lengths.
* MCMCTree's figures move in the third decimal from run to run (0.721 in the
  archived run; 0.731 and 0.756 in two earlier captures) because it samples; every
  other row is bit-reproducible (§3). The archived `scorecard.tsv` is the
  authoritative record of what was last run.
* **`lsd2` is unchanged between variants** because IQ-TREE 2 re-optimises branch
  lengths under `GTR+G` on the fixed topology; the input lengths are not used as
  data. Its MAE of 1.88 Ma is the estimator's own error.
* **`r8s`'s range violations are an engine property, not a silent failure.**
  pyr8s implements `CONSTRAIN min_age/max_age` as a *relaxed barrier penalty*,
  so NPRS may end outside the requested range, while `FIXAGE` is enforced exactly
  (measured: `Homininae` fixage 9 Ma → 9.000 Ma in both variants). The experiment
  that establishes this is shipped as
  [`scripts/probe_pyr8s_constraints.py`](../scripts/probe_pyr8s_constraints.py).
  PhyloDater now (a) warns before the run that pyr8s ranges are advisory and
  (b) verifies each requested range against the returned ages and records a
  warning naming every unfulfilled calibration. The benchmark's gate is
  "honoured **or disclosed**", and the disclosed column above is what the
  adapter reported.
* **`pathd8` root behaviour improved during this campaign.** The adapter used to
  drop any non-point root calibration ("PATHd8 has no min/max directive for the
  root"), which left the root at 153 Ma against a requested [100, 110]. PATHd8
  does accept `minage`/`maxage` on the root MRCA — verified directly against the
  binary — so the range is now written and the root is returned at 110 Ma, inside
  the requested interval (variant B; MAE 4.49 Ma over all 12 nodes).

### 2.3 Post-type-refactor timing re-check (no performance regression)

The type-debt work (P0) touched every adapter's constructor and 286 function
definitions (the changed `def` lines of annotation pass `cb4fd01`), so the wall
times were re-measured afterwards — same host, same conda environment, three
independent variant-A runs plus two variant-B runs:

| method | variant A wall (s), three runs | variant B wall (s), two runs | archived baseline (A / B) |
|--------|--------------------------------|------------------------------|---------------------------|
| `mcmctree` | 15.07 / 13.86 / 14.41 | 14.60 / 14.31 | 15.07 / 14.60 |
| `lsd2` | 4.39 / 4.40 / 3.86 | 4.35 / 4.13 | 4.39 / 4.35 |
| `r8s` | 4.31 / 4.68 / 4.54 | 4.56 / 4.44 | 4.31 / 4.56 |
| `treepl` | **64.50 / 67.81 / 59.85** | 55.78 / 63.50 | 64.50 / 55.78 |
| `pathd8` | 4.35 / 4.36 / 4.20 | 4.35 / 4.25 | 4.35 / 4.35 |
| `wlogdate` | 5.81 / 5.28 / 5.48 | 6.03 / 6.01 | 5.81 / 6.03 |
| `mdcat` | 10.51 / 10.29 / 10.05 | 18.37 / 17.74 | 10.51 / 18.37 |

`treepl` remains the dominant cost (55–68 s per run; the cross-validation search
over smoothing values, not the wrapper, is what takes the time), and its spread is
the same before and after the refactor — no measurable regression. Accuracy is
unchanged in every deterministic method: the variant-A MAE over
`lsd2` / `r8s` / `treepl` / `pathd8` / `wlogdate` / `mdcat` reproduces bit-for-bit
at 1.881 / 0.183 / 1.099 / 0.509 / 0.000 / 0.081 Ma.

## 3. Determinism

Two independent full benchmark runs with `--seed 42`, compared node by node
(`node_ages.tsv`):

| method | repeated runs |
|--------|---------------|
| `lsd2`, `r8s`, `treepl`, `pathd8`, `wlogdate`, `mdcat` | bit-identical node ages |
| `mcmctree` | differs by ≤ 0.3 Ma (≤ 0.3 % of root age) — a Bayesian sampler is not bit-reproducible from a PRNG seed alone |

Re-verified after the type refactor by hashing each method's `node_ages.tsv` rows
across two variant-A runs: `lsd2`, `r8s`, `treepl`, `pathd8`, `wlogdate` and
`mdcat` were **byte-for-byte identical**; only `mcmctree` moved
(MAE 0.705 vs 0.741 Ma).

MCMCTree quality is therefore judged on its convergence diagnostics, which
PhyloDater records: `max PSRF = 0.99988`, `min ESS = 423` over 2 chains ×
4 001 samples (thresholds 1.05 / 100), verdict `converged`.

## 4. Automated test suites

Three counts are reported, deliberately, because they answer different questions and
they are not interchangeable:

| accounting | collected | passed | skipped | failed |
|---|---|---|---|---|
| **release snapshot** — §1's validated environment (Ubuntu 22.04, Python 3.12.14, conda env `phylodater`, all seven engines present) at the `v0.1.0` tag; the figure quoted by the manuscript and the SI | 1,489 | 1,487 | 2 | 0 |
| **v0.1.0 clean tree, bare environment** — `git archive v0.1.0` unpacked into an empty directory (no uncommitted working-tree changes, no fix branch), measured on macOS arm64, Python 3.11.15, pytest 9.1.1 with **no external engine installed and no `python` alias on PATH** (only `python3`); 94 s; the quotable coverage figure in §4.2 comes from this run | 1,495 (1) | 1,430 | 36 | **29** |
| **current tree** — the review-fix branch working tree, measured on macOS / Python 3.11.15 with pytest 9.1.1 | 1,550 | 1,514 | 36 | 0 |

(1) The 1,495 is the total of the three outcome counts of that run
(29 failed + 1,430 passed + 36 skipped), not an independently collected figure.

```bash
python -B -m pytest tests -q              # run everything
python -B -m pytest tests -q --collect-only   # counts only, no execution
```

> **The collected count is environment-sensitive; treat it as a measured value, not
> a constant.** Parametrised and marker-gated cases come and go with the
> interpreter minor version and which optional backends are importable — `ete3`
> (not installable for Python ≥ 3.13), `ete4`, `plotly` and `psutil` each change
> what can be collected. An intermediate measurement of this same tree, before the
> review fixes added unit tests, gave 1,495 items (unit 1,114 / functional 199 /
> integration 24 / validation 158); the current tree gives 1,550 (unit **1,169** /
> functional 199 / integration 24 / validation 158). Only `tests/unit` moved. When
> quoting a count, quote the interpreter and the command with it.

> **The 29 failures in the middle row are a test-harness defect, not 29 broken
> features** — see the reproducibility disclosure in §4.2. They disappear on the
> review-fix branch (`94b9781`), which has **not** been released, so the public
> `v0.1.0` archive still behaves this way.

### 4.1 Breakdown

* `tests/unit` + `tests/functional` + `tests/integration` — the developer suite:
  models, constraint rendering, adapters' input writing/parsing, CLI, reports,
  security, visualisation, checkpointing.
* `tests/validation` — the release-validation suite: per-adapter contract tests,
  CLI subcommand tests (`dating`, `check`, `fold`), error-path tests, security
  tests (path traversal, YAML safety, input safety), and end-to-end runs of every
  installed engine on the benchmark fixtures.
* Engine-dependent tests skip with an explicit reason when a binary is absent, so
  the same command works on a machine without the tools installed.

### 4.2 Coverage

```bash
python -B -m pytest tests -q --cov=phylodater --cov-report=term
```

Four coverage readings exist for this repository. **Only reading 1 is quoted by the
manuscript and the SI**; readings 2–4 are tabulated so that anyone who lands on a
different number can see exactly which source tree and which data file produced it.

| # | reading | statements | missing | coverage | measured on | status |
|---|---|---|---|---|---|---|
| **1** | **v0.1.0 clean tree** | 16,729 | 5,483 | **67%** | `git archive v0.1.0` unpacked into an empty directory — no uncommitted working-tree changes, no fix branch — on macOS arm64, Python 3.11.15, pytest 9.1.1; **no external engine installed, no `python` alias on PATH (only `python3`)**; the run is the middle row of §4 (29 failed / 1,430 passed / 36 skipped, 94 s) | **the value the manuscript and the SI quote** — it confirms the "approximately 67%" originally written for v0.1.0 |
| **2** | review-fix branch | 16,892 | 5,325 | 68.48% | working tree of `fix/review-2026-10-08` on the same machine and interpreter; 1,550 collected / 1,514 passed / 36 skipped / 0 failed (re-run on the final branch tree, same figures) | true, but **a branch value**: it describes no released artefact and must be labelled as such |
| **3** | stale `.coverage` read against the pre-fix tree | 16,742 | 5,631 | 66.4% | the git-ignored `.coverage` data file shipped in the working tree, reported *after* the fact | **invalid** — data file and source tree do not match; recorded only to explain where the number logged during the review came from |
| **4** | stale `.coverage` read against the branch tree | 16,892 | 7,054 | 58% | the very same stale data file, reported against the branch tree | **invalid** — same mismatch; the lines the branch added have no recorded hits at all |

| item | value |
|---|---|
| quotable measurement | reading 1: **67%** on the public `v0.1.0` clean tree |
| gate | `fail_under = 60` in `pyproject.toml` — satisfied by both quotable readings (1 and 2). The invalid reading 4 would fall below it, which is another reason it is not a coverage result |
| command | `python -B -m pytest tests -q --cov=phylodater --cov-report=term` (`python3 -B -m pytest …` where no `python` alias exists) |

Four caveats, none of which may be dropped when this number is quoted:

1. **It is a lower bound, not the true figure.** `tests/functional` and
   `tests/validation` drive the CLI through `subprocess`, and no coverage
   subprocess hook is configured anywhere in the repo (`COVERAGE_PROCESS_START` is
   never set), so child-process execution is scored as unexecuted. Concretely,
   `phylodater/infrastructure/self_test.py` reports **0%** and
   `phylodater/cli/main.py` **36%**, even though `phylodater check` passes 26/26 —
   that code only ever runs in a subprocess. Enabling subprocess coverage and
   `coverage combine` would raise the total.
2. **Readings 3 and 4 are not coverage results and may not be quoted.** They come
   from the checked-out `.coverage` file, which is git-ignored (the `.coverage*`
   entry in `.gitignore`), is a *partial* recording that predates the fixes in this
   branch, and is being reported against source trees it was never produced by.
   The statement base also moves (16,742 vs 16,892) purely because different code is
   in the tree, so the same stale file yields 66.4% on one tree and 58% on another.
   The cluster's earlier "correct to 66%" edit to the manuscript was made on reading 3
   and was therefore an **incorrect correction**: reading 1 supersedes it, and the
   manuscript's "approximately 67%" is right for v0.1.0. Always regenerate
   (`--cov` above) instead of re-reporting a stale data file; if you want to keep the
   historical artifact in place, point the data file elsewhere first
   (`export COVERAGE_FILE=$PWD/../.coverage.regen`).
3. **The number moves with the environment, like the collected count in §4.** The
   comment beside `fail_under` in `pyproject.toml` quotes 67.87%; that was measured
   in §1's all-seven-engine Linux environment at the `v0.1.0` tag, and it is not
   reproducible on the current tree without that environment. Quote a coverage number
   together with its tree *and* its interpreter, as reading 1 does.
4. **A disclosed reproducibility defect sits under reading 1.** The published
   `v0.1.0` — Zenodo archive and GitHub tag — cannot reproduce a clean test run on any
   machine that has only `python3` and no `python` alias: **28 of the 29 failures**
   come from the validation harness itself, `tests/validation/helpers.py::run_cli()`,
   which hard-codes the interpreter as the literal string `'python'` (line 75 of
   `v0.1.0`, the `FileNotFoundError: No such file or directory: 'python'` being raised
   at line 76 where `subprocess.run()` is called). They are distributed as 18 in
   `tests/validation/functional`, 9 in `tests/validation/security`, 1 in
   `tests/validation/integration`, and they fail *before* any assertion, so they
   measure nothing about the package. The remaining failure is an `ete4`
   incompatibility: `tests/unit/test_viz_geo.py::TestRootAgeInference::test_documented_example_with_a_real_ete_backend`
   builds `Tree(newick, format=1)`, which ete4's `Tree()` does not accept.
   Branch commit `94b9781` replaces the literal with `sys.executable` (matching what
   `tests/functional/conftest.py` already did) and adds the missing
   `requires_mcmctree` marker, which takes the branch to 0 failed — but that commit
   **has not been released**, so the archived public version remains affected. Until a
   release carries it, a reviewer running §8 on such a machine should expect those 29
   failures and should invoke `python3 -B -m pytest tests -q --cov=phylodater`.

## 5. Self-test of the installed package

```
$ python -m phylodater check
...
总计: 26 项检查, 26 通过, 0 失败
✓ 所有检查通过！
```

Those last two lines are the program's verbatim Chinese output (the summary is printed
in Chinese by design and is quoted here unmodified); they read "Total: 26 checks, 26
passed, 0 failed" and "✓ All checks passed!".

The 26 checks cover the required Python dependencies (biopython, dendropy, numpy,
pandas, pyyaml, matplotlib), the four optional engine dependencies — `logdate`
(wLogDate), `emd` (MD-Cat), `pyr8s` (r8s) and `ete3`, probed but never required for
the Python layer to work — Newick parsing with and without embedded taxonomy,
both taxonomy layouts (embedded GTDB and the tabular semicolon form, including
correct missing-rank handling: `species=None` is the contract, not a bug),
monophyly detection (non-monophyletic, monophyletic and the LUCA special
identifier), deep validation (bracket balance, duplicate ids), the LUCA/ROOT
special identifiers, the 500 MB input-size guard, malicious-character rejection,
multi-tree NEXUS handling, external taxonomy tables, and the library-mode
cross-validation API.

## 6. Defects found by this campaign and fixed

| # | defect | symptom | fix | pinned by |
|---|--------|---------|-----|-----------|
| 1 | LSD2 default substitution model was the amino-acid model `LG+G` | IQ-TREE 2: `ERROR: File not found LG`; the method could never run on DNA | model is resolved from the alignment type (`GTR+G` for nucleotide, `LG+G` for protein); an obviously mismatched explicit model raises `ConfigurationError` | `TestLSD2ModelResolution` (5 tests) |
| 2 | MD-Cat `--CI` passed as one comma-joined token | upstream `int("100,0.025,0.975")` → `ValueError`, every default run failed | pass `--CI "n lower upper"` (single whitespace-separated argument, as upstream's `.split()` expects) | `TestMDCatCliContract` |
| 3 | MD-Cat CI defaulted to 100 bootstrap replicates | hours of unannounced work per run (measured: > 10 min for a *single* replicate on 6 taxa) | default `ci_nboots = 0`, matching upstream; opt-in documented | `test_ci_parameters` |
| 4 | PAML gate demanded the literal version `4.10.8` | the documented `conda install -c bioconda paml` (4.10.10) was refused: "Environment not available" | compatibility decided by major.minor series; same-series patch differences pass with an informational log; cross-series still blocked | 2 new tests in `tests/unit/adapters/test_mcmctree_method.py` |
| 5 | `runtime_metadata.json` never recorded engine versions | shipped results had no engine provenance | each adapter records what `validate_environment()` detected; the pipeline writes it | benchmark scorecard + adapter tests |
| 6 | The zero-dependency Newick walker rejected a leading `[&R]` | MCMCTree and treePL aborted with "Cannot parse Newick for tree operations" on trees exported by BEAST/IQ-TREE/DendroPy; ancestry validation silently degraded | leading comments are skipped by the walker, and a shared `strip_leading_newick_comments()` normalises them for ete3 | `TestLeadingCommentTolerance` (5 tests) |
| 7 | PATHd8 dropped non-point root calibrations | root age ignored: 153 Ma returned against a requested [100, 110] Ma | root ranges are now rendered as `minage`/`maxage` on the root MRCA (behaviour verified against the binary) | `tests/unit/adapters/test_pathd8_method.py` |
| 8 | pyr8s range calibrations were neither enforced nor reported | ages outside the requested range reached the report without a word | pre-run disclosure + post-run per-calibration verification in the result warnings | `tests/unit/adapters/test_r8s_pyr8s_method.py` |
| 9 | three unit tests asserted on the process-wide `sys.modules` | they failed purely because of test order (a visualisation test imports ete3) | `module_imports_cleanly` fixture probes a fresh interpreter | the three converted tests |
| 10 | raw MCMCTree smoke test wrote `BDparas = 1 1 0.1` | PAML ≥ 4.10 aborts: "BDparas: expect flag for birth-death process prior" | control file now includes the construction flag, mirroring the adapter | `tests/integration/full_pipeline/test_all_dating_tools.py` |
| 11 | the validation harness hard-coded its own interpreter as the literal `'python'` (`tests/validation/helpers.py::run_cli()`), and `TestMultiMethod` carried only `slow` instead of `requires_mcmctree` | on a machine with only `python3` — no `python` alias — 28 `tests/validation` cases raised `FileNotFoundError: No such file or directory: 'python'` before reaching any assertion; without PAML the MCMCTree integration case reported a failure instead of a skip | `run_cli()` now uses `sys.executable` (as `tests/functional/conftest.py` already did) and the missing availability marker is added. **Not released**: the fix lives in branch commit `94b9781`, so the public `v0.1.0` archive is still affected (§4.2 caveat 4) | branch commit `94b9781` (unreleased) |

## 7. Engine- and environment-independent limitations

* treePL's self-reported banner is `treePL version 1.0` even for the 2.6.3 build;
  version provenance therefore also relies on the recorded executable path.
* PATHd8 prints no version at all.
* The Docker image ships 3 of the 7 engines (EXTERNAL_TOOLS.md §7).
* Accuracy values above are properties of *this* dataset (13 taxa, 3 kb, strict
  clock). They must not be quoted as general accuracy claims for the methods.
* A method returning MAE 0 on a noise-free input tree (variant A, `wlogdate`) is a
  property of the input, not of the estimator — variant B is what a user should
  expect on an ML tree.

## 8. Reproducing everything in this document

```bash
conda activate phylodater            # env with all seven engines (§1.1)

# variant A
python scripts/run_reference_benchmark.py \
    --output test-data/benchmark/results

# variant B (ML branch lengths)
python scripts/run_reference_benchmark.py --tree-file tree_ml.nwk \
    --output test-data/benchmark/results/ml_tree

# engine-semantics experiment behind §2.2
python scripts/probe_pyr8s_constraints.py

# automated suites + self test
python -m pytest tests -q
phylodater check
```

> If you are reproducing **§4.2 reading 1** (the 67% the manuscript quotes), do it on
> a clean unpacking of the tag (`git archive v0.1.0 | tar -x`) in an environment with
> no external engine and no `python` alias, and use `python3 -B -m pytest tests -q
> --cov=phylodater --cov-report=term`. Expect `29 failed, 1,430 passed, 36 skipped`
> and `TOTAL 16,729 statements / 5,483 missing = 67%`; the 29 failures are the harness
> defect disclosed in §4.2 caveat 4, not regressions.

Artefacts committed with the package:
`test-data/benchmark/results/{scorecard.tsv,node_ages.tsv,benchmark_results.json}`
and the same three files under `results/ml_tree/`. The `runs/` sub-directories
(each method's dated tree, comparison table, `runtime_metadata.json`, log) are
regenerated on demand and are git-ignored.

---

*See also: [TESTING.md](TESTING.md), [EXTERNAL_TOOLS.md](EXTERNAL_TOOLS.md),
[Chinese version](TEST_RESULTS_CN.md).*
