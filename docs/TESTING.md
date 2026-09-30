> **Language:** English (canonical) · [中文](TESTING_CN.md)

# Testing Guide

This document specifies **how PhyloDater is tested**, why each test is designed
the way it is, and how to reproduce every result quoted in
[TEST_RESULTS.md](TEST_RESULTS.md). It is part of the released package: the test
data, the harness, the expectations and the measured results ship together so that
a reviewer can re-run the validation without reconstructing it.

Scope: functional validation of all seven supported dating methods, plus the
unit / integration / security layers that protect them.

---

## 1. What "supported" means, and what must be tested

PhyloDater claims seven `dating` methods. Each one is a thin, standardised facade
over a different engine, so the test matrix has to cover three independent things:

1. **the facade** — does PhyloDater build valid engine input, run the engine, and
   parse its output into the common result model?
2. **the semantics** — does a calibration keep its meaning when it is translated
   into the engine's own syntax (and is every loss of meaning reported)?
3. **the estimate** — is the resulting chronogram *numerically* close to the truth?

Points 1 and 2 can be checked with mocks; point 3 cannot. It needs a dataset with
known divergence times, which is why the package ships a simulated ground-truth
benchmark (§3). Every one of the seven methods is exercised against the real
engine; none of the seven is tested only through a stub.

## 2. Test layers

```
tests/
├── unit/          29 files - pure Python: models, constraints, config, services, adapters' logic
├── functional/    12 files - CLI behaviour, validation, reports, security, visualisation
├── integration/    4 files - multi-step pipelines, plotting, raw-engine smoke tests
└── validation/    29 files - release-validation suite (real engines, end-to-end)
    ├── unit/ adapters/ functional/ integration/ security/
    ├── fixtures/        static inputs (trees, alignments, calibrations, configs)
    ├── helpers.py       engine-availability probes + CLI runner
    └── conftest.py      markers: slow, requires_mcmctree, requires_treepl, ...
```

| layer | external software | default `pytest tests` run |
|-------|-------------------|----------------------------|
| `unit`, `functional` | never invoked | always executed |
| `integration` | only for the "raw engine smoke" tests | executed when the binary exists, otherwise skipped |
| `validation` | real end-to-end runs | executed when the engine is installed, otherwise skipped with an explicit reason |

Tests that need an engine declare it, so a machine without PAML reports *skipped*,
never *passed*:

```python
@pytest.mark.slow
@pytest.mark.requires_mcmctree
class TestMCMCTreeBenchmark: ...
```

`helpers.external_available()` maps each method name onto what really has to be
present: `lsd2 → iqtree2/iqtree`, `r8s → pyr8s package or r8s binary`,
`wlogdate → logdate package or launch_wLogDate.py`, `mdcat → emd package or
md_cat.py`, `pathd8 → PATHd8`, `treepl → treePL`, `mcmctree → mcmctree + baseml`.

## 3. Reference benchmark dataset

`test-data/benchmark/` is a **simulated** dataset with known answers:

| file | content |
|------|---------|
| `tree.nwk` | rooted binary tree, 13 taxa, branch lengths = substitutions/site |
| `alignment.fasta` | 13 × 3000 nt multiple alignment |
| `calibrations.yaml` | four calibrations on four distinct nodes |
| `truth.yaml` | generation parameters + true age of all 12 internal nodes |
| `generate_benchmark_dataset.py` | the generator (deterministic, numpy-only) |
| `results/` | `scorecard.tsv`, `node_ages.tsv`, `benchmark_results.json` |

Design decisions, and why they are the scientifically defensible ones:

* **Strict clock, JC69, equal base frequencies.** Each method is then evaluated
  where its own assumptions hold. Deviations are attributable to the estimator,
  not to model violation. Rate autocorrelation, among-lineage rate variation or a
  non-representable substitution model would make "how far from the truth"
  uninterpretable.
* **Fossil-calibrated real data cannot be used as ground truth.** The methods would
  be scored against the same calibrations they were given. The simulated dataset
  breaks that circularity; the real-data examples stay in `examples/`.
* **Topology and ages are synthetic but plausible** (a eutherian scheme: 13 taxa,
  root 105 Ma). They are a numerical reference, not a claim about evolution.
* **Rate `1.5e-3` substitutions/site/Myr over 3000 sites** gives ≈ 470 expected
  differences root-to-tip: enough signal for ML branch lengths, far from saturation
  at the deepest node.
* **Seed `20260929`, stdlib + numpy only.** Re-running the generator reproduces
  byte-identical files, so the dataset's hash stays valid:

  ```bash
  python test-data/benchmark/generate_benchmark_dataset.py --force
  sha256sum test-data/benchmark/{tree.nwk,alignment.fasta,truth.yaml,calibrations.yaml}
  ```

* **Only point and two-sided-range calibrations.** `Homininae` fixed 9 Ma,
  `Primates` uniform 68–82, `Glires` uniform 60–80, `Root` uniform 100–110 (all
  centred on the truth). One-sided bounds (`maximum`, `soft_lower`) are *not* used
  here because engines that accept only point calibrations (wLogDate, MD-Cat)
  degrade them to the midpoint of the implied range, which can produce a chronology
  inconsistent with the remaining calibrations and abort inside the engine. Those
  constraint types are exercised deliberately in the constraint-coverage dataset
  under `examples/benchmark/degradation/`, where the degradation *warning* is the
  object under test.
* **No two calibrations resolve to the same node**, and no internal calibration
  lands on the root: MCMCTree writes the root age into its control file
  (`RootAge`) rather than into the tree, and refuses to run when a requested
  calibration cannot be found in the tree it was given.

## 4. Scoring methodology

`scripts/run_reference_benchmark.py` runs the CLI for each method and then checks,
in this order:

| # | check | how it is decided |
|---|-------|-------------------|
| 1 | completion | CLI exit code 0 **and** `runtime_metadata.json` reports `status: completed` |
| 2 | tips complete | the dated tree's taxon set equals the input's, no renames or drops |
| 3 | topology identical | set of internal clades (leaf sets) of output = input; PolyPhy-style clade matching, so node *labels* are irrelevant |
| 4 | ages non-negative | every internal node age ≥ 0 |
| 5 | chronology monotone | no child node older than its parent |
| 6 | calibrations honoured | fixed: \|estimate − age\| ≤ 1 Ma; range: within `[min, max]` widened by 5 % of the range width |
| 7 | accuracy | MAE, RMSE, max abs error and max relative error against `truth.yaml` over all 12 internal nodes |

Gate values used by the shipped benchmark: **MAE ≤ 10 Ma** and **worst-node
relative error ≤ 0.5**. They are deliberately loose enough to absorb statistical
noise on a 3 kb alignment and tight enough to fail a method that silently
mis-scales its output (the unit bug in check 8 below would, for example).

How node ages are recovered from each engine's output — this matters because a
wrong reading here would invalidate the whole benchmark:

* **Chronogram decoding.** Dated trees store time as branch *differences*, so
  `age(node) = max_root_to_tip_distance − distance_from_root(node)`. Using the raw
  distance from the root returns "time elapsed since the root" — 0 for the root and
  largest for the youngest nodes, i.e. exactly inverted ages.
* **Unit inference.** MCMCTree writes Ga, every other engine writes Ma. The scale
  factor is inferred from the root calibration (nearest of 1, 10⁻³, 10³, 10⁶)
  instead of being hard-coded, so a future engine with a different unit fails
  visibly rather than silently.
* **Clade matching.** Nodes are identified by their descendant leaf set, not by
  node names: several engines strip or rename internal labels
  (`IntNode7`, `I9`, none at all).
* Annotation stripping. `[&date=…]`, `[t=…, mu=…]`, `[&95%HPD={…}]` and internal
  labels are removed before re-parsing, and both DendroPy 4.x and 5.x label
  accessors are supported.

## 5. Reproducing the tests

```bash
# 1. Environment (all seven engines) - see EXTERNAL_TOOLS.md
conda activate phylodater

# 2. Self-test of the Python layer (26 checks, no external engine needed)
phylodater check

# 3. Full automated suite (1 489 collected: 1 487 passed, 2 skipped; engine-dependent
#    tests skip when missing)
python -m pytest tests -q

# 4. Validation suite only (real engines, end-to-end)
python -m pytest tests/validation -q

# 5. Reference benchmark with ground-truth scoring (regenerates results/)
python scripts/run_reference_benchmark.py \
    --dataset test-data/benchmark --output test-data/benchmark/results
```

Useful switches: `--methods` to select a subset, `--method-args METHOD=key=value,...`
to override one method's benchmark parameters (this is how the MD-Cat
sensitivity check below was run), `--timeout` per method, `--seed` (default 42)
forwarded to PhyloDater, `-m "not slow"` to deselect the engine-dependent tests.

Runtime on the reference machine (48 physical cores): ≈ 2 min for the whole
seven-method benchmark, of which treePL takes ≈ 1 min (cross-validation search)
and MCMCTree ≈ 15 s with the reduced MCMC settings below.

The release gates are not only the test suite: `black --check`, `isort --check`,
`ruff check`, `python scripts/check_docs.py` and `mypy phylodater` (0 errors, with
`disallow_untyped_defs` and `warn_return_any` enabled) all have to pass, and both
benchmark variants (reference tree and ML-length tree) have to score 7/7 PASS. The
verbatim output of one such pass is archived at
[`test-data/benchmark/results/validation_log.txt`](../test-data/benchmark/results/validation_log.txt).

## 6. Engine parameters used by the benchmark

Method-specific parameters are reduced from the tool defaults so that the matrix
finishes in minutes. Reductions are recorded here and in the scorecard, never
hidden:

| method | parameters used | default | why |
|--------|-----------------|---------|-----|
| `mcmctree` | `clock=2,num_runs=2,burnin=2000,nsample=4000,sampfreq=10` | `clock=2,num_runs=2,burnin=20000,nsample=50000,sampfreq=5` | 13 taxa/3 kb needs no 50 k samples; convergence is *measured*, not assumed (PSRF ≤ 1.05, ESS ≥ 100) |
| `mdcat` | `ncat=10,nrep=20,max_iter=50` | `ncat=50,nrep=100,max_iter=100` | EM restarts dominate the runtime. **Measured:** the reduced settings reproduce the tool-default estimates exactly on this dataset (all 12 nodes identical, MAE 0.081 Ma both ways, 10 s vs 46 s), so the reduction costs no accuracy here |
| `lsd2`, `r8s`, `treepl`, `pathd8`, `wlogdate` | PhyloDater defaults | — | no reduction needed |

`mdcat.ci_nboots` stays at its default `0` (confidence intervals off), matching
upstream: a *single* bootstrap rep took minutes on this dataset, so a default of
100 would have made the method unusable.

PATHd8 cross-check: the binary compiled from upstream source (Dockerfile recipe,
`gcc -O3 PATHd8.c -lm`) produced a dated tree identical to the reference binary
already in the environment, which is how the engine itself — not just PhyloDater's
use of it — is validated.

## 7. Determinism

Two independent benchmark runs with the same `--seed` were compared node by node:

* `lsd2`, `r8s`, `treepl`, `pathd8`, `wlogdate`, `mdcat`: **bit-identical** node
  estimates across runs.
* `mcmctree`: differs between runs by ≤ 0.3 Ma (≤ 0.3 % of the root age). This is
  expected behaviour for a Bayesian sampler; a fixed PRNG seed does not make PAML's
  MCMC bit-reproducible. The correct quality gate for MCMCTree is therefore the
  convergence diagnostic that PhyloDater records (`max PSRF`, `min ESS`), not
  run-to-run equality.

## 8. Semantic-degradation reporting

Constraint translation is where dating wrappers do real damage silently, so every
degradation must (a) be visible in the log, (b) be recorded in the result, and
(c) be covered by a unit test. Examples that the validation suite checks:

| input constraint | engine | what happens | warning |
|------------------|--------|--------------|---------|
| `soft_lower 60` | lsd2 | hard one-sided bound `NA:-60` | `SemanticDegradationWarning` |
| `maximum 130` | wlogdate / mdcat | point estimate at 65 Ma (midpoint of `[0,130]`) | `SemanticDegradationWarning`, incl. "the upper bound itself is NOT enforced" |
| `uniform 8–10` | mdcat | midpoint 9 Ma | `SemanticDegradationWarning` |
| `gamma` / `skew_normal` / `skew_t` | r8s, pathd8, treepl, lsd2 | normal-approximation interval | `SemanticDegradationWarning` naming the approximation |
| any soft bound | mcmctree | exact `B()/RN()/ST()` prior | none (faithful) |

A method also fails closed when calibrations are *lost* rather than degraded: the
MCMCTree adapter compares the constraints it asked for with the ones actually
present in `mcmctree.tree` and aborts instead of running a silently under-constrained
analysis.

## 9. Regression tests attached to fixed defects

Each defect found by this test work is pinned by a test, so the fix cannot rot:

| defect | symptom before the fix | pinned by |
|--------|------------------------|-----------|
| LSD2 substitution-model default was an amino-acid model (`LG+G`) for a nucleotide workflow | IQ-TREE 2 aborts with `ERROR: File not found LG`, the method never runs | `TestLSD2ModelResolution` (5 cases) |
| MD-Cat `--CI` was passed as a comma-joined token | upstream `int("100,0.025,0.975")` raises `ValueError`; default CI made the method fail every run | `TestMDCatCliContract` |
| MD-Cat CI defaulted to 100 bootstrap replicates | hours of implicit runtime | `TestMDCatConfig::test_ci_parameters` |
| PAML version gate required an exact string (`4.10.8`) | bioconda's `paml` 4.10.10 — the documented install route — was refused outright | version-series compatibility tests in `tests/unit/adapters/test_mcmctree_method.py` |
| `runtime_metadata.json` never recorded engine versions | shipped results had no engine provenance | adapter version recording + benchmark scorecard |
| three tests asserted on the *process-wide* `sys.modules` | they failed only because of test order (a visualisation test imports ete3), not because of product behaviour | `module_imports_cleanly` fixture (subprocess probe) |
| raw MCMCTree smoke test wrote `BDparas = 1 1 0.1` | PAML ≥ 4.10 aborts: "BDparas: expect flag for birth-death process prior" | fixed ctl string in `tests/integration/full_pipeline/test_all_dating_tools.py` |

## 10. Known limitations (accepted, documented, not hidden)

* **Docker image covers 3 of 7 engines** (MCMCTree, LSD2, PATHd8). See
  [EXTERNAL_TOOLS.md](EXTERNAL_TOOLS.md) §7.
* **`r8s` runs through `pyr8s`, NPRS only.** `method=PL` is ignored with an
  `EnvironmentWarning`.
* **treePL has no conda/PyPI package today**, so its installation is environment
  specific; the tests skip cleanly when the binary is absent.
* **CI intervals are engine-specific.** `comparison_plot.png`/`comparison_table.tsv`
  mark cells without an interval rather than drawing them as zero-length bars.
* **Accuracy numbers are dataset-specific.** The MAE values in
  [TEST_RESULTS.md](TEST_RESULTS.md) describe this 13-taxon, 3 kb, strict-clock
  simulation. They are not a general accuracy claim about the methods.

---

*See also: [TEST_RESULTS.md](TEST_RESULTS.md), [EXTERNAL_TOOLS.md](EXTERNAL_TOOLS.md),
[Chinese version](TESTING_CN.md).*
