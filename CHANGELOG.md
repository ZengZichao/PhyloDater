> **Language:** English (canonical) · [中文](CHANGELOG_CN.md)

# Changelog

All notable changes to **PhyloDater** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html),
as stated in [`README.md`](README.md) ("Version Management") and in
[`pyproject.toml`](pyproject.toml) (`Development Status :: 4 - Beta`).

> **How this file was written.** It was reconstructed on 2026-10-08 from what the
> repository itself documents — [`docs/TEST_RESULTS.md`](docs/TEST_RESULTS.md),
> [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md),
> [`docs/EXTERNAL_TOOLS.md`](docs/EXTERNAL_TOOLS.md),
> [`docs/TESTING.md`](docs/TESTING.md) and the release notes in
> [`README.md`](README.md). The published history of this repository is a single
> squashed commit, so **no interim release dates or version increments are
> claimed here**: everything below describes the 0.1.0 release and the defects its
> own validation campaign records as fixed.

## [0.1.0] - 2026-09-30

First public release. Date from `CITATION.cff` (`date-released: 2026-09-30`,
`version: 0.1.0`); the release-gate campaign behind the numbers below is dated
2026-09-29 (`docs/TEST_RESULTS.md:10`, verdict **PASS**) and its archived command
output was captured 2026-09-30 (`test-data/benchmark/results/validation_log.txt:3`).

### Added

- **Seven engine adapters behind one interface** (`mcmctree`, `lsd2`, `treepl`,
  `pathd8`, `r8s`, `wlogdate`, `mdcat`), each registered through
  `DatingMethodRegistry` and sharing one abstract contract:
  `validate_environment` / `prepare_inputs` / `execute` / `parse_results`
  (`docs/ARCHITECTURE.md`, "Adapter Layer"). `r8s` and its Python port `pyr8s`
  are served by a single adapter registered as `r8s`.
- **Eight calibration constraint types** plus their abstract base
  (`models/constraints.py`): fixed age, uniform, soft lower bound, maximum age,
  soft bounds, gamma prior, skew-normal and skew-t — with per-engine rendering
  (`to_software_format`) and an explicit degradation report whenever a target
  engine cannot represent a prior exactly
  (`docs/ARCHITECTURE.md`, "Constraint Type Hierarchy"; `docs/TESTING.md` §8).
- **Automatic node resolution** for calibration taxa: special identifiers
  (LUCA/LBCA/LACA/ROOT) followed by five progressive matching levels, from exact
  tip names to taxonomy-aware and diversity-maximising candidate selection
  (`services/calibration_resolver.py`).
- **Parallel multi-method pipeline**: thread-pool execution of several methods
  over one input (`core/pipeline.py`), capped at 8 workers and configurable via
  `--threads` / `PHYLODATER_THREADS`; per-method seeds derived deterministically
  from a single `--seed` so stochastic engines are reproducible
  (`docs/ARCHITECTURE.md`, "Data Flow"; `docs/TESTING.md` §7).
- **Checkpointing and resume** for both the pipeline and MCMCTree's two-stage
  workflow (`--resume`, `--clear-checkpoint`), plus graceful `SIGINT` handling
  with resource cleanup (exit code 130).
- **Three CLI subcommands** — `dating`, `check` (the 26-item self-test) and
  `fold` — with documented exit codes 0/1/2/3/130 (`README.md`, "CLI
  Architecture", "Exit Codes").
- **Standardised reporting artefacts**: `comparison_table.txt`,
  `comparison_table.tsv`, `summary.txt`, `comparison_plot.png`, per-method dated
  trees (`{method}_dated_tree.nexus` / `.nhx`), `phylodater.log`, and
  `runtime_config.yaml` / `runtime_metadata.json`. The runtime metadata now
  records the engine version each adapter detected, so a dated tree carries the
  provenance of the binary that produced it (`docs/EXTERNAL_TOOLS.md` §6).
- **Visualisation layer** (`phylodater.viz`): `DatingFigure` for dated trees with
  geological-timescale colour/size mapping, clade highlighting and collapsing,
  plus a pluggable backend registry — `matplotlib` by default and an interactive
  `plotly` backend (opt-in via the `interactive` extra, `pyproject.toml`) — and a
  bundled ICS geological timescale table.
- **Input validation and hardening**: bracket-balance and duplicate-name tree
  checks, sequence alphabet and length consistency, encoding disclosure with
  BOM/GBK repair hints, zero-width and Unicode-bidi character rejection, a
  500 MB input-size guard, path-traversal confinement of outputs, and
  taxonomy circular-dependency detection — all reported with actionable
  messages (`README.md` "Input Validation", `SECURITY.md`,
  `test-data/README.md` "boundary" inputs, `docs/TEST_RESULTS.md` §4.1 and §5).
- **PAML version awareness**: 4.8 and 4.10.x control-file dialects, with
  compatibility decided by the major.minor series and cross-series use refused
  rather than silently downgraded (`README.md` "PAML 4.8 Support").
- **MCMC convergence diagnostics** for multi-run MCMCTree: automatic Gelman–Rubin
  PSRF and ESS computation against the 1.05 / 100 thresholds
  (`README.md` "MCMC Convergence Diagnostics").
- **Distribution and governance**: MIT licence, bilingual (English/Chinese)
  README, QUICKSTART, CONTRIBUTING, CODE_OF_CONDUCT, SECURITY and
  THIRD_PARTY_LICENSES files, Sphinx API sources under `docs/api/`, tutorials,
  `environment.yml` for Conda, a four-stage `Dockerfile`
  (`builder`/`runtime`/`complete`/`tutorial`) and `docker-compose.yml`, plus
  `CITATION.cff` / `CITATION.md`.
- **Test and benchmark infrastructure**: 1,489 collected tests (1,487 passed,
  2 skipped) across unit, functional, integration and release-validation tiers,
  a 26-item `phylodater check` self-test, coverage measured at 67.87% with a
  `fail_under = 60` gate, a 13-taxon / 3 kb reference benchmark with archived
  scorecards (`test-data/benchmark/results/`), and a CI workflow running lint,
  docs-policy, mypy and the test matrix (`docs/TEST_RESULTS.md` §4-§5,
  `.github/workflows/ci.yml`).

### Fixed

Ten defects were found and fixed by the 0.1.0 validation campaign itself; each is
recorded with its symptom, its fix and the regression test now pinning it in
[`docs/TEST_RESULTS.md`](docs/TEST_RESULTS.md) §6:

1. LSD2 defaulted to the amino-acid model `LG+G`, so IQ-TREE 2 could never run on
   DNA; the model is now resolved from the alignment type.
2. MD-Cat's `--CI` argument was passed as one comma-joined token that upstream
   cannot parse; it is now a single whitespace-separated argument.
3. MD-Cat CI defaulted to 100 bootstrap replicates (hours of unannounced work);
   the default is now `ci_nboots = 0`, matching upstream.
4. The PAML gate demanded the literal version `4.10.8` and refused the bioconda
   `paml` 4.10.10; compatibility is now decided by major.minor series.
5. `runtime_metadata.json` recorded no engine versions; every adapter now stores
   what `validate_environment()` detected.
6. The zero-dependency Newick walker rejected a leading `[&R]`, breaking trees
   exported by BEAST/IQ-TREE/DendroPy; leading comments are now skipped and
   normalised.
7. PATHd8 dropped non-point root calibrations (153 Ma returned against a
   requested [100, 110] Ma); root ranges are now rendered as
   `minage`/`maxage`.
8. pyr8s range calibrations were neither enforced nor reported; the adapter now
   discloses that pyr8s ranges are advisory before the run and verifies each
   requested range afterwards.
9. Three unit tests asserted on the process-wide `sys.modules` and failed purely
   because of test order; they now probe a fresh interpreter.
10. The raw MCMCTree smoke test wrote a `BDparas` line PAML ≥ 4.10 rejects; the
    control file now includes the construction flag.

### Known limitations (documented, not hidden)

- The Docker `complete` stage ships **3 of the 7 engines** (MCMCTree via PAML
  4.10.8, LSD2 via IQ-TREE 2.3.4, PATHd8); treePL, r8s/pyr8s, wLogDate and MD-Cat
  must be added by hand (`Dockerfile`, `docs/EXTERNAL_TOOLS.md` §7).
- `treePL` prints the banner `treePL version 1.0` even in the 2.6.3 build, and
  PATHd8 prints no version at all, so engine provenance also relies on the
  recorded executable path (`docs/TEST_RESULTS.md` §7).
- Cross-platform engine performance was not measured; the recorded wall times and
  accuracies belong to the single Linux host of `docs/TEST_RESULTS.md` §1 and to
  the reference dataset, and must not be quoted as general claims.
- `pyr8s` implements NPRS only and treats range calibrations as advisory
  (`docs/EXTERNAL_TOOLS.md` §5).
