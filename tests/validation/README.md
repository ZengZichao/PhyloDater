> **Language:** English (canonical) · [中文](README_CN.md)

# Release-Validation Test Suite

This package validates PhyloDater the way a reviewer would: by running the real
dating engines end to end and checking the numbers they return, not only the
Python plumbing around them.

Everything here is confined to this directory: tests never modify files outside
`tests/validation/` (scratch output goes to `tests/validation/_tmp/`, which is
git-ignored).

## Layout

```
tests/validation/
├── conftest.py     shared fixtures, custom markers, scratch-directory handling
├── helpers.py      engine-availability probes, CLI runner, fixture builders
├── unit/           calibration models/loading/resolution, tree model & topology,
│                   checkpointing, configuration, logging, method registry
├── functional/     CLI subcommands (dating / check / fold), dry runs, error
│                   handling, auto-calibration
├── adapters/       adapter contracts, MCMCTree adapter, other adapters'
│                   environment detection
├── integration/    end-to-end runs of every installed engine, multi-method
│                   scheduling, checkpoint/resume
├── security/       output confinement, path traversal, YAML and input safety
└── fixtures/       static inputs: trees, alignments, calibrations, configs,
                    taxonomies, benchmark
```

## How this relates to the rest of `tests/`

| tree | purpose | external software |
|------|---------|-------------------|
| `tests/unit`, `tests/functional`, `tests/integration` | developer suite: fast, mostly mocked | only for raw-engine smoke tests |
| `tests/validation` | release validation: real engines, end-to-end | required; tests skip with a reason when absent |
| `test-data/benchmark` + `scripts/run_reference_benchmark.py` | **accuracy** against simulated ground truth | required |

The measured outcome of this suite is documented in
[docs/TEST_RESULTS.md](../../docs/TEST_RESULTS.md) — there is deliberately no
second, easily-stale results file inside this directory.

## Running

```bash
python -m pytest tests/validation -v          # everything, engines as available
python -m pytest tests/validation/unit -v     # one category
python -m pytest tests/validation -m slow     # only tests that run real engines
python -m pytest tests/validation -m "not slow"
python -m pytest tests/validation -m security
PHYLODATER_KEEP_TMP=1 python -m pytest tests/validation/integration -v   # keep scratch dirs
```

`python -m pytest tests` (the default `testpaths` in `pyproject.toml`) already
includes this package.

## Markers

| marker | meaning |
|--------|---------|
| `slow` | runs a real external engine (seconds to minutes) |
| `security` | path safety, injection, input-validation tests |
| `requires_mcmctree` … `requires_mdcat` | one per dating engine; skipped with an explicit reason when the engine is missing |

Availability is probed by `helpers.external_available()`, which knows that a
method may be backed by a Python package rather than a binary:

| method | what must exist |
|--------|-----------------|
| `mcmctree` | `mcmctree` (+ `baseml`) |
| `lsd2` | `iqtree2` or `iqtree` (LSD2 is built into IQ-TREE 2) |
| `r8s` | `pyr8s` package **or** `r8s` binary |
| `wlogdate` | `logdate` package **or** `launch_wLogDate.py` |
| `mdcat` | `emd` package **or** `md_cat.py` / `mdcat` |
| `pathd8` | `PATHd8` binary |
| `treepl` | `treePL` binary |

Engine installation and provenance: [docs/EXTERNAL_TOOLS.md](../../docs/EXTERNAL_TOOLS.md).
Methodology and acceptance gates: [docs/TESTING.md](../../docs/TESTING.md).
