> **Language:** English (canonical) · [中文](README_CN.md)

# Documentation

## Documentation language policy

Every document in this repository exists in **two** versions:

| rule | meaning |
|------|---------|
| English is canonical | the English text lives at the bare file name (`README.md`, `docs/TESTING.md`) and is the version the project treats as authoritative |
| Chinese is a translation twin | the Chinese text uses the `_CN` suffix (`README_CN.md`, `docs/TESTING_CN.md`) |
| English is displayed first | the language switch line at the top of every document lists English before Chinese, in both files |
| Both versions change together | a change to behaviour that is documented must be applied to the English file **and** its `_CN` twin in the same commit |
| Machine-readable metadata | `CITATION.cff` is English-only by format; `CITATION_CN.md` is its Chinese companion |
| Issue templates | `.github/ISSUE_TEMPLATE/*.md` stay one file each, with English headings and Chinese guidance under every heading — GitHub lists one entry per file, so a second `_CN` file would duplicate the picker |

## User documentation

| document | content |
|----------|---------|
| [../README.md](../README.md) | feature overview, installation, quick example, method/parameter reference |
| [../QUICKSTART.md](../QUICKSTART.md) | five-minute walkthrough, common workflows, troubleshooting |
| [EXTERNAL_TOOLS.md](EXTERNAL_TOOLS.md) | how each dating engine was obtained, version provenance, source builds, Docker coverage |
| [TESTING.md](TESTING.md) | how the software is tested: layers, ground-truth benchmark design, scoring rules, acceptance gates, reproduction |
| [TEST_RESULTS.md](TEST_RESULTS.md) | measured results for all seven engines, accuracy vs simulated truth, determinism, defects found and fixed |
| [ARCHITECTURE.md](ARCHITECTURE.md) | package layout and module responsibilities |
| [ARCHITECTURE_DECISIONS.md](ARCHITECTURE_DECISIONS.md) | the ADR log: what was decided, why, and what was rejected |
| [ERROR_CODES.md](ERROR_CODES.md) | exit codes and error classes with their meanings |
| [CONDA_PACKAGING.md](CONDA_PACKAGING.md) | packaging PhyloDater (and its engines) for conda |
| [../CONTRIBUTING.md](../CONTRIBUTING.md) | development setup, code style, test and documentation requirements |
| [../SECURITY.md](../SECURITY.md) | vulnerability reporting and scope |
| [../THIRD_PARTY_LICENSES.md](../THIRD_PARTY_LICENSES.md) | dependency and engine licences |

## Data shipped for testing

| path | content |
|------|---------|
| [`../test-data/benchmark/`](../test-data/benchmark/) | ground-truth accuracy benchmark (simulated tree/alignment, true node ages, calibrations, generator, archived results) |
| [`../test-data/input/`](../test-data/input/) | valid / invalid / boundary scenario inputs for manual end-to-end checks |
| [`../examples/`](../examples/) | small example inputs used by the quick start |
| [`../tests/`](../tests/) | automated tests: `unit`, `functional`, `integration`, and the `validation` release suite |

## API reference (Sphinx)

The API reference is generated from docstrings with Sphinx + numpydoc.

```bash
pip install -e ".[dev]"
sphinx-build -b html docs/api docs/api/_build
```

Sources: `api/index.rst` (toctree), `api/conf.py` (configuration), one `.rst` per
package module. Adding a module means adding its `phylodater.<module>.rst` and
listing it in `api/index.rst`.

## Contributing a document

1. Write or update the English file first.
2. Apply the same structural change to the `_CN` file in the same commit.
3. Keep the two-line language switch banner at the top (English listed first).
4. Link new documents from the table above and from the root README so they are
   discoverable from both entry points.

The policy is machine-checked: `python scripts/check_docs.py` fails if an internal
link breaks, a document loses its twin, a `*_EN.md` duplicate reappears, or a
language switch line stops listing English first. CI runs it as a lint step.
