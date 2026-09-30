> **Language:** English (canonical) · [中文](EXTERNAL_TOOLS_CN.md)

# External Dating Engines

PhyloDater is a orchestration layer: every `dating` method delegates the actual
estimation to an external engine. This document records **exactly** how each
engine was obtained and verified for the released test results in
[TEST_RESULTS.md](TEST_RESULTS.md), because a dating result is only as
reproducible as the binary that produced it.

All commands below were executed inside the conda environment `phylodater` with
its own `micromamba` (`micromamba 2.9.0`). Package availability was verified
against the live channel metadata on 2026-09-29; where an engine is **no longer**
distributed through conda or PyPI this is stated explicitly rather than
hidden behind a command that fails.

## 1. Engine map

| `--method` | Engine actually executed | Distribution status (2026-09-29) | Verified install route |
|-----------|--------------------------|----------------------------------|------------------------|
| `mcmctree` | PAML `mcmctree` (+ `baseml` for the rate-estimation stage) | on bioconda | `micromamba install -c bioconda paml` (4.10.10) |
| `lsd2` | IQ-TREE 2 (`--date`, least-squares dating) | on bioconda | `micromamba install -c bioconda iqtree` (3.1.3) |
| `wlogdate` | `launch_wLogDate.py` from the `logdate` package | on PyPI **and** bioconda | `pip install wlogdate` (1.0.2) |
| `r8s` | `pyr8s` Python API (NPRS) or the `r8s` binary | `pyr8s`: GitHub only (not PyPI, not bioconda); `r8s`: SourceForge/GitHub only | `pip install git+https://github.com/iTaxoTools/pyr8s.git` |
| `mdcat` | `md_cat.py` from MD-Cat | GitHub only (not PyPI, not bioconda) | `pip install git+https://github.com/uym2/MD-Cat.git` |
| `treepl` | `treePL` binary | **no longer published** on bioconda/conda-forge/PyPI | source build, or materialise the historic build + `patchelf` (see §3) |
| `pathd8` | `PATHd8` binary | **no longer published** on bioconda/conda-forge/PyPI | source build from the official site (see §4) |

> **Channels.** `treepl`, `pathd8` and `r8s` are not resolvable from
> `bioconda`/`conda-forge` today (`micromamba search` returns *no entries*, and
> `https://api.anaconda.org/search?name=<pkg>` has no bioconda hit), so
> `conda install -c bioconda treepl pathd8 r8s` — which appears in older PhyloDater
> documentation and in many pipelines — cannot work. Do not treat that command as
> the supported route.

## 2. Python-level dependencies

```bash
micromamba install -n phylodater -c conda-forge -c bioconda \
    python=3.12 paml iqtree nlopt=2.7.1 patchelf
pip install wlogdate
pip install git+https://github.com/iTaxoTools/pyr8s.git
pip install git+https://github.com/uym2/MD-Cat.git
```

`nlopt=2.7.1` is required only by the `treePL` binary (it links `libnlopt.so.0`);
NLopt ≥ 2.10 ships `libnlopt.so.1` instead, which does not satisfy that binary.

## 3. treePL

`treePL` is the only engine whose install cannot be expressed as one package
command today. Two supported options:

**Option A — build from upstream source (recommended):**

```bash
git clone https://github.com/blackrim/treePL.git
cd treePL/deps && tar xzf nlopt-*.tar.gz   # or use the conda nlopt from §2
cd ../src && ./configure && make && make install
```

treePL needs NLopt for optimisation; ADOL-C is optional (parallel cross-validation).

**Option B — reuse the last conda build and repair its search path:**

The `treepl 2.6.3 (h3fd9d12_0)` package puts the executable at the *prefix root*
with `RPATH = $ORIGIN/lib64:$ORIGIN/lib`, so simply copying it into `<prefix>/bin`
breaks the library lookup:

```
error while loading shared libraries: libnlopt.so.0
```

Fix it by pointing `RPATH` at the environment's `lib` directory:

```bash
patchelf --set-rpath '$ORIGIN/../lib' "$CONDA_PREFIX/bin/treepl"
```

Verified in the test environment:

| item | value |
|------|-------|
| executable | `$CONDA_PREFIX/bin/treepl` |
| SHA-256 (before patchelf) | `50d2894a61c3f2be78c7f5bbfa9986f506496d854754d1b14986eafd6eaf5962` |
| provenance | identical to the `treePL` file recorded in `info/paths.json` of the conda package `treepl-2.6.3-h3fd9d12_0` |
| self-reported banner | `treePL version 1.0` (upstream prints 1.0 even in the 2.6.3 build) |
| required shared library | `libnlopt.so.0` (nlopt 2.7.1) |

PhyloDater calls the binary as `treepl` (`software.treepl.treepl_bin`, overridable
with `--method-args treepl_bin=/abs/path/treepl`).

## 4. PATHd8

PATHd8 ships as portable C source from the official site
<https://www2.math.su.se/PATHd8/>; the whole program is one translation unit that
includes the files in `headers/`:

```bash
wget https://www2.math.su.se/PATHd8/PATHd8.zip && unzip PATHd8.zip
cc PATHd8.c -lm -O3 -o PATHd8
install -m 0755 PATHd8 "$CONDA_PREFIX/bin/"
```

This was verified on gcc 11 (Ubuntu 22.04): the compilation is warning-free and the
built binary reproduces the reference binary's estimates exactly on the benchmark
dataset (identical dated tree, see [TESTING.md](TESTING.md) §6).

PATHd8 prints no version banner and accepts no `--version`, so PhyloDater records
the resolved executable path instead; pin the version yourself if you need strict
provenance.

## 5. r8s versus pyr8s

PhyloDater prefers `pyr8s` when the Python package is importable and falls back to
the `r8s` binary otherwise. Consequences that must be understood before quoting
results:

* `pyr8s` implements **NPRS only**. When it is the active backend, the configured
  `method` (`PL`, `LNL`, …) is ignored and an `EnvironmentWarning` is logged.
* The original `r8s` (Sanderson 2003) is distributed from
  <https://sourceforge.net/projects/r8s/> and mirrored at
  <https://github.com/R8S/r8s>. The URL previously printed in this repository
  (`phymaster.massey.ac.nz`) does not resolve.
* `pyr8s` is not published on PyPI or on any conda channel; install it from source.

## 6. What PhyloDater records automatically

For every run, `runtime_metadata.json` now contains the engine identification that
`validate_environment()` discovered:

```json
"software_versions": {"mcmctree": "MCMCTREE in paml version 4.10.10, 29 Jan 2026"}
```

Examples produced by the shipped benchmark:

| method | recorded version string |
|--------|-------------------------|
| `mcmctree` | `MCMCTREE in paml version 4.10.10, 29 Jan 2026` |
| `lsd2` | `IQ-TREE version 3.1.3 for Linux x86 64-bit built Jul 26 2026` |
| `r8s` | `pyr8s 0.3.1 (NPRS only)` |
| `treepl` | `treePL version 1.0` |
| `wlogdate` | `wLogDate 1.0.2` |
| `mdcat` | `1.0.1` |
| `pathd8` | `PATHd8 (no version banner): /…/bin/PATHd8` |

## 7. Docker image coverage

[`Dockerfile`](../Dockerfile) does **not** provide all seven engines. It builds
from upstream tarballs:

| engine | in the image? | source used by the Dockerfile |
|--------|---------------|-------------------------------|
| MCMCTree | yes | PAML 4.10.8 release tarball from GitHub (`abacus-gene/paml`) |
| LSD2 | yes | IQ-TREE 2.3.4 Linux tarball |
| PATHd8 | yes | `PATHd8.zip` from `www2.math.su.se`, compiled with `gcc -O3` |
| treePL | no | — |
| r8s / pyr8s | no | — |
| wLogDate | no | — |
| MD-Cat | no | — |

The image therefore supports 3 of the 7 methods. To use the others inside a
container, add the routes from §2 and §5 (`apt` packages `libnlopt0` and
`libnlopt-dev` are already installed for a future treePL build). Update this
table whenever the Dockerfile changes.

---

*See also: [TESTING.md](TESTING.md) (how the engines are validated),
[TEST_RESULTS.md](TEST_RESULTS.md) (measured results),
[Chinese version](EXTERNAL_TOOLS_CN.md).*
