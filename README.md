> **Language:** English (canonical) · [中文](README_CN.md)

# PhyloDater

**A Multi-Software Parallel Platform for Phylogenetic Molecular Dating**

[![Python Version](https://img.shields.io/badge/python-3.9%20--%203.14-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> **Version**: 0.1.0 | **License**: MIT License

[中文文档](README_CN.md)

## Overview

PhyloDater is a unified platform for phylogenetic molecular dating that integrates multiple mainstream dating tools, providing a consistent interface and standardized result output.

> **Note**: This software processes phylogenetic trees (not genomic coordinates) and does not involve 0-based/1-based coordinate system issues.

### Standards and Formats

- **Tree formats**: [Newick](https://en.wikipedia.org/wiki/Newick_format) (`.nwk`), [Nexus](https://en.wikipedia.org/wiki/Nexus_file) (`.nex`)
- **Sequence formats**: [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)
- **Configuration**: YAML 1.2

### Key Features

- **Multi-method support**: MCMCTree, LSD2, r8s, treePL, PATHd8, wLogDate, MD-Cat
- **Subcommand architecture**: `dating` (divergence dating), `check` (self-test), `fold` (tree folding)
- **Taxonomy parsing**: GTDB embedded format and tabular semicolon format
- **Deep input validation**: bracket balance, negative branch lengths, duplicate node names, malicious character detection
- **Real-time logging**: ISO 8601 timestamps, streaming output, file persistence
- **Monophyly detection**: Automatic identification of LUCA/LACA/LBCA and arbitrary taxa

### Supported Dating Methods

| Method | Description | `--method-args` Parameters |
|--------|-------------|---------------------------|
| **MCMCTree** | Bayesian dating from PAML package | `clock`, `num_runs`, `burnin`, `nsample`, `sampfreq`, `skip_rate_estimation`, `paml_version`, `paml_path`, `mcmctree_bin` |
| **LSD2** | Fast dating in IQ-TREE2 | `model`, `date_ci`, `clock_sd`, `date_options`, `date_outlier` |
| **r8s** | Classic molecular clock tool | `method`, `smoothing`, `algorithm` |
| **treePL** | Penalized likelihood dating | `initial_smooth`, `cvstart`, `cvstop`, `cviter`, `cvmultstep` |
| **PATHd8** | Molecular clock calibration | (no extra parameters) |
| **wLogDate** | Log-date dating | `backward_time`, `num_replicates`, `root_time`, `leaf_time`, `sequence_length`, `max_iter`, `pseudocount`, `zero_branch_len` |
| **MD-Cat** | Matrix calibration dates | `ncat`, `nrep`, `max_iter`, `backward_time`, `as_date`, `use_direct_import`, `ci_nboots`, `ci_plower`, `ci_pupper`, `annotate_level` |

### Analysis Workflow

```
Input Files              Processing & Validation                                        Analysis                Output
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Tree file (.nwk)  ─┐
                    ├─→ QC Processing(optional) ─→ Rerooting ─→ Input Validation ─→ Calibration ─┐
Alignment (.fa)   ─┤                                                    │                   │                   │
                    │                                                    ↓                   ↓                   ↓
Calibration (.yaml)┘                                                 Deep Validation       MRCA定位          Dating Method
                                                                       Bracket balance        Monophyly          MCMCTree
                                                                       Negative branches      Special IDs        treePL
                                                                       Duplicate names        Taxonomy           r8s/PATHd8/...
                                                                       Cross-validation                              │
                                                                                                                      ↓
                                                                                                              Results & Report
```

#### Method-Specific Parameter Defaults

| Method | Parameter | Default | Description |
|--------|-----------|---------|-------------|
| MCMCTree | `clock` | `2` | 1=strict, 2=independent relaxed, 3=autocorrelated relaxed |
| | `num_runs` | `2` | MCMC parallel runs (for convergence testing) |
| | `burnin` | `20000` | MCMC burn-in generations |
| | `nsample` | `50000` | MCMC samples |
| | `sampfreq` | `5` | MCMC sampling frequency |
| | `skip_rate_estimation` | `false` | Skip stage-0 rate estimation (`true` = skip and use the default rgene_gamma prior); see "PAML 4.8 Support" for version-specific behavior |
| | `paml_version` | `4.10.8` | Expected PAML **release series**; a detected 4.10.x passes, a 4.8 binary does not |
| | `paml_path` | - | PAML installation directory (recommended); uses `<paml_path>/bin/mcmctree` and isolates matching `baseml`/`codeml` |
| | `mcmctree_bin` | `mcmctree` | MCMCTree executable path or command name |
| r8s | `method` | `NPRS` | Algorithm: `NPRS` (default, works with both backends) / `PL` / `LF` (native `r8s` binary only) |
| | `smoothing` | `None` | Smoothing parameter (None = not written to NEXUS) |
| | `algorithm` | `TN` | Tree search algorithm |
| | `r8s_bin` | `r8s` | Executable name/path for r8s or pyr8s |
| treePL | `initial_smooth` | `100.0` | Initial smoothing parameter |
| | `cvstart` | `1000.0` | Cross-validation start value |
| | `cvstop` | `0.1` | Cross-validation stop value |
| | `cviter` | `10` | CV iteration count |
| | `cvmultstep` | `0.1` | CV multiplier step |
| LSD2 | `model` | _auto_ | Substitution model; resolved from the alignment type when unset (`GTR+G` for nucleotides, `LG+G` for proteins). Set explicitly to override, e.g. `HKY+G` |
| | `date_ci` | `100` | Date confidence interval |
| | `clock_sd` | `0.2` | Clock rate standard deviation |
| MD-Cat | `ncat` | `50` | Number of rate categories |
| | `nrep` | `100` | Random initialization count |
| | `max_iter` | `100` | Maximum EM iterations |
| | `ci_nboots` | `0` | Bootstrap replicates for confidence intervals; `0` = off (as upstream). Each replicate is expensive — enable deliberately |

wLogDate defaults: `backward_time=True`, `num_replicates=1`, `max_iter=50000`, `pseudocount=0.01`.

## Installation

### Via pip

```bash
pip install phylodater
```

### From Source

```bash
git clone https://github.com/ZengZichao/phylodater.git
cd phylodater
pip install -e .
```

### Prerequisites

Every dating method delegates to an external engine, which must be installed
separately (it is not a pip dependency). The table below was verified against the
live conda/PyPI indexes on 2026-09-29; where no package exists, the source route is
listed instead of a command that would fail.

| Software | Method | Where to get it | Installable via |
|----------|--------|-----------------|-----------------|
| PAML (`mcmctree`, `baseml`) | MCMCTree | https://abacus.gene.ucl.ac.uk/software/paml.html | bioconda `paml` |
| IQ-TREE 2 (LSD2 built in) | LSD2 | https://github.com/iqtree/iqtree2 | bioconda `iqtree` |
| wLogDate | wLogDate | https://pypi.org/project/wlogdate/ | PyPI `wlogdate` or bioconda `wlogdate` |
| pyr8s | r8s | https://github.com/iTaxoTools/pyr8s | source (`pip install git+…`) — not on PyPI/conda |
| r8s (original) | r8s | https://sourceforge.net/projects/r8s/ | source build — not on PyPI/conda |
| treePL | treePL | https://github.com/blackrim/treePL | **source build** — no longer published on bioconda/conda-forge/PyPI |
| PATHd8 | PATHd8 | https://www2.math.su.se/PATHd8/ | **source build** (`cc PATHd8.c -lm -O3`) — no longer published on conda/PyPI |
| MD-Cat | MD-Cat | https://github.com/uym2/MD-Cat | source (`pip install git+…`) — not on PyPI/conda |

Per-engine versions measured in the release tests, checksums and exact build
commands: [docs/EXTERNAL_TOOLS.md](docs/EXTERNAL_TOOLS.md).

> **Note on r8s vs pyr8s**: if the Python package `pyr8s` is importable, PhyloDater
> uses it through its native API — which implements **NPRS only**, so a configured
> `method: PL` is ignored with a warning, and `CONSTRAIN`-style ranges are advisory
> rather than enforced (measured; see [docs/TEST_RESULTS.md](docs/TEST_RESULTS.md) §2.2).
> Without `pyr8s`, PhyloDater falls back to the command-line `r8s` binary given by
> `r8s_bin`.

### Using Conda

```bash
conda create -n phylodater-env python=3.12
conda activate phylodater-env
conda install -c bioconda -c conda-forge paml iqtree   # MCMCTree + LSD2
pip install phylodater
# remaining engines (verified routes, see docs/EXTERNAL_TOOLS.md):
pip install wlogdate
pip install git+https://github.com/iTaxoTools/pyr8s.git
pip install git+https://github.com/uym2/MD-Cat.git
# treePL and PATHd8: build from source (docs/EXTERNAL_TOOLS.md §3-§4)
```

### Using Docker

```bash
# Build image
docker build -t phylodater .

# Run container
docker run -v $(pwd):/data phylodater dating \
    -t /data/tree.nwk \
    -s /data/alignment.fasta \
    -c /data/calibrations.yaml \
    -o /data/results/ \
    --method mcmctree
```

> **What the image actually contains.** The image builds **3 of the 7** engines —
> MCMCTree (PAML), LSD2 (IQ-TREE 2) and PATHd8 — and must not be described as a
> full seven-engine platform. See
> [docs/EXTERNAL_TOOLS.md](docs/EXTERNAL_TOOLS.md) §7 for the per-engine table.
> Asking for a method whose engine is absent does **not** fail silently and does
> **not** fabricate a result: the run logs `<engine> not found` with installation
> hints, then `... environment not available, skipping`, then
> `N methods failed: <method> (Environment not available)` and
> `No results to report`. (Verified 2026-09-30 by running the CLI from a clean
> sdist install with no engines on `PATH`; the same behaviour is what the
> `phylodater check` optional-dependency rows report.) `docker build` itself was
> **not** executed in this round's validation — no Docker/Podman runtime was
> available on the validation host — so treat the image as unverified here rather
> than as tested.

## Quick Start

### 1. Self-Test (Verify Installation)

```bash
phylodater check
```

Example output — the table is emitted by the tool in Chinese and is shown verbatim
(measured on the release environment, 2026-09-30):
```
检查项                                      状态       信息
----------------------------------------------------------------------
依赖: biopython                            [PASS]   v1.88
可选依赖: logdate (wLogDate)                 [PASS]   vunknown
可选依赖: emd (MD-Cat)                       [PASS]   vunknown
可选依赖: pyr8s                              [PASS]   vunknown
可选依赖: ete3                               [PASS]   v3.1.3
...
总计: 26 项检查, 26 通过, 0 失败

✓ 所有检查通过！
```

### 2. Prepare Input Files

Three input files are required:

#### Tree File

Supported formats: `.nwk`, `.newick`, `.tree`, `.tre`, `.nex`, `.nexus`, `.nhx` (annotated Newick)

> **Note**: The input tree must be a **rooted** tree. Unrooted trees are now rejected during validation (`phylodater check` or before dating). If your tree is unrooted, root it first (e.g., on an appropriate outgroup) before running PhyloDater.

#### Alignment File

Primary supported formats: FASTA (`.fasta`, `.fa`, `.fna`, `.faa`), FASTQ (`.fastq`, `.fq`)

Detection-only formats (validation only, not used for analysis): PHYLIP (`.phy`), Clustal (`.aln`), Stockholm (`.sto`)

> **Note**: PHYLIP, Stockholm and other alignment formats are detected but not fully supported for analysis. Please convert to FASTA first.

#### Calibration Configuration File (YAML)

Example calibration configuration (calibrations.yaml):

```yaml
calibrations:
  - name: "Primates"
    node:
      type: "mrca_pair"
      mrca_pair: ["human", "cat"]
    constraint:
      type: "uniform"
      min: 70.0
      max: 90.0

  - name: "Rodents"
    node:
      type: "mrca_pair"
      mrca_pair: ["mouse", "rat"]
    constraint:
      type: "soft_lower"
      min: 60.0

  - name: "Root"
    is_root: true
    constraint:
      type: "maximum"
      max: 150.0
```

> **Node definition types**: `auto` (infer from name/taxonomy), `mrca` (list of taxa), `mrca_pair` (two taxa defining the node).

#### Supported Constraint Types

| Type | Description |
|------|-------------|
| `fixed` | Fixed age |
| `uniform` | Uniform distribution `[min, max]` |
| `soft_lower` | Soft lower bound |
| `maximum` | Upper bound constraint |
| `soft_bounds` | Soft bounds |
| `gamma` | Gamma prior |
| `skew_normal` | Skew-normal distribution |
| `skew_t` | Skew-t distribution |

#### Method-Specific Calibration Notes

- **MCMCTree root node**: MCMCTree's `RootAge` can only be set as an upper bound, not a true fixed root age. If the root uses a `fixed` constraint, PhyloDater will raise an error before execution. Use `uniform` or `maximum` for the MCMCTree root, or use treePL/PATHd8/LSD2 if an exact fixed root is required.
- **Per-method calibration files**: Different dating tools support different constraint types. You can prepare different YAML calibration files for different methods and pass them via `-c/--calibrations`.

### 3. Run Dating Analysis

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --method mcmctree
```

### 4. Use Method-Specific Parameters

```bash
# MCMCTree: Specify PAML 4.8 version (also set paml_path to prevent mixing mcmctree with baseml/codeml from another PAML version)
phylodater dating \
    -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ \
    --method mcmctree \
    --method-args "paml_path=/opt/paml-4.8,paml_version=4.8,num_runs=4"

# treePL: Specify smoothing parameter
phylodater dating \
    -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ \
    --method treepl \
    --method-args "initial_smooth=100"
```

### 5. Rerooting Options

```bash
phylodater dating \
    -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ \
    --method pathd8 \
    --reroot midpoint
```

Available strategies: `outgroup`, `midpoint`, `mad`, `none`

### 6. Auto-Add Calibration Points

```bash
phylodater dating \
    -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ \
    --method mcmctree \
    --auto-calibrate "Cyanobacteriota:2500" "Archaea:3800"
```

## CLI Architecture

```
phylodater [common options] <subcommand> [subcommand options]

Subcommands:
  dating          Run dating analysis
  check (self-test)  Self-test and dependency verification
  fold            Fold tree by taxonomic groups
```

### Common Parameters (Shared by All Subcommands)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--log-level` | Log level (DEBUG/INFO/WARNING/ERROR/CRITICAL) | `INFO` |
| `--log-file` | Log file path (UTF-8, real-time flush) | None |
| `--taxonomy-file` | Taxonomy mapping file | None |
| `--taxonomy-levels` | Custom taxonomy level prefixes (e.g., `k:kingdom ss:subspecies`) | None |
| `--taxonomy-delimiter-mode` | Format A parsing strategy (`reverse`/`greedy`/`segment`) | `reverse` |
| `--taxonomy-source-priority` | Taxonomy source priority (`embedded`/`table`) | `table` |
| `--table-sep` | Force taxonomy table separator | Auto-detect |
| `--ignore-malformed` | Skip malformed input (e.g., unsupported formats with warning) | No |
| `--version` | Show version info (with Git hash, dependency versions) | - |

### dating Subcommand Parameters

| Parameter Group | Parameter | Description | Default |
|----------------|-----------|-------------|---------|
| Input Files | `-t/--tree` | Tree file path | Required |
| | `-s/--sequence` | Sequence file path | Required |
| | `-c/--calibrations` | Calibration config file | Required |
| | `-o/--output` | Output directory | Required |
| | `--config` | Config file path (YAML) | None |
| | `--generate-config` | Generate default config file to given path and exit | None |
| Dating Method | `--method` | Dating method name | Required |
| | `--method-args` | Method-specific parameters (key=value) | None |
| Tree Processing | `--reroot` | Rerooting strategy (`outgroup`/`midpoint`/`mad`/`none`) | `none` |
| | `--outgroup` | Outgroup name | None |
| | `--multi-tree-mode` | Multi-tree mode (`split`/`first`/`last`/`random`/`ask`/`error`) | `error` |
| Calibration | `--auto-calibrate` | Auto-add calibration points | None |
| Execution | `-j/--threads` | Parallel threads | `4` |
| | `--seed` | Global random seed | Auto-generated |
| | `--dry-run` | Dry-run mode (validation only) | No |
| | `--mol-type` | Sequence molecule type (`DNA`/`RNA`/`protein`) | Auto-detect |
| | `--no-cross-check` | Disable tree/sequence cross-validation | No |
| | `--skip-length-check` | Skip negative branch length check | No |
| | `--low-memory` | Low memory mode (single thread, reduced MCMC params) | No |
| Checkpoint | `--resume` | Resume from checkpoint (skip completed steps) | No |
| | `--clear-checkpoint` | Clear checkpoint files before running | No |
| Output Control | `--force` | Overwrite existing output files | No |
| | `--no-clobber` | Skip existing output files | No |

#### fold Subcommand Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-t/--tree` | Input tree file | Required |
| `-o/--output` | Output directory | Required |
| `--fold-rank` | Taxonomic rank to fold to (`domain`/`phylum`/`class`/`order`/`family`/`genus`) | `family` |

> **Note**: The `fold` subcommand is currently a **stub implementation**. It accepts parameters and logs the requested rank, but does not actually perform tree folding and will output a "not yet implemented" warning. Do not use it for real analysis.

> **Custom software paths**: Provide per-tool binary paths via the YAML config `software_paths` section (e.g. `mcmctree_bin`, `iqtree_bin`, `r8s_bin`, `treepl_bin`, `pathd8_bin`) or the environment-specific `PATH`.
>
> **Coexisting PAML versions**: PAML 4.8 and 4.10.8 binaries share the same names (`mcmctree`, `codeml`, `baseml`). To avoid selecting the wrong version:
> 1. Install each version in a separate directory (e.g. `/opt/paml-4.8`, `/opt/paml-4.10.8`) and set `paml_path` in `software_paths` to the target root; the adapter will automatically use `<paml_path>/bin/mcmctree`.
> 2. Alternatively, provide the full path via `mcmctree_bin`, e.g. `/opt/paml-4.10.8/bin/mcmctree`.
> 3. Before execution PhyloDater parses the actual PAML version from the executable banner and compares its **release series** with the configured `paml_version`: any 4.10.x satisfies the default `4.10.8`, while a cross-series mismatch (e.g. a 4.8 binary) is refused with a clear error, preventing control-file mismatches.

## Taxonomy Parsing

PhyloDater supports two taxonomy name formats:

### Format A (Embedded)

```
GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_c_Cyanobacteriia_o_Cyanobacteriales_f_Prochloraceae_g_Prochloron
```

Delimiters: `_d_` (domain), `_p_` (phylum), `_c_` (class), `_o_` (order), `_f_` (family), `_g_` (genus)

Parsing strategies (controlled by `--taxonomy-delimiter-mode`):
- `reverse` (default): right-to-left matching, avoids ambiguity from underscores in taxonomy names
- `greedy`: greedy matching
- `segment`: segment from the first `_d_`

### Format B (Tabular Semicolon)

```
d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__Synechococcales;f__Synechococcaceae;g__Synechococcus;s__
```

**Tolerance rules**:
- Taxonomy values must not contain `;` or `__`, otherwise ERROR and skip the line
- `s__` with no content is parsed as `None`, logged as DEBUG

### External Taxonomy Table

When leaf nodes lack complete taxonomy names, supplement via table file:

```
# Two-column format (tab or comma separated)
GB_GCA_000252485.1	d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;...
```

Use `--taxonomy-file` parameter to specify the table file.

**Label consistency check**: Leaves in the table not found in the tree trigger WARNING; tree leaves without table records and unparseable IDs get "unknown" taxonomy.

**Source priority**: When a node has both embedded and table sources, control priority via `--taxonomy-source-priority` (default `table`). Conflicts trigger WARNING with details.

## Monophyly Detection

PhyloDater supports automatic monophyly detection by specifying taxon names:

### General Taxon Detection

Given a taxonomic label (e.g., "Cyanobacteriota"):
1. Collect all leaf nodes belonging to that taxon
2. Compute the Most Recent Common Ancestor (MRCA) of these leaves
3. Get all descendant leaves of the MRCA
4. If the descendant set equals the step 1 set → monophyletic; otherwise → non-monophyletic

Non-monophyletic groups trigger ERROR and refuse operation, listing extra leaves not belonging to the taxon.

### Special Identifiers

- **LUCA**: MRCA of all Bacteria and Archaea
- **LBCA**: MRCA of all Bacteria
- **LACA**: MRCA of all Archaea

Special identifiers are case-sensitive (uppercase only). Taxon existence is verified before computation.

## PAML 4.8 Support

PhyloDater applies a unified three-stage workflow to both PAML 4.8 and 4.10.8:

| Stage | Function | Description |
|-------|----------|-------------|
| Stage 0 | Rate estimation (baseml/codeml, clock=1) | Enabled by default. Fits a global clock, extracts the overall substitution rate from the "Substitution rate is per time unit" output, and uses it as the `rgene_gamma` prior mean. The result is saved as `mcmctree_rate_estimation.txt` in the output directory. |
| Stage 1 | Hessian matrix (usedata=3) | Generate branch lengths, gradient, Hessian |
| Stage 2 | MCMC sampling (usedata=2) | Estimate divergence times |

### Skipping rate estimation (optional)

To restore the legacy behavior of using the default weakly informative prior
(`rgene_gamma` mean = rate_alpha/20, i.e. the default G(2, 20)), set
`skip_rate_estimation: true`:

```bash
phylodater dating -t tree.nwk -s aln.fasta -c cal.yaml -o out/ \
    --method mcmctree \
    --method-args "skip_rate_estimation=true"
```

It can also be set in the `software.mcmctree` block of a YAML config file.
The default is `false` (rate estimation is performed).

### Version-specific behavior of the rate output

| PAML version | "Substitution rate is per time unit" output | Notes |
|--------------|---------------------------------------------|-------|
| ≤ 4.9h (incl. 4.8) | Printed | Normal behavior |
| 4.9i | **Missing** | A confirmed bug (see PAML changelog, v4.9j entry). PhyloDater does not support 4.9i; on this version rate estimation falls back to the default prior. |
| 4.9j / 4.10.x (incl. 4.10.8) | Printed | The bug was fixed in 4.9j; behavior matches 4.8. The only difference is cosmetic: 4.10.8 prints an extra blank line between the rate heading and the value, which PhyloDater's extraction logic already handles. |

Verified empirically: on identical data, baseml (clock=1) of PAML 4.8 and 4.10.8
produce the same rate estimate (differences only in the 5th decimal place,
e.g. 4.695519 vs 4.695467).

**Note**: the rate line is only printed when the tree contains an absolute age
annotation (`@`, i.e. fossil/root age constraints). PhyloDater automatically
writes a representative root-calibration age (stored in Ma, converted to Ga to
match in-tree calibrations) as `@age` into the rate-estimation tree. If no root
calibration and no `root_age` are configured, rate estimation cannot proceed and
falls back to the default prior with a WARNING.

```bash
phylodater dating -t tree.nwk -s aln.fasta -c cal.yaml -o out/ \
    --method mcmctree \
    --method-args "paml_version=4.8"
```

## MCMC Convergence Diagnostics

For MCMCTree with `num_runs >= 2`, convergence diagnostics are automatically performed:

- **PSRF** (Potential Scale Reduction Factor): Gelman-Rubin statistic, ideal < 1.2
- **ESS** (Effective Sample Size): effective sample size, recommended ≥ 100

Convergence results are logged. PSRF ≥ 2.0 triggers error; 1.2 ≤ PSRF < 2.0 triggers warning.

## Input Validation

### Tree File Validation

- Bracket balance check
- Non-negative branch lengths (negative → CRITICAL termination)
- Duplicate node name detection (tip duplicates → ERROR, internal duplicates → WARNING)
- Malicious character injection detection (control characters, Unicode bidi text)

### Sequence File Validation

- Format auto-detection (FASTA/FASTQ)
- Sequence alphabet verification (DNA/RNA/protein)
- Sequence ID uniqueness check
- Sequence length consistency check

### Cross-Validation

When both tree and sequence files are provided, leaf label bidirectional matching is performed by default. Disable with `--no-cross-check`.

### Adversarial Input Protection

- **Malicious character injection**: detection of control characters and Unicode bidi override characters
- **Taxonomy circular dependencies**: detection of circular references in tables (e.g., `d__A;p__B` and `d__B;p__A`)
- **Empty files**: termination when tree or sequence file is 0 bytes

## Multi-Tree Handling

By default, tree files are expected to contain a single tree. If multiple trees are detected:

| Mode | Behavior |
|------|----------|
| `error` (default) | Raise an error and stop (non-interactive default) |
| `split` | Process each tree separately |
| `first` | Process only the first tree |
| `last` | Process only the last tree |
| `random` | Randomly select one tree |
| `ask` | Interactive selection prompt, then exit with code 2 |

## Logging System

### Log Format

```
2025-03-21T10:15:30.123 | INFO     | module/func Parsing tree file
2025-03-21T10:15:30.456 | WARNING  | module/func Large tree detected (15000 tips)
```

### Log File

Use `--log-file` parameter to write logs to file simultaneously:

```bash
phylodater dating -t tree.nwk -s aln.fasta -c cal.yaml -o out/ \
    --method mcmctree \
    --log-file phylodater.log
```

## Output Structure

After completion, the output directory contains:

```
results/
├── summary.txt             # Execution summary
├── comparison_table.txt    # Text format comparison table
├── comparison_table.tsv    # TSV format comparison table
├── comparison_plot.png     # Comparison plot (generated even for single-method runs; most useful when multiple methods)
├── phylodater.log          # Runtime log
├── runtime_metadata.json   # Runtime metadata snapshot
├── runtime_config.yaml     # Runtime configuration snapshot
├── mcmctree/               # MCMCTree results
│   ├── run1/
│   │   ├── FigTree.tre     # Dated tree (with calibration annotations)
│   │   └── mcmc.txt        # MCMC samples
│   └── run2/
├── lsd2/                   # LSD2 results
├── treepl/                 # treePL results
└── pathd8/                 # PATHd8 results
```

### Calibration Annotation Format

In dated tree files, calibrated nodes are annotated with `[&height=AGE,height_95%_HPD={LOWER,UPPER}]`.

## Configuration File

Use `--config` parameter to specify a YAML configuration file as an alternative to command-line arguments:

```yaml
software_paths:
  # For coexisting PAML versions, specify paml_path and the adapter will use
  # <paml_path>/bin/mcmctree, <paml_path>/bin/codeml, <paml_path>/bin/baseml
  paml_path: /opt/paml-4.10.8
  # To pin a specific version (e.g. when PATH contains a same-named binary):
  # mcmctree_bin: /opt/paml-4.10.8/bin/mcmctree
  iqtree_bin: /path/to/iqtree2
  r8s_bin: /path/to/r8s
  treepl_bin: /path/to/treepl
  pathd8_bin: /path/to/PATHd8

software:
  common:
    nthreads: 8
    seed: 42

  mcmctree:
    clock: 2
    num_runs: 4
    burnin: 20000
    nsample: 50000

  r8s:
    method: pl
    smoothing: 1000

  treepl:
    initial_smooth: 100

  lsd2:
    model: LG+G
```

Configuration priority: defaults < YAML file < command-line arguments. Values are type-checked and validated against each method's schema.

## Exit Codes

| Code | Meaning | Description |
|------|---------|-------------|
| `0` | Success | Analysis completed normally |
| `1` | Runtime error | Missing dependencies, self-test failure, etc. |
| `2` | Usage error | Invalid command-line arguments |
| `3` | Data error | Input file format or content error |
| `130` | Interrupted | Ctrl+C termination |

## Large-Scale Data Handling

For trees with more than 10,000 tips:

- INFO log prompts about required computational resources
- Use `--low-memory` to reduce resource usage (single thread, reduced MCMC parameters)
- Use `--threads` to increase parallel threads
- Recommended at least 16GB memory

## User Interruption & Resource Cleanup

Press `Ctrl+C` for graceful termination:

- Outputs current progress (number of trees processed, etc.)
- Closes all open file handles
- Deletes incomplete temporary outputs
- Exits with code `130`

## Cross-Platform & Encoding

- **Encoding**: All file I/O uses UTF-8, falls back to utf-8-sig on illegal bytes
- **Path handling**: Uses `pathlib.Path`, compatible with Windows/macOS/Linux
- **Line endings**: Auto-normalized on read, unified to `\n` (Unix style) on write
- **Temporary files**: Uses `tempfile` module, ensures cleanup on abnormal exit

## Version Management

- Version follows semantic versioning `MAJOR.MINOR.PATCH`; the release version is pinned statically to `0.1.0` (in `pyproject.toml` and `phylodater/_version.py`) — no Git-derived `.dev`/hash iteration strings are baked into the package
- Inside a Git checkout, `--version` additionally reports the current short commit hash for support purposes
- CI automatically runs `phylodater check` to verify all test cases

`--version` output example:
```
PhyloDater 0.1.0
License: MIT License
Updated: 2026-07-12
git: abc1234
deps: BioPython=1.87, ETE3=3.1.3, DendroPy=5.0.8, NumPy=2.4.6
```

### Tested Environments

The release validation ran on Python 3.12.14 in a conda environment holding all
seven engines; see [docs/TEST_RESULTS.md](docs/TEST_RESULTS.md) §1 for the full
inventory (OS, library versions, engine versions and their provenance). The suite
is portable: `python -m pytest tests` passes with fewer engines installed, because
each engine-dependent test skips with an explicit reason instead of failing.

> **Coverage note**: `pytest tests --cov=phylodater` measures **67.87%** on the
> validated environment (and 67.89% with `-m "not slow"`, so the figure does not
> depend on the engines being installed). The CI gate in `pyproject.toml` is
> `fail_under = 60` — that margin over the measured value absorbs platform
> differences in CI; it is a floor, not a target.

## Example Data

The `examples/` directory contains example files for testing:

```
examples/
├── tree.nwk            # Example tree file (Newick format)
├── alignment.fasta     # Example alignment file (FASTA format)
└── calibrations.yaml   # Example calibration config file
```

> **Data source**: Example data is simulated data for functional testing only, not involving any real biological samples or private data.

Run example:
```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o results/ \
    --method mcmctree
```

## Project Structure

```
phylodater/
├── phylodater/
│   ├── adapters/            # Dating method adapters (one per engine, DatingMethod[ConfigT])
│   ├── cli/                 # Command-line interface (subcommand architecture)
│   ├── core/                # Core interfaces (method_interface, registry), exceptions,
│   │                        #   pipeline, batch processing, comparison reporting
│   ├── infrastructure/      # Logging, configuration, process runner, checkpoints,
│   │                        #   signal handling, plugins, self-test
│   ├── models/              # Data models (tree, constraints, calibration, results)
│   ├── services/            # Taxonomy parsing, calibration loading/resolution,
│   │                        #   validation, rooting, name mapping
│   ├── viz/                 # Visualization (dating figures, tree/geo plots, backends)
│   └── data -> viz/data/    # bundled geological timescale table (CSV package data)
├── examples/                # Mini example dataset (+ truth.yaml, checksums)
├── test-data/
│   └── benchmark/           # Ground-truth reference dataset + generator
│       ├── tree.nwk / tree_ml.nwk / alignment.fasta
│       ├── truth.yaml / calibrations.yaml / checksums.sha256
│       └── results/         # archived scorecard.tsv, node_ages.tsv,
│                            #   benchmark_results.json (+ ml_tree/ variant),
│                            #   PYTHON_SUPPORT_MATRIX.txt
├── tests/
│   ├── unit/                # Fast, mocked
│   ├── functional/          # CLI / pipeline behaviour
│   ├── integration/         # full-pipeline scenarios
│   └── validation/          # real-engine end-to-end validation (markers: slow,
│                            #   requires_mcmctree, requires_lsd2, ...)
├── scripts/
│   ├── run_reference_benchmark.py   # the official, ground-truth-scored benchmark
│   ├── check_docs.py                # bilingual policy + link gate
│   ├── generate_checksums.py        # dataset checksum manifests
│   ├── probe_pyr8s_constraints.py   # pyr8s capability probe (results in docs)
│   ├── tier1_manual_comparison.py   # manual-vs-wrapper experiment
│   └── degradation_experiment.py    # constraint-degradation experiment
├── tutorials/               # Jupyter quick-start (EN + _CN notebook)
└── docs/                    # Documentation (English canonical + _CN Chinese)
    ├── ARCHITECTURE.md / ARCHITECTURE_DECISIONS.md / ERROR_CODES.md
    ├── EXTERNAL_TOOLS.md / TESTING.md / TEST_RESULTS.md / CONDA_PACKAGING.md
```

## Testing and Validation

| what | where |
|------|-------|
| How the software is tested, and why (methodology, gates, reproducing commands) | [docs/TESTING.md](docs/TESTING.md) |
| Measured results: all seven engines, accuracy against simulated ground truth, determinism, defects found and fixed | [docs/TEST_RESULTS.md](docs/TEST_RESULTS.md) |
| Ground-truth dataset + generator + archived results | [`test-data/benchmark/`](test-data/benchmark/) |
| Reproduce the whole benchmark | `python scripts/run_reference_benchmark.py` |
| How each external engine was obtained (versions, checksums, source builds) | [docs/EXTERNAL_TOOLS.md](docs/EXTERNAL_TOOLS.md) |
| Release-validation suite (real engines, end to end) | [`tests/validation/`](tests/validation/) |

```bash
phylodater check                     # 26-item self test, no external engine needed
python -m pytest tests -q            # full automated suite (1489 tests)
python scripts/run_reference_benchmark.py   # seven engines vs ground truth
```

## Dependencies

| Library | Version | License |
|---------|---------|---------|
| BioPython | >=1.79,<2.0 | BSD-like |
| ETE3 | >=3.1.2,<4.0 (optional) | GPLv3 |
| DendroPy | >=4.5.0 | BSD 3-Clause |
| NumPy | >=1.21.0 | BSD 3-Clause |
| Pandas | >=1.3.0 | BSD 3-Clause |
| PyYAML | >=6.0,<7.0 | MIT |
| SciPy | >=1.9.0,<2.0 | BSD 3-Clause |
| Matplotlib | >=3.5.0 | PSF-based |
| psutil | >=5.9.0,<6.0 | BSD 3-Clause |

> **Note**: ETE3 is an optional dependency (declared in `extras_require["ete3"]`). It is installed via pip by default. All other libraries are core dependencies.

## Documentation

- [Quick Start Guide](QUICKSTART.md)
- [Documentation index](docs/README.md)
- [Testing guide](docs/TESTING.md) and [test results](docs/TEST_RESULTS.md)
- [External dating engines](docs/EXTERNAL_TOOLS.md)
- [API Reference](docs/api/)
- [Architecture](docs/ARCHITECTURE.md) · [design decisions](docs/ARCHITECTURE_DECISIONS.md) · [error codes](docs/ERROR_CODES.md)
- [Contributing Guide](CONTRIBUTING.md)

Every document exists in English (canonical, shown first) and Chinese
(`*_CN.md`).

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for how to contribute.

- **Bug Reports**: [GitHub Issues](https://github.com/ZengZichao/phylodater/issues)
- **Feature Requests**: [GitHub Discussions](https://github.com/ZengZichao/phylodater/discussions)
- **Contact**: Zichao Zeng (曾子超) — zengzichao@sjtu.edu.cn · ORCID: [0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

## Citation

If you use PhyloDater in your research, please cite:

> Zeng, Zichao. (2026). PhyloDater: A multi-software parallel platform for phylogenetic molecular dating. GitHub repository. https://github.com/ZengZichao/phylodater (ORCID: 0000-0001-6553-970X)

BibTeX format:
```bibtex
@software{zeng2026phylodater,
  author = {Zeng, Zichao},
  title = {PhyloDater: A multi-software parallel platform for phylogenetic molecular dating},
  year = {2026},
  url = {https://github.com/ZengZichao/phylodater},
  version = {0.1.0},
  note = {ORCID: 0000-0001-6553-970X}
}
```

Full citation guidance (machine-readable CFF, engine citations, and the fact that
no article DOI exists yet) lives in [CITATION.md](CITATION.md).
