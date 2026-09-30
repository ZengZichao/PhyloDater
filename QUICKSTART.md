> **Language:** English (canonical) · [中文](QUICKSTART_CN.md)

# PhyloDater Quick Start Guide

> **Last Updated**: 2026-07-12

This guide will help you run your first PhyloDater analysis in 5 minutes.

## 1. Installation (1 minute)

### Using pip

```bash
pip install phylodater
```

### Installing External Software (using conda)

```bash
conda create -n phylodater python=3.12
conda activate phylodater
conda install -c bioconda -c conda-forge paml iqtree   # MCMCTree + LSD2
pip install phylodater
pip install wlogdate
pip install git+https://github.com/iTaxoTools/pyr8s.git     # r8s backend (NPRS)
pip install git+https://github.com/uym2/MD-Cat.git          # MD-Cat
# treePL and PATHd8 must be built from source: they are no longer
# published on bioconda/conda-forge/PyPI. See docs/EXTERNAL_TOOLS.md §3-§4.
```

> `conda install -c bioconda r8s treepl pathd8 pyr8s` does **not** work today: none
> of those four packages exists on the channel any more. Do not copy that line into
> a pipeline or Dockerfile.

### Verify Installation

```bash
phylodater check
```

## 2. Prepare Input Files (2 minutes)

### Example Files

You can use the example data in the `examples/` directory:

```
examples/
├── tree.nwk           # Phylogenetic tree
├── alignment.fasta    # Sequence alignment
└── calibrations.yaml  # Calibration configuration
```

> **Input tree must be rooted**: unrooted trees are rejected during validation. Root your tree (e.g., on an outgroup) before running.

### Minimal Configuration File

Example **calibrations.yaml**:

```yaml
calibrations:
  - name: "Primates"
    mrca_pair: ["human", "cat"]
    constraint:
      type: "uniform"
      min: 70.0
      max: 90.0

  - name: "Rodents"
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

> **Note**: Use `mrca_pair` to specify MRCA by two leaf node names, or `is_root: true` for root node calibration.

## 3. Run Analysis (1 minute)

### Basic Command (subcommand architecture)

```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o my_results/ \
    --method pathd8
```

### Using MCMCTree with Convergence Diagnostics

```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o my_results/ \
    --method mcmctree \
    --method-args "num_runs=4"
```

### Rerooting Options

```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o my_results/ \
    --method pathd8 \
    --reroot midpoint
```

Available strategies: `outgroup`, `midpoint`, `mad`, `none`.

### Using PAML 4.8

```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o my_results/ \
    --method mcmctree \
    --method-args "paml_path=/opt/paml-4.8,paml_version=4.8"
```

### Specifying Software Paths

PAML version and binary paths are passed via `--method-args`, not as top-level CLI options.

Because PAML 4.8 and 4.10.8 binaries share the same names (`mcmctree`, `codeml`, `baseml`), install each version in a separate directory and set `paml_path`:

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --method mcmctree \
    --method-args "paml_path=/opt/paml-4.10.8,paml_version=4.10.8"
```

PhyloDater parses the `mcmctree` banner version before execution and compares the **release series** with `paml_version`: a detected 4.10.10 satisfies the default `4.10.8` (control-file semantics are identical within the 4.10 series), while a 4.8 binary against a 4.10 expectation is refused with a clear error, because those two really do differ. When `paml_path` is set, PhyloDater also temporarily prepends `<paml_path>/bin` to `PATH` so that the `baseml`/`codeml` binaries spawned internally by `mcmctree` come from the same PAML version. If a same-named binary already exists in `PATH`, force a specific executable with `mcmctree_bin=/opt/paml-4.10.8/bin/mcmctree`.

## 4. View Results (1 minute)

After running, check the output directory:

```bash
ls my_results/
```

Main output files:

- `summary.txt` - Execution summary
- `comparison_table.tsv` - Node age comparison table (TSV)
- `comparison_table.txt` - Node age comparison table (text)
- `pathd8/timetree.nex` - Dated tree (PATHd8 NEXUS format)
- `mcmctree/run1/FigTree.tre` - Dated tree (MCMCTree)
- `runtime_metadata.json` - Runtime metadata

## Common Command Examples

### Run Multiple Methods for Comparison

```bash
# Run different methods separately
phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ --method lsd2
phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ --method treepl
phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ --method pathd8
```

> Each method writes to its own subdirectory under `results/`.

### Specify Substitution Model (LSD2)

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --method lsd2 \
    --method-args "model=WAG+G"
```

### Specify Initial Smoothing (treePL)

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --method treepl \
    --method-args "initial_smooth=100"
```

### Use Configuration File

Create `config.yaml`:

```yaml
software_paths:
  paml_path: /path/to/paml
  iqtree_bin: /path/to/iqtree2

software:
  mcmctree:
    clock: 2
    num_runs: 2
    burnin: 10000
    nsample: 20000

  lsd2:
    model: GTR+G   # unset = chosen from the alignment (GTR+G nucleotide, LG+G protein)
```

Run:

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --config config.yaml \
    --method mcmctree
```

## Supported Constraint Types

PhyloDater supports multiple calibration types:

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

> **Note**: MCMCTree's `RootAge` can only be set as an upper bound, not a true fixed root age. If the root uses `fixed`, PhyloDater will raise an error before execution. Use `uniform` or `maximum` for the MCMCTree root, or use treePL/PATHd8/LSD2 if an exact fixed root is required. Different methods may use different calibration files.

## Sequence Formats

FASTA and FASTQ are supported for analysis. PHYLIP, Stockholm, and Clustal formats can be detected but are not supported for analysis; please convert to FASTA first.

## Next Steps

- Read [README.md](README.md) for full features
- Check the `examples/` directory for sample data
- Refer to `docs/api/` for API usage

## FAQ

### Q: External software not found?

Ensure external tools are in your PATH, or specify paths via `--method-args` (e.g., `paml_path=/path/to/paml`, `iqtree_bin=/path/to/iqtree2`).

### Q: Calibration point format?

PhyloDater supports multiple types: `fixed`, `uniform`, `soft_lower`, `maximum`, `soft_bounds`, `gamma`, `skew_normal`, `skew_t`. Use `min`/`max` for `uniform` constraints.

### Q: Which sequence formats are supported?

FASTA and FASTQ are supported for analysis. PHYLIP, Stockholm, and Clustal can be detected but should be converted to FASTA.

### Q: How to check MCMC convergence?

Run MCMCTree with `--method-args "num_runs=4"`. The software automatically calculates PSRF and ESS for convergence diagnostics.

### Q: How to reroot the tree?

Use `--reroot` with one of: `outgroup`, `midpoint`, `mad`, `none`.

### Q: How to specify PAML version?

Use `--method-args "paml_version=4.8"`. This affects MCMCTree control file format and rate estimation workflow.
