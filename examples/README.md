> **Language:** English (canonical) · [中文](README_CN.md)

# examples

Small, self-consistent inputs for the quick start (`QUICKSTART.md`,
`tutorials/01_quickstart.ipynb`). Ten taxa, an 800 bp alignment, three
calibrations — small enough for a first run to finish in seconds with any engine.

| file | content |
|------|---------|
| `tree.nwk` | rooted binary Newick, branch lengths = expected substitutions/site under a strict clock |
| `alignment.fasta` | 10 × 800 nt alignment simulated under JC69 on that tree |
| `calibrations.yaml` | fixed age (`Homininae` 9 Ma) + uniform range (`Primates` 70–80) + uniform root range (100–110) |
| `taxonomy.tsv` | separate example of the tabular taxonomy input (`--taxonomy-file`), GTDB-style strings |
| `benchmark/` | extra inputs used by `scripts/tier1_manual_comparison.py` and the constraint-degradation experiments |
| `checksums.sha256` | SHA-256 of the files in this directory |

## Provenance

The tree, alignment and calibrations are **generated, not hand-written**:

```bash
python test-data/benchmark/generate_benchmark_dataset.py \
    --profile mini --out examples --force
python scripts/generate_checksums.py --recursive examples
```

They come from the same simulator that produces the ground-truth accuracy dataset
in [`../test-data/benchmark/`](../test-data/benchmark/), so the ages implied by
the branch lengths are known (root 105 Ma, `Homininae` 9 Ma, `Primates` 75 Ma).
`taxonomy.tsv` is unrelated to the tree above: it only demonstrates the
`--taxonomy-file` input format, whose taxa names are accession-style on purpose.

Two properties matter for a first run and are guaranteed by the generator:

* **the tree is binary** — a polytomy makes several engines refuse or degrade;
* **at least one calibration is a fixed age** — PATHd8 rejects a calibration set
  without any `fixage` point, and MD-Cat/wLogDate can only express point
  calibrations at all.

## Verified

All seven methods complete successfully on these inputs in the reference
environment (see [docs/TEST_RESULTS.md](../docs/TEST_RESULTS.md) §1.1 for the
engines and their provenance):

```bash
for m in mcmctree lsd2 r8s treepl pathd8 wlogdate mdcat; do
  phylodater dating -t examples/tree.nwk -s examples/alignment.fasta \
      -c examples/calibrations.yaml -o out_$m --method "$m" --seed 42
done
```

These files are examples, **not** the accuracy reference: quantitative scoring uses
[`../test-data/benchmark/`](../test-data/benchmark/), which ships the true node ages.
