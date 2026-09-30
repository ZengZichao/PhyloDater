#!/usr/bin/env python3
"""Semantic-degradation experiment: MCMCTree node ages under gamma priors vs
their 95%-CI uniform degradations (paper Fig. 5B / abstract claim).

Design:
- A 10 kb JC69 alignment with a strict clock is simulated along the six-hominoid
  benchmark tree (fixed seed), so that MCMCTree runs are well-behaved.
- Two calibration files are generated for the same tree:
    gamma   : HumanChimp G(alpha=50, beta=50/6), GreatApes G(alpha=50, beta=50/9),
              Root maximum 15 Ma;
    uniform : the same nodes with U[lower, upper] where the bounds are the
              95% credible interval of each gamma (scipy.stats.gamma.ppf(0.025/0.975,
              alpha, scale=1/beta) * scale + offset), i.e. exactly the
              degradation rule implemented in phylodater.models.constraints;
              Root unchanged.
- MCMCTree is run through the PhyloDater CLI once per (calibration file, seed)
  combination, seeds 42 and 43, num_runs=1.
- Node-age deviations (uniform - gamma) / gamma are computed per node and seed.

Tool invocation (run from examples/benchmark/degradation/ after `prepare`;
PY is the phylodater conda environment python):

    $PY -m phylodater dating -t tree.nwk -s alignment.fasta \
        -c calibrations_gamma.yaml -o gamma_seed42 --method mcmctree \
        --seed 42 --method-args "clock=2,num_runs=1,burnin=2000,sampfreq=2,nsample=5000"
    (repeat with seed 43, and with calibrations_uniform.yaml into uniform_seed42/43)

Outputs:
- examples/benchmark/degradation/            inputs and per-run PhyloDater outputs
- examples/benchmark/degradation_results.json  per-node, per-seed deviations
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCH = PROJECT_ROOT / "examples" / "benchmark"
DEG = BENCH / "degradation"

NEWICK = (
    "((human:0.03,(chimp:0.01,bonobo:0.01):0.02):0.05,"
    "(gorilla:0.06,(orangutan:0.04,sumatran:0.04):0.02):0.02);"
)
TAXA = ["human", "chimp", "bonobo", "gorilla", "orangutan", "sumatran"]
SEQ_LEN = 10000
SIM_SEED = 42
# Substitution rate (subs/site/My). The human-(chimp,bonobo) path is 0.05
# subs/site, so the simulated HumanChimp age is c. 6 Ma and the root c. 12 Ma.
RATE = 0.05 / 6.0

GAMMA = {
    "HumanChimp": {"alpha": 50.0, "beta": 50.0 / 6.0},
    "GreatApes": {"alpha": 50.0, "beta": 50.0 / 9.0},
}
ROOT_MAX = 15.0


def simulate_alignment() -> str:
    """Simulate a JC69 strict-clock alignment along NEWICK (fixed seed)."""
    rng = np.random.default_rng(SIM_SEED)
    bases = np.array(["A", "C", "G", "T"])

    # Parse the newick into (parent, child, branch) triples via Bio.Phylo
    from io import StringIO

    from Bio import Phylo

    tree = Phylo.read(StringIO(NEWICK), "newick")
    sequences = {}  # clade id -> array of base indices
    root_seq = rng.integers(0, 4, size=SEQ_LEN)
    sequences[id(tree.root)] = root_seq

    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        parent_seq = None
        # find parent sequence by walking up
        path = tree.get_path(clade)
        parent = tree.root if len(path) == 1 else path[-2]
        parent_seq = sequences[id(parent)]
        branch = clade.branch_length
        x = np.exp(-4.0 * branch / 3.0)
        p_same = 0.25 + 0.75 * x
        _p_change_each = (1.0 - x) / 4.0
        # JC69 transition: with prob p_same stay, else uniform over the other 3
        u = rng.random(SEQ_LEN)
        change = u > p_same
        child_seq = parent_seq.copy()
        # draw replacement bases uniformly from the 4 and re-draw if identical
        new_bases = rng.integers(0, 4, size=int(change.sum()))
        same_mask = new_bases == parent_seq[change]
        # resample identical draws once (approximation <1e-18 error per site)
        n_resample = int(same_mask.sum())
        if n_resample:
            redraw = rng.integers(0, 4, size=n_resample)
            still_same = redraw == parent_seq[change][same_mask]
            redraw[still_same] = (
                redraw[still_same] + rng.integers(1, 4, size=int(still_same.sum()))
            ) % 4
            new_bases[same_mask] = redraw
        child_seq[change] = new_bases
        sequences[id(clade)] = child_seq

    lines = []
    for tip in tree.get_terminals():
        name = tip.name
        seq = sequences[id(tip)]
        lines.append(f">{name}")
        lines.append("".join(bases[seq]))
    return "\n".join(lines) + "\n"


GAMMA_YAML = """calibrations:
  - name: "HumanChimp"
    mrca_pair: ["human", "chimp"]
    constraint:
      type: "gamma"
      alpha: 50.0
      beta: 8.333333333333334
  - name: "GreatApes"
    mrca_pair: ["human", "gorilla"]
    constraint:
      type: "gamma"
      alpha: 50.0
      beta: 5.555555555555555
  - name: "Root"
    is_root: true
    constraint:
      type: "maximum"
      max: 15.0
"""


def uniform_yaml() -> str:
    """Uniform degradations: 95% credible intervals of the gammas above,
    computed with the same formula as phylodater.models.constraints."""
    lines = ["calibrations:"]
    for name, p in GAMMA.items():
        scale = 1.0 / p["beta"]
        lo = stats.gamma.ppf(0.025, p["alpha"], scale=scale)
        hi = stats.gamma.ppf(0.975, p["alpha"], scale=scale)
        pair_taxon = "chimp" if name == "HumanChimp" else "gorilla"
        lines += [
            f'  - name: "{name}"',
            f'    mrca_pair: ["human", "{pair_taxon}"]',
            "    constraint:",
            '      type: "uniform"',
            f"      min: {lo:.6f}",
            f"      max: {hi:.6f}",
        ]
    lines += [
        '  - name: "Root"',
        "    is_root: true",
        "    constraint:",
        '      type: "maximum"',
        f"      max: {ROOT_MAX}",
        "",
    ]
    return "\n".join(lines)


def prepare() -> None:
    DEG.mkdir(parents=True, exist_ok=True)
    (DEG / "tree.nwk").write_text(NEWICK + "\n")
    (DEG / "alignment.fasta").write_text(simulate_alignment())
    (DEG / "calibrations_gamma.yaml").write_text(GAMMA_YAML)
    (DEG / "calibrations_uniform.yaml").write_text(uniform_yaml())
    print(f"inputs prepared under {DEG}")
    print("next: run the four MCMCTree analyses as documented in the docstring")


def read_tsv_means(run_dir: Path) -> dict:
    tsv = run_dir / "comparison_table.tsv"
    rows = {}
    with open(tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        mean_col = next(c for c in reader.fieldnames if c.endswith("_mean"))
        for row in reader:
            rows[row["Node"]] = row[mean_col]
    return rows


def analyze() -> None:
    seeds = (42, 43)
    deviations = {}  # node -> {seed: percent}
    ages = {}
    for node_label in (
        "HumanChimp",
        "GreatApes",
        "Root",
        "internal_node_1",
        "internal_node_2",
        "internal_node_3",
    ):
        deviations[node_label] = {}
        ages[node_label] = {}
        for seed in seeds:
            gamma_v = read_tsv_means(DEG / f"gamma_seed{seed}").get(node_label, "NA")
            uni_v = read_tsv_means(DEG / f"uniform_seed{seed}").get(node_label, "NA")
            if gamma_v in ("NA", "", None) or uni_v in ("NA", "", None):
                deviations[node_label][seed] = None
                continue
            g, u = float(gamma_v), float(uni_v)
            ages[node_label][seed] = {"gamma": g, "uniform": u}
            deviations[node_label][seed] = (u - g) / g * 100.0

    flat = [
        abs(v) for node in deviations.values() for v in node.values() if v is not None
    ]
    summary = {
        "experiment": "gamma_to_uniform_degradation_mcmctree",
        "description": (
            "MCMCTree node-age deviations when gamma priors are replaced by "
            "their 95%-CI uniform degradations (the rule implemented in "
            "phylodater.models.constraints); simulated JC69 data, 10 kb, "
            "two seeds, one MCMC chain per run."
        ),
        "rate_subs_site_My": RATE,
        "sim_seed": SIM_SEED,
        "gamma_params": GAMMA,
        "uniform_intervals": {
            name: {
                "lower": float(
                    stats.gamma.ppf(0.025, p["alpha"], scale=1.0 / p["beta"])
                ),
                "upper": float(
                    stats.gamma.ppf(0.975, p["alpha"], scale=1.0 / p["beta"])
                ),
            }
            for name, p in GAMMA.items()
        },
        "node_ages": ages,
        "deviation_percent": deviations,
        "mean_abs_deviation_percent": float(np.mean(flat)) if flat else None,
        "max_abs_deviation_percent": float(np.max(flat)) if flat else None,
    }
    out_json = BENCH / "degradation_results.json"
    out_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(
        f"mean|dev| = {summary['mean_abs_deviation_percent']:.2f}%  "
        f"max|dev| = {summary['max_abs_deviation_percent']:.2f}%"
    )
    print(f"written: {out_json}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("stage", choices=["prepare", "analyze"])
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare()
    else:
        analyze()


if __name__ == "__main__":
    main()
