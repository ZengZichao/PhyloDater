#!/usr/bin/env python3
"""PhyloDater reference-benchmark harness.

Runs every supported dating method over the accuracy benchmark dataset in
``test-data/benchmark/`` (a strict-clock simulation with **known** divergence
times), then scores each method on:

* completion        - the CLI exits 0 and reports ``status: completed``
* tip completeness  - all input taxa are present in the dated tree
* topology          - the set of internal clades of the dated tree equals the
                      set of internal clades of the input tree
* chronology        - node ages are non-negative and every parent is older than
                      each of its children
* calibration       - every calibration is honoured within a stated tolerance
                      (the value the method was *given*, not the truth)
* accuracy          - mean / root-mean-square / maximum absolute error against
                      the simulated truth, in Ma and as a relative error
* reproducibility   - optional ``--repeat`` runs must give the same estimates
                      for deterministic methods

Outputs (under ``--output``):

``runs/<method>/``        full PhyloDater output directory for each method
``node_ages.tsv``         one row per method x node (truth, estimate, errors)
``scorecard.tsv``         one row per method (all checks + accuracy summary)
``benchmark_results.json`` machine-readable version of everything above

Usage:
    python scripts/run_reference_benchmark.py \
        --dataset test-data/benchmark --output test-data/benchmark/results \
        --methods lsd2 pathd8 treepl r8s wlogdate mdcat mcmctree
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import dendropy

REPO_ROOT = Path(__file__).resolve().parent.parent


def _neutral_output_dir(run_dir: Path, out: Path) -> str:
    """Return a machine-neutral path for the committed benchmark artifacts.

    ``benchmark_results.json`` is checked into the repository, so an absolute
    ``run_dir`` would leak the operator's home/checkout path (non-reproducible
    across machines and a privacy leak). Prefer a path relative to the results
    directory, then to the repo root, and finally just the run directory name.
    """
    for base in (out, REPO_ROOT):
        try:
            return str(Path(run_dir).resolve().relative_to(Path(base).resolve()))
        except ValueError:
            continue
    return Path(run_dir).name


# Methods supported by the CLI (``phylodater dating --method``).
ALL_METHODS = ["mcmctree", "lsd2", "r8s", "treepl", "pathd8", "wlogdate", "mdcat"]

# Per-method CLI arguments used for the shipped benchmark. They are deliberately
# modest so that the whole matrix finishes in minutes on a laptop, while still
# exercising the real engine (no mocks, no dry-runs). Anything reduced from the
# tool default is documented in docs/TEST_RESULTS.md.
DEFAULT_METHOD_ARGS: Dict[str, List[str]] = {
    "mcmctree": [
        "--method-args",
        "clock=2,num_runs=2,burnin=2000,nsample=4000,sampfreq=10",
    ],
    "mdcat": ["--method-args", "ncat=10,nrep=20,max_iter=50"],
}

# Tolerances used by the scorecard.
FIXED_TOLERANCE_MA = 1.0  # point calibration: |estimate - age| <= 1 Ma
RANGE_TOLERANCE_FRACTION = 0.05  # range calibration: 5 % of the range width
ACCURACY_MAE_LIMIT_MA = 10.0  # per-method MAE limit against the truth
ACCURACY_MAX_REL_LIMIT = 0.5  # per-method worst-node relative-error limit


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #


@dataclass
class NodeAge:
    node: str
    mrca_pair: Tuple[str, str]
    true_age: float
    estimate: Optional[float]
    ci_lower: Optional[float] = None
    ci_upper: Optional[float] = None


@dataclass
class MethodRun:
    method: str
    output_dir: Path
    exit_code: int = -1
    wall_seconds: float = 0.0
    status: str = "not-run"
    software_version: str = ""
    node_count: int = 0
    warnings: List[str] = field(default_factory=list)
    checks: Dict[str, object] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    node_ages: List[NodeAge] = field(default_factory=list)
    disclosed_limitations: List[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.errors and all(
            bool(v) for k, v in self.checks.items() if k != "accuracy"
        )


# --------------------------------------------------------------------------- #
# Dataset helpers
# --------------------------------------------------------------------------- #


def load_truth(dataset: Path) -> Dict[str, Dict]:
    """Read ``truth.yaml`` into ``{node: {mrca_pair, true_age_ma}}``."""
    import yaml

    data = yaml.safe_load((dataset / "truth.yaml").read_text(encoding="utf-8"))
    nodes = {}
    for name, entry in (data.get("nodes") or {}).items():
        pair = tuple(entry["mrca_pair"])
        nodes[name] = {"mrca_pair": pair, "true_age_ma": float(entry["true_age_ma"])}
    return nodes


def load_dataset_meta(dataset: Path) -> Dict:
    import yaml

    data = yaml.safe_load((dataset / "truth.yaml").read_text(encoding="utf-8"))
    return data.get("dataset") or {}


def tip_name(taxon) -> str:
    """Taxon label that works on DendroPy 4.x (``name``) and 5.x (``label``)."""
    return (getattr(taxon, "label", None) or getattr(taxon, "name", None) or "").strip()


def leaf_names(tree: dendropy.Tree) -> List[str]:
    return sorted(tip_name(taxon) for taxon in tree.taxon_namespace)


def input_clades(tree_path: Path) -> Tuple[Set[Tuple[str, ...]], List[str]]:
    """Internal clades (as sorted leaf tuples) and the tip list of the input tree."""
    tree = dendropy.Tree.get(
        data=tree_path.read_text(encoding="utf-8"),
        schema="newick",
        preserve_underscores=True,
        rooting="force-rooted",
    )
    return clades_of(tree), leaf_names(tree)


def clades_of(tree: dendropy.Tree) -> Set[Tuple[str, ...]]:
    clades: Set[Tuple[str, ...]] = set()

    def leaves_of(node) -> List[str]:
        if node.is_leaf():
            return [tip_name(node.taxon)]
        out: List[str] = []
        for child in node.child_nodes():
            out.extend(leaves_of(child))
        return out

    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        clades.add(tuple(sorted(leaves_of(node))))
    return clades


def read_dated_tree(run_dir: Path, method: str) -> Optional[dendropy.Tree]:
    """Load the standardised dated tree written by PhyloDater for ``method``."""
    candidates = sorted(run_dir.glob(f"{method}_dated_tree.*"))
    for path in candidates:
        text = path.read_text(encoding="utf-8")
        if path.suffix.lower() in (".nexus", ".nex", ".nhx"):
            text = _strip_nexus_annotations(text)
        text = _strip_internal_labels(text)
        try:
            return dendropy.Tree.get(
                data=text,
                schema="newick",
                preserve_underscores=True,
                rooting="force-rooted",
            )
        except Exception:
            continue
    return None


_ANNOTATION_RE = re.compile(r"\[[^\]]*\]")
_DATE_ATTR_RE = re.compile(r"\[[^\]]*\]")


def _strip_nexus_annotations(text: str) -> str:
    """Reduce a NEXUS/NHX time tree to a bare labelled Newick string."""
    if text.lstrip().upper().startswith("#NEXUS"):
        match = re.search(r"(?im)^\s*TREE\s+\S+\s*=\s*(.+)$", text)
        if match:
            text = match.group(1)
    text = text.strip()
    if text.upper().startswith("[&R]"):
        text = text[4:].strip()
    return text


def _strip_internal_labels(text: str) -> str:
    """Remove ``)Label`` and ``[...]`` annotations, keeping tip names and lengths."""

    def repl(match: re.Match) -> str:
        return match.group(1)

    cleaned = _ANNOTATION_RE.sub("", text)
    cleaned = re.sub(r"\)\s*[A-Za-z0-9_.\-]+(?=:|,|\)|;)", r")", cleaned)
    return cleaned


def node_ages_from_tree(tree: dendropy.Tree) -> Dict[Tuple[str, ...], float]:
    """``{clade: age}`` for every internal node of a chronogram.

    A dated tree writes time as branch *differences*: the distance from the root
    to a tip equals the root age when all tips are contemporaneous. Node ages are
    therefore measured **back from the present plane**, i.e.
    ``age(node) = max_root_to_tip_distance - distance_from_root(node)``; taking
    the raw distance from the root would return "time elapsed since the root"
    (0 for the root, largest for the youngest nodes), i.e. exactly inverted ages.
    Using the maximum (rather than assuming ultrametricity) keeps the estimate
    well defined for slightly non-ultrametric outputs.
    """
    heights: Dict[Tuple[str, ...], float] = {}
    tip_distances: List[float] = []

    def walk(node, accumulated: float) -> List[str]:
        height = accumulated + (node.edge.length or 0.0)
        if node.is_leaf():
            tip_distances.append(height)
            return [tip_name(node.taxon)]
        leaves: List[str] = []
        for child in node.child_nodes():
            leaves.extend(walk(child, height))
        heights[tuple(sorted(leaves))] = height
        return leaves

    seed = tree.seed_node
    for child in seed.child_nodes():
        walk(child, 0.0)
    heights[tuple(leaf_names(tree))] = 0.0

    origin = max(tip_distances) if tip_distances else 0.0
    return {clade: origin - height for clade, height in heights.items()}


def infer_unit_factor(root_age: float, expected_root: float) -> float:
    """Choose the branch-length unit (Ma / Ga / ka->Ma) closest to the calibration.

    MCMCTree writes its chronogram in Ga while every other backend writes Ma;
    rather than hard-coding that per method, the scale is inferred from the root
    calibration, which every dataset in this harness provides.
    """
    if root_age <= 0:
        return 1.0
    factors = [1.0, 1000.0, 1e-3, 1e6]
    best = min(factors, key=lambda f: abs(math_log_ratio(root_age * f, expected_root)))
    return best


def math_log_ratio(value: float, reference: float) -> float:
    if value <= 0 or reference <= 0:
        return float("inf")
    import math

    return abs(math.log(value / reference))


# --------------------------------------------------------------------------- #
# Calibration checks
# --------------------------------------------------------------------------- #


def load_calibrations(dataset: Path) -> List[Dict]:
    import yaml

    data = yaml.safe_load((dataset / "calibrations.yaml").read_text(encoding="utf-8"))
    return data.get("calibrations") or []


def calibration_target(calibration: Dict) -> Tuple[str, float, float, float]:
    """Return (kind, value, low, high) for a calibration entry.

    kind is ``fixed`` (point) or ``range`` (two-sided) -- the benchmark dataset
    only uses those two so that every backend can express the constraint without
    semantic degradation.
    """
    constraint = calibration["constraint"]
    ctype = constraint.get("type")
    if ctype == "fixed":
        age = float(constraint["age"])
        return "fixed", age, age, age
    if ctype == "uniform":
        low, high = float(constraint["min"]), float(constraint["max"])
        return "range", (low + high) / 2.0, low, high
    raise ValueError(f"unsupported calibration type in benchmark dataset: {ctype}")


def mrca_clade(
    pair: Tuple[str, ...], heights: Dict[Tuple[str, ...], float]
) -> Optional[Tuple[str, ...]]:
    """Smallest recorded clade containing both leaves of ``pair``."""
    contained = [clade for clade in heights if set(pair) <= set(clade)]
    if not contained:
        return None
    return min(contained, key=len)


# --------------------------------------------------------------------------- #
# Running PhyloDater
# --------------------------------------------------------------------------- #


def run_method(
    method: str,
    dataset: Path,
    tree_file: str,
    run_dir: Path,
    python: str,
    timeout: int,
    extra_args: Sequence[str],
    seed: int,
) -> MethodRun:
    record = MethodRun(method=method, output_dir=run_dir)
    if run_dir.exists():
        shutil_rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        python,
        "-m",
        "phylodater",
        "dating",
        "-t",
        str(dataset / tree_file),
        "-s",
        str(dataset / "alignment.fasta"),
        "-c",
        str(dataset / "calibrations.yaml"),
        "-o",
        str(run_dir),
        "--method",
        method,
        "--seed",
        str(seed),
        *extra_args,
    ]
    started = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        record.exit_code = proc.returncode
    except subprocess.TimeoutExpired:
        record.exit_code = -9
        record.errors.append(f"timed out after {timeout}s")
        return record
    record.wall_seconds = round(time.monotonic() - started, 2)

    metadata_path = run_dir / "runtime_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        method_meta = (metadata.get("methods") or {}).get(method, {})
        record.status = method_meta.get("status", "unknown")
        record.software_version = (method_meta.get("software_versions") or {}).get(
            method, ""
        )
        record.node_count = method_meta.get("node_count", 0) or 0
        record.warnings = list(method_meta.get("warnings") or [])
    if record.exit_code != 0:
        record.errors.append(f"CLI exit code {record.exit_code}")
    if record.status != "completed":
        record.errors.append(f"pipeline status '{record.status}'")
    return record


def shutil_rmtree(path: Path) -> None:
    import shutil

    shutil.rmtree(path, ignore_errors=True)


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #


def score_method(
    record: MethodRun,
    dataset: Path,
    truth: Dict[str, Dict],
    input_clade_set: Set[Tuple[str, ...]],
    tips: List[str],
) -> None:
    checks: Dict[str, object] = {}
    tree = read_dated_tree(record.output_dir, record.method)
    if tree is None:
        record.errors.append("no dated tree produced")
        record.checks = checks
        return

    dated_tips = leaf_names(tree)
    checks["tips_complete"] = dated_tips == tips
    if not checks["tips_complete"]:
        record.errors.append(
            f"tip mismatch: missing={set(tips) - set(dated_tips)} "
            f"extra={set(dated_tips) - set(tips)}"
        )

    heights = node_ages_from_tree(tree)
    root_clade = tuple(sorted(dated_tips))
    root_calibration = next(
        (c for c in load_calibrations(dataset) if c.get("is_root")), None
    )
    if root_calibration is not None:
        _, mid, _, _ = calibration_target(root_calibration)
        factor = infer_unit_factor(heights.get(root_clade, 0.0), mid)
    else:
        factor = 1.0
    heights = {clade: age * factor for clade, age in heights.items()}

    # Topology: the dated tree must contain exactly the input clades.
    dated_clades = {clade for clade in heights if len(clade) > 1}
    checks["topology_identical"] = dated_clades == input_clade_set
    if not checks["topology_identical"]:
        record.errors.append(
            "topology changed: "
            f"lost={sorted(input_clade_set - dated_clades)} "
            f"gained={sorted(dated_clades - input_clade_set)}"
        )

    # Chronology: non-negative, parent older than children.
    ages = [age for clade, age in heights.items() if len(clade) > 1]
    checks["ages_nonnegative"] = bool(ages) and all(age >= 0 for age in ages)
    monotone = True
    for clade, age in heights.items():
        for child_clade, child_age in heights.items():
            if set(child_clade) < set(clade) and child_age > age + 1e-6:
                monotone = False
    checks["chronology_monotone"] = monotone
    if not checks["ages_nonnegative"]:
        record.errors.append("negative node age in output")
    if not monotone:
        record.errors.append("a child node is older than its parent")

    # Calibration honouring.
    #
    # The gate is PhyloDater's own contract, not the engines' internals: a requested
    # calibration must either be satisfied by the result **or** be explicitly
    # disclosed by the adapter as not satisfied. Engines differ in what a range
    # means (pyr8s treats CONSTRAIN as a relaxed barrier penalty, so NPRS can land
    # outside it), and silently returning an out-of-range age is the actual defect.
    calib_ok = True
    disclosed: List[str] = []
    undisclosed: List[str] = []
    warning_text = " ".join(str(w) for w in record.warnings)
    for calibration in load_calibrations(dataset):
        kind, mid, low, high = calibration_target(calibration)
        pair = (
            root_clade
            if calibration.get("is_root")
            else tuple(sorted(calibration["mrca_pair"]))
        )
        clade = mrca_clade(pair, heights)
        if clade is None:
            record.errors.append(f"calibrated node not found: {calibration['name']}")
            calib_ok = False
            continue
        estimate = heights[clade]
        if kind == "fixed":
            tol = FIXED_TOLERANCE_MA
        else:
            tol = max(1.0, (high - low) * RANGE_TOLERANCE_FRACTION)
        inside = (low - tol) <= estimate <= (high + tol)
        if not inside:
            marker = f"'{calibration['name']}' was NOT satisfied"
            alternative = f"'{calibration['name']}' was NOT written"
            if marker in warning_text or alternative in warning_text:
                disclosed.append(
                    f"{calibration['name']} ({estimate:.2f} Ma outside [{low}, {high}])"
                )
            else:
                undisclosed.append(
                    f"{calibration['name']}: estimate {estimate:.2f} Ma outside "
                    f"[{low}, {high}] Ma (+-{tol:.2f}) and NOT disclosed"
                )
                calib_ok = False
    checks["calibrations_honoured"] = calib_ok
    record.checks = checks
    record.disclosed_limitations = list(disclosed)

    # Accuracy against the simulated truth, matched by clade.
    rows: List[NodeAge] = []
    for name, entry in truth.items():
        clade = mrca_clade(tuple(entry["mrca_pair"]), heights)
        estimate = heights.get(clade) if clade else None
        rows.append(
            NodeAge(
                node=name,
                mrca_pair=entry["mrca_pair"],
                true_age=entry["true_age_ma"],
                estimate=None if estimate is None else round(estimate, 4),
            )
        )
    record.node_ages = rows

    scored = [r for r in rows if r.estimate is not None]
    checks["all_nodes_recovered"] = len(scored) == len(rows)
    if not scored:
        record.errors.append("no truth node could be matched in the output")
        record.checks = checks
        return

    errors = [abs(r.estimate - r.true_age) for r in scored]
    rel = [abs(r.estimate - r.true_age) / r.true_age for r in scored if r.true_age > 0]
    mae = sum(errors) / len(errors)
    rmse = (sum(e * e for e in errors) / len(errors)) ** 0.5
    max_rel = max(rel) if rel else float("inf")
    record.checks = checks
    record.checks["accuracy"] = {
        "nodes_scored": len(scored),
        "mae_ma": round(mae, 3),
        "rmse_ma": round(rmse, 3),
        "max_abs_error_ma": round(max(errors), 3),
        "max_relative_error": round(max_rel, 4),
        "median_relative_error": round(sorted(rel)[len(rel) // 2], 4) if rel else None,
    }
    if not checks["all_nodes_recovered"]:
        record.errors.append(
            f"{len(rows) - len(scored)} calibrated clades missing from the output"
        )
    if mae > ACCURACY_MAE_LIMIT_MA:
        record.errors.append(
            f"MAE {mae:.2f} Ma exceeds the {ACCURACY_MAE_LIMIT_MA} Ma limit"
        )
    if max_rel > ACCURACY_MAX_REL_LIMIT:
        record.errors.append(
            f"worst-node relative error {max_rel:.2f} exceeds {ACCURACY_MAX_REL_LIMIT}"
        )


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #


def write_reports(
    runs: List[MethodRun], out: Path, dataset_meta: Dict, tree_file: str
) -> None:
    scorecard = out / "scorecard.tsv"
    lines = [
        "\t".join(
            [
                "method",
                "status",
                "exit_code",
                "wall_seconds",
                "software_version",
                "nodes_recovered",
                "mae_ma",
                "rmse_ma",
                "max_abs_error_ma",
                "max_relative_error",
                "checks",
                "warnings",
                "disclosed_limitations",
                "errors",
            ]
        )
    ]
    detail = [
        "method\tnode\tmrca_pair\ttrue_age_ma\testimate_ma\tabs_error_ma\trel_error"
    ]
    payload = {
        "dataset": dataset_meta,
        "tree_variant": tree_file,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runs": [],
    }
    for run in runs:
        accuracy = run.checks.get("accuracy") or {}
        check_text = ";".join(
            f"{k}={'PASS' if v else 'FAIL'}"
            for k, v in run.checks.items()
            if k != "accuracy"
        )
        lines.append(
            "\t".join(
                str(x)
                for x in [
                    run.method,
                    run.status,
                    run.exit_code,
                    run.wall_seconds,
                    run.software_version or "n/a",
                    accuracy.get("nodes_scored", 0),
                    accuracy.get("mae_ma", "NA"),
                    accuracy.get("rmse_ma", "NA"),
                    accuracy.get("max_abs_error_ma", "NA"),
                    accuracy.get("max_relative_error", "NA"),
                    check_text,
                    len(run.warnings),
                    "; ".join(run.disclosed_limitations) or "none",
                    " | ".join(run.errors) if run.errors else "none",
                ]
            )
        )
        for row in run.node_ages:
            abs_err = (
                "NA"
                if row.estimate is None
                else f"{abs(row.estimate - row.true_age):.4f}"
            )
            rel_err = (
                "NA"
                if row.estimate is None or row.true_age == 0
                else f"{abs(row.estimate - row.true_age) / row.true_age:.4f}"
            )
            detail.append(
                "\t".join(
                    [
                        run.method,
                        row.node,
                        f"{row.mrca_pair[0]}|{row.mrca_pair[1]}",
                        f"{row.true_age:.3f}",
                        "NA" if row.estimate is None else f"{row.estimate:.4f}",
                        abs_err,
                        rel_err,
                    ]
                )
            )
        payload["runs"].append(
            {
                "method": run.method,
                "status": run.status,
                "exit_code": run.exit_code,
                "wall_seconds": run.wall_seconds,
                "software_version": run.software_version,
                "checks": {k: v for k, v in run.checks.items() if k != "accuracy"},
                "accuracy": run.checks.get("accuracy"),
                "warnings": run.warnings,
                "disclosed_limitations": run.disclosed_limitations,
                "errors": run.errors,
                "node_ages": [
                    {
                        "node": r.node,
                        "mrca_pair": list(r.mrca_pair),
                        "true_age_ma": r.true_age,
                        "estimate_ma": r.estimate,
                    }
                    for r in run.node_ages
                ],
                "output_dir": _neutral_output_dir(run.output_dir, out),
            }
        )
    scorecard.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out / "node_ages.tsv").write_text("\n".join(detail) + "\n", encoding="utf-8")
    (out / "benchmark_results.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def print_summary(runs: List[MethodRun]) -> int:
    header = (
        f"{'method':<10} {'status':<10} {'wall':>8} {'MAE(Ma)':>9} "
        f"{'maxRel':>8} {'nodes':>6}  checks"
    )
    print(header)
    print("-" * len(header))
    failures = 0
    for run in runs:
        accuracy = run.checks.get("accuracy") or {}
        check_names = [k for k, v in run.checks.items() if k != "accuracy"]
        failed = [k for k in check_names if not run.checks[k]]
        verdict = "PASS" if run.passed else "FAIL"
        if not run.passed:
            failures += 1
        print(
            f"{run.method:<10} {verdict:<10} {run.wall_seconds:>8.2f} "
            f"{accuracy.get('mae_ma', float('nan')):>9} "
            f"{accuracy.get('max_relative_error', 'NA'):>8} "
            f"{accuracy.get('nodes_scored', 0):>6}  "
            f"{'all passed' if not failed else 'FAILED: ' + ','.join(failed)}"
        )
        for err in run.errors:
            print(f"{'':<10} ! {err}")
        for note in run.disclosed_limitations:
            print(f"{'':<10} ~ disclosed engine limitation: {note}")
    return failures


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--dataset", type=Path, default=REPO_ROOT / "test-data" / "benchmark"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="results directory (default: <dataset>/results)",
    )
    parser.add_argument(
        "--methods", nargs="*", default=ALL_METHODS, choices=ALL_METHODS
    )
    parser.add_argument(
        "--tree-file",
        default="tree.nwk",
        help="tree inside --dataset: 'tree.nwk' (expected clock lengths) or "
        "'tree_ml.nwk' (ML-estimated lengths, the realistic variant)",
    )
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--method-args",
        nargs="*",
        default=[],
        metavar="METHOD=key=value,key=value",
        help="override the benchmark parameters for one method, e.g. "
        "--method-args mdcat:ncat=50,nrep=100,max_iter=100 "
        "mcmctree:clock=1 (use this to quantify sensitivity to the reduced defaults)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    method_args = dict(DEFAULT_METHOD_ARGS)
    for spec in args.method_args:
        if ":" not in spec:
            parser.error(f"--method-args expects METHOD=key=value,... got {spec!r}")
        method, pairs = spec.split(":", 1)
        if method not in ALL_METHODS:
            parser.error(f"unknown method {method!r} in --method-args")
        method_args[method] = ["--method-args", pairs]

    dataset = args.dataset.resolve()
    out = (args.output or (dataset / "results")).resolve()
    out.mkdir(parents=True, exist_ok=True)

    truth = load_truth(dataset)
    dataset_meta = load_dataset_meta(dataset)
    tree_path = dataset / args.tree_file
    if not tree_path.exists():
        print(
            f"error: {tree_path} not found (run "
            f"'python test-data/benchmark/generate_benchmark_dataset.py --force "
            f"--include-ml-tree' to create tree_ml.nwk)",
            file=sys.stderr,
        )
        return 2
    input_clade_set, tips = input_clades(tree_path)

    runs_root = out / "runs"
    runs_root.mkdir(parents=True, exist_ok=True)

    runs: List[MethodRun] = []
    for method in args.methods:
        print(f"[benchmark] running {method} ...", flush=True)
        record = run_method(
            method=method,
            dataset=dataset,
            tree_file=args.tree_file,
            run_dir=runs_root / method,
            python=args.python,
            timeout=args.timeout,
            extra_args=method_args.get(method, []),
            seed=args.seed,
        )
        if record.exit_code == 0 and record.status == "completed":
            score_method(record, dataset, truth, input_clade_set, tips)
        runs.append(record)

    write_reports(runs, out, dataset_meta, args.tree_file)
    failures = print_summary(runs)
    print(f"\n[benchmark] reports written to {out}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
