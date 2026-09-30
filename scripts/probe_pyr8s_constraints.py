"""Characterise pyr8s' handling of CONSTRAIN (range) vs FIXAGE (point) calibrations.

Runs the same ML-branch-length benchmark tree through pyr8s in three configurations
and prints the ages of the calibrated nodes, to establish whether a ``constrain``
range is enforced, advisory, or ignored. This is the evidence behind the
"pyr8s ranges are advisory" disclosure in ``phylodater/adapters/r8s_pyr8s_method.py``
and docs/TEST_RESULTS.md.

Usage:  python scripts/probe_pyr8s_constraints.py     (requires the pyr8s package)
"""

import re
import tempfile
from pathlib import Path

from pyr8s.parse import from_file_nexus

REPO_ROOT = Path(__file__).resolve().parent.parent
BASE = REPO_ROOT / "test-data" / "benchmark"

TREE = re.sub(r"^\s*(?:\[[^\]]*\]\s*)+", "", (BASE / "tree_ml.nwk").read_text().strip())

HEAD = """#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=13;
    TAXLABELS
        Dog Rabbit Mouse Rat Lemur Marmoset Macaque Gibbon Orangutan Gorilla Human Chimp Bonobo
    ;
END;

BEGIN TREES;
    TREE tree1 = {tree}
END;
"""

RATES = """
BEGIN RATES;
    BLFORMAT LENGTHS=PERSITE NSITES=3000;
    MRCA Root Dog Rabbit;
{body}
    DIVTIME METHOD={method} ALGORITHM=POWELL;
    DESCRIBE PLOT=CHRONOGRAM;
END;
"""

ALL = """    constrain taxon=Root min_age=100;
    constrain taxon=Root max_age=110;
    MRCA Homininae Human Gorilla;
    fixage taxon=Homininae age=9;
    MRCA Primates Lemur Human;
    constrain taxon=Primates min_age=68;
    constrain taxon=Primates max_age=82;
    MRCA Glires Rabbit Mouse;
    constrain taxon=Glires min_age=60;
    constrain taxon=Glires max_age=80;"""

FIXAGE_ONLY = """    MRCA Homininae Human Gorilla;
    MRCA Primates Lemur Human;
    MRCA Glires Rabbit Mouse;
    fixage taxon=Primates age=75;
    fixage taxon=Glires age=70;
    fixage taxon=Homininae age=9;"""

SCENARIOS = [
    ("NPRS + constrain ranges", ALL, "NPRS"),
    ("NPRS + fixage only", FIXAGE_ONLY, "NPRS"),
    ("LNL + constrain ranges", ALL, "LNL"),
]

BOUNDS = {
    "Root": (100, 110),
    "Primates": (68, 82),
    "Glires": (60, 80),
    "Homininae": (9, 9),
}

for title, body, method in SCENARIOS:
    print("=" * 72)
    print("scenario:", title)
    nexus = HEAD.format(tree=TREE) + RATES.format(body=body, method=method)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "in.nex"
        path.write_text(nexus, encoding="utf-8")
        try:
            analysis = from_file_nexus(str(path), run=True)
            chronogram = analysis.results.chronogram
        except Exception as exc:
            print("  FAILED:", type(exc).__name__, exc)
            continue

    ages = {}

    def leaves(node):
        if node.is_leaf():
            return [node.taxon.label]
        out = []
        for child in node.child_nodes():
            out += leaves(child)
        return out

    tip_max = [0.0]

    def walk(node, acc):
        h = acc + (node.edge.length or 0.0)
        if node.is_leaf():
            tip_max[0] = max(tip_max[0], h)
        else:
            ages[tuple(sorted(leaves(node)))] = h
        for child in node.child_nodes():
            walk(child, h)

    walk(chronogram.seed_node, 0.0)
    computed = {clade: tip_max[0] - h for clade, h in ages.items()}

    def show(name, wanted_pair):
        best = None
        for clade, age in computed.items():
            if set(wanted_pair) <= set(clade):
                if best is None or len(clade) < len(best[0]):
                    best = (clade, age)
        if best is None:
            print(f"  {name:11s} not resolved")
            return
        low, high = BOUNDS[name]
        value = best[1]
        status = "OK" if low - 1 <= value <= high + 1 else "VIOLATED"
        print(f"  {name:11s} = {value:8.3f}   requested [{low}, {high}]   {status}")

    show("Root", ("Dog", "Human"))
    show("Primates", ("Lemur", "Human"))
    show("Glires", ("Rabbit", "Mouse"))
    show("Homininae", ("Human", "Gorilla"))
