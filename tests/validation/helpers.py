"""Shared helpers for the comprehensive test suite."""

import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

from phylodater.infrastructure.safe_io import safe_writer

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

SUITE_ROOT = Path(__file__).resolve().parent
# The suite lives at ``<project root>/tests/validation``.
PROJECT_ROOT = SUITE_ROOT.parent.parent
FIXTURES_DIR = SUITE_ROOT / "fixtures"


# --------------------------------------------------------------------------- #
# External binary availability
# --------------------------------------------------------------------------- #


def external_available(name: str) -> bool:
    """Check whether an external dating tool is available.

    The check maps the PhyloDater method name to the concrete binary / Python
    package that backs it:

    - ``lsd2``     -> IQ-TREE binary (``iqtree2`` / ``iqtree``)
    - ``pathd8``   -> PATHd8 executable (``PATHd8`` / ``pathd8``)
    - ``r8s``      -> ``pyr8s`` Python package (native NPRS API) or the ``r8s`` binary
    - ``wlogdate`` -> ``logdate`` Python package or the ``launch_wLogDate.py`` script
    - ``mdcat``    -> ``emd`` Python package or the ``mdcat``/``md_cat.py`` executable
    """
    if name == "wlogdate":
        try:
            import logdate  # noqa: F401

            return True
        except ImportError:
            return shutil.which("launch_wLogDate.py") is not None
    if name == "mdcat":
        try:
            from emd.emd_normal_lib import MDCat  # noqa: F401

            return True
        except ImportError:
            return shutil.which("mdcat") is not None
    if name == "r8s":
        try:
            from pyr8s.core import RateAnalysis  # noqa: F401

            return True
        except ImportError:
            return shutil.which("r8s") is not None
    if name == "lsd2":
        return shutil.which("iqtree2") is not None or shutil.which("iqtree") is not None
    if name == "pathd8":
        return shutil.which("PATHd8") is not None or shutil.which("pathd8") is not None
    return shutil.which(name) is not None


# --------------------------------------------------------------------------- #
# CLI helpers
# --------------------------------------------------------------------------- #


def run_cli(argv: List[str], cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    """Run the phylodater CLI in a subprocess and return the CompletedProcess."""
    env = os.environ.copy()
    env["PYTHONPATH"] = str(PROJECT_ROOT)
    cmd = ["python", "-m", "phylodater"] + argv
    return subprocess.run(
        cmd,
        cwd=str(cwd or PROJECT_ROOT),
        capture_output=True,
        text=True,
        env=env,
        timeout=300,
    )


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #


def write_text(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def make_small_tree() -> str:
    return "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"


def make_medium_tree() -> str:
    return (
        "((((A:0.05,B:0.05):0.1,C:0.15):0.1,"
        "((D:0.05,E:0.05):0.1,F:0.15):0.1):0.2,"
        "((G:0.05,H:0.05):0.1,(I:0.05,J:0.05):0.1):0.3);"
    )


def make_unrooted_tree() -> str:
    return "(A:0.1,B:0.2,C:0.3,D:0.4);"


def make_negative_branch_tree() -> str:
    return "((A:0.1,B:-0.2):0.3,(C:0.4,D:0.5):0.6);"


def make_alignment(tip_names: List[str], length: int = 60, char: str = "A") -> str:
    lines = []
    for name in tip_names:
        lines.append(f">{name}")
        lines.append(char * length)
    return "\n".join(lines) + "\n"


def make_calibration(
    name: str,
    constraint_type: str,
    tip_names: Optional[List[str]] = None,
    **kwargs,
) -> dict:
    entry = {"name": name, "constraint": {"type": constraint_type, **kwargs}}
    if tip_names:
        entry["node"] = {"type": "mrca", "taxa": tip_names}
    return entry


def make_calibration_file(
    path: Path,
    entries: List[dict],
) -> Path:
    import yaml

    path.parent.mkdir(parents=True, exist_ok=True)
    with safe_writer(path, encoding="utf-8") as fh:
        yaml.dump(
            {"calibrations": entries}, fh, default_flow_style=False, sort_keys=False
        )
    return path


# --------------------------------------------------------------------------- #
# Assertion helpers
# --------------------------------------------------------------------------- #


def assert_file_contains(path: Path, substring: str) -> None:
    text = path.read_text(encoding="utf-8")
    assert substring in text, f"Expected {path} to contain {substring!r}"


def assert_file_exists(path: Path) -> None:
    assert path.exists(), f"Expected file to exist: {path}"
