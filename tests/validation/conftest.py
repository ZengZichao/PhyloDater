"""Pytest shared configuration and fixtures for the comprehensive test suite."""

import os
import shutil
import sys
from pathlib import Path

import pytest

# --------------------------------------------------------------------------- #
# Path setup
# --------------------------------------------------------------------------- #
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# --------------------------------------------------------------------------- #
# Directories
# --------------------------------------------------------------------------- #
SUITE_ROOT = Path(__file__).resolve().parent
FIXTURES_DIR = SUITE_ROOT / "fixtures"
TMP_ROOT = SUITE_ROOT / "_tmp"


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line("markers", "slow: tests that run real external software")
    config.addinivalue_line(
        "markers", "security: tests for path safety, injection, and input validation"
    )
    config.addinivalue_line("markers", "requires_mcmctree: requires MCMCTree/PAML")
    config.addinivalue_line("markers", "requires_treepl: requires treePL")
    config.addinivalue_line("markers", "requires_pathd8: requires PATHd8")
    config.addinivalue_line("markers", "requires_lsd2: requires LSD2")
    config.addinivalue_line("markers", "requires_r8s: requires r8s")
    config.addinivalue_line("markers", "requires_wlogdate: requires wLogDate")
    config.addinivalue_line("markers", "requires_mdcat: requires MDCAT")


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #


@pytest.fixture
def tmp_test_dir(request):
    """Provide an isolated per-test temporary directory under the test suite.

    The directory is created before the test and removed afterwards unless the
    test fails and the environment variable PHYLODATER_KEEP_TMP is set.
    """
    test_name = request.node.name.replace("/", "_").replace("\\", "_")
    tmp_dir = (
        TMP_ROOT
        / f"{request.node.nodeid.split('::')[0].replace('/', '_').replace('.', '_')}_{test_name}_{os.getpid()}"
    )
    tmp_dir.mkdir(parents=True, exist_ok=True)
    # Restrict permissions like the application does
    try:
        os.chmod(str(tmp_dir), 0o700)
    except OSError:
        pass
    yield tmp_dir
    keep = os.environ.get("PHYLODATER_KEEP_TMP", "0") == "1"
    if not keep:
        shutil.rmtree(tmp_dir, ignore_errors=True)


@pytest.fixture(scope="session")
def fixtures_dir():
    """Return the path to the static fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture(scope="session")
def small_tree_path(fixtures_dir):
    return fixtures_dir / "trees" / "small_sample.nwk"


@pytest.fixture(scope="session")
def taxonomy_tree_path(fixtures_dir):
    return fixtures_dir / "trees" / "taxonomy_tree.nwk"


@pytest.fixture(scope="session")
def taxonomy_tree_alignment_path(fixtures_dir):
    return fixtures_dir / "alignments" / "taxonomy_tree.fasta"


@pytest.fixture(scope="session")
def taxonomy_tree_calibration_path(fixtures_dir):
    return fixtures_dir / "calibrations" / "taxonomy_tree.yaml"


@pytest.fixture(scope="session")
def medium_tree_path(fixtures_dir):
    return fixtures_dir / "trees" / "medium_sample.nwk"


@pytest.fixture(scope="session")
def small_alignment_path(fixtures_dir):
    return fixtures_dir / "alignments" / "small_sample.fasta"


@pytest.fixture(scope="session")
def small_tree_alignment_path(fixtures_dir):
    return fixtures_dir / "alignments" / "small_tree.fasta"


@pytest.fixture(scope="session")
def valid_calibration_path(fixtures_dir):
    return fixtures_dir / "calibrations" / "valid_uniform.yaml"


@pytest.fixture(scope="session")
def small_tree_calibration_path(fixtures_dir):
    return fixtures_dir / "calibrations" / "small_tree_uniform.yaml"


@pytest.fixture(scope="session")
def fixed_calibration_path(fixtures_dir):
    return fixtures_dir / "calibrations" / "valid_fixed.yaml"


@pytest.fixture(scope="session")
def benchmark_tree_path(fixtures_dir):
    """7-taxon primate tree used by the all-dating-tools benchmark tests."""
    return fixtures_dir / "benchmark" / "tree.nwk"


@pytest.fixture(scope="session")
def benchmark_alignment_path(fixtures_dir):
    return fixtures_dir / "benchmark" / "alignment.fasta"


@pytest.fixture(scope="session")
def benchmark_calibration_path(fixtures_dir):
    return fixtures_dir / "benchmark" / "calibrations.yaml"


@pytest.fixture(scope="session")
def taxonomy_table_path(fixtures_dir):
    return fixtures_dir / "taxonomies" / "taxonomy.tsv"


@pytest.fixture(scope="session")
def empty_calibration_path(fixtures_dir):
    return fixtures_dir / "calibrations" / "empty.yaml"


@pytest.fixture(scope="session")
def default_config_path(fixtures_dir):
    return fixtures_dir / "configs" / "default.yaml"


# --------------------------------------------------------------------------- #
# Helper marker application
# --------------------------------------------------------------------------- #


def pytest_runtest_setup(item):
    """Skip tests whose required external binary is unavailable."""
    from tests.validation.helpers import external_available

    marker_map = {
        "requires_mcmctree": "mcmctree",
        "requires_treepl": "treepl",
        "requires_pathd8": "pathd8",
        "requires_lsd2": "lsd2",
        "requires_r8s": "r8s",
        "requires_wlogdate": "wlogdate",
        "requires_mdcat": "mdcat",
    }
    for marker, binary in marker_map.items():
        if item.get_closest_marker(marker) and not external_available(binary):
            pytest.skip(f"required external binary '{binary}' not available")
