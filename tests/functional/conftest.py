"""
功能流程测试共享 fixtures

提供测试数据路径、CLI 运行辅助等共享资源。
"""

import subprocess
import sys
from pathlib import Path

import pytest

# 功能流程测试根目录(实际数据位于项目根目录的 test-data/)
#   test-data/input/valid/   ← 正常数据
#   test-data/input/invalid/  ← 异常数据
#   test-data/input/boundary/ ← 边界数据
#   test-data/output/         ← 期望输出模板/测试输出目录
FUNCTIONAL_TEST_DIR = Path(__file__).parent.parent.parent / "test-data"
INPUT_DIR = FUNCTIONAL_TEST_DIR / "input"
OUTPUT_DIR = FUNCTIONAL_TEST_DIR / "output"


@pytest.fixture(scope="session")
def functional_input_dir() -> Path:
    """返回功能流程测试输入数据目录"""
    return INPUT_DIR


@pytest.fixture(scope="session")
def functional_output_dir() -> Path:
    """返回功能流程测试输出目录"""
    OUTPUT_DIR.mkdir(exist_ok=True)
    return OUTPUT_DIR


@pytest.fixture
def normal_tree(functional_input_dir: Path) -> Path:
    return functional_input_dir / "valid" / "tree.nwk"


@pytest.fixture
def normal_alignment(functional_input_dir: Path) -> Path:
    return functional_input_dir / "valid" / "alignment.fasta"


@pytest.fixture
def normal_calibration(functional_input_dir: Path) -> Path:
    return functional_input_dir / "valid" / "calibrations.yaml"


def run_phylodater_cli(args: list, timeout: int = 60) -> subprocess.CompletedProcess:
    """运行 phylodater CLI 命令并返回结果"""
    cmd = [sys.executable, "-m", "phylodater"] + args
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
