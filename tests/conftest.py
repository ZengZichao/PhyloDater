"""
Pytest 配置和共享 fixtures

此文件包含所有测试共享的 fixtures 和配置
"""

import sys
from pathlib import Path

import pytest

# 添加项目根目录到路径
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture
def test_data_dir():
    """返回测试数据目录"""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def mock_tree_newick():
    """返回模拟的树 Newick 字符串"""
    return "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"


@pytest.fixture
def mock_alignment_fasta():
    """返回模拟的 FASTA 比对字符串"""
    return """>human
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chimp
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>mouse
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>rat
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""


@pytest.fixture
def mock_calibrations_yaml():
    """返回模拟的校准点配置 YAML 字符串"""
    return """
calibrations:
  - name: "Primates"
    mrca_pair: ["human", "chimp"]
    constraint:
      type: "uniform"
      min_age: 6.0
      max_age: 8.0

  - name: "Rodents"
    mrca_pair: ["mouse", "rat"]
    constraint:
      type: "soft_lower"
      min: 10.0

  - name: "Root"
    is_root: true
    constraint:
      type: "maximum"
      max: 100.0
"""


@pytest.fixture
def temp_work_dir(tmp_path):
    """创建临时工作目录"""
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    return work_dir


@pytest.fixture
def temp_output_dir(tmp_path):
    """创建临时输出目录"""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir


@pytest.fixture
def module_imports_cleanly():
    """在干净的解释器进程里导入一个模块，返回「它有没有拉起某模块」的探针。

    形如 ``assert "ete3" not in sys.modules`` 的断言检查的是**整个进程**的模块表，
    一旦同一批运行里先跑了别的测试（例如可视化测试真的导入了 ete3），断言就会
    因为测试顺序而不是产品行为失败。可选后端「不被模块导入阶段拉起」是产品属性，
    因此必须在新进程里验证。
    """

    def probe(module_name: str, dependency: str) -> bool:
        """Return True when ``module_name`` imports without pulling in ``dependency``."""
        import subprocess

        code = (
            "import sys\n"
            f"import {module_name}\n"
            f"sys.exit(1 if '{dependency}' in sys.modules else 0)\n"
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode in (
            0,
            1,
        ), f"could not import {module_name} in a clean interpreter: {result.stderr}"
        return result.returncode == 0

    return probe


@pytest.fixture(autouse=True)
def reset_registry():
    """
    在每个测试后重置 DatingMethodRegistry

    防止测试之间的状态污染
    """
    from phylodater.core.method_interface import DatingMethodRegistry

    # 保存原始状态
    original_registry = DatingMethodRegistry._registry.copy()

    yield

    # 恢复原始状态
    DatingMethodRegistry._registry.clear()
    DatingMethodRegistry._registry.update(original_registry)
