"""
测试数据生成器

提供用于单元测试的模拟数据生成工具
"""

import tempfile
from pathlib import Path
from typing import List, Optional


def create_mock_tree_newick() -> str:
    """创建模拟的系统发育树（Newick 格式）"""
    return "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"


def create_mock_alignment_fasta() -> str:
    """创建模拟的序列比对（FASTA 格式）"""
    return """>human
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chimp
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>mouse
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>rat
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""


def create_mock_calibrations_yaml() -> str:
    """创建模拟的校准点配置（YAML 格式）。

    键名必须跟现行 schema（phylodater/services/calibration_loader.py 的
    ``CONSTRAINT_SPECS``，参考样本见 test-data/benchmark/calibrations.yaml）：
    uniform 的必填键是 ``min`` / ``max``，不是早期的 ``min_age`` / ``max_age``。
    夹具若写着旧键名，只会教出新的误用法：那些键会被当成"不认识的键"告警，
    然后因缺必填键而整条条目不予解析。
    """
    return """
calibrations:
  - name: "Primates"
    mrca_pair: ["human", "chimp"]
    constraint:
      type: "uniform"
      min: 6.0
      max: 8.0

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


def create_mock_config_yaml() -> str:
    """创建模拟的配置文件（YAML 格式）。

    只写现行配置层真正认识的键（否则 ``Configuration`` 会按
    ``Unknown parameter`` 告警丢弃，同样误导后来人）：

    * ``lsd2.model`` 默认已改成 ``None``，含义是"按比对类型自动选模型"
      （核酸 GTR+G / 蛋白 LG+G），所以下面不写死 ``LG+G``；
    * 线程数是 ``common.nthreads``，``lsd2`` 下没有 ``threads`` 这个字段。
    """
    return """
software_paths:
  paml_path: /opt/paml4.10.8
  iqtree_bin: /opt/iqtree2

software:
  common:
    nthreads: 4

  mcmctree:
    clock: 2
    num_runs: 2
    burnin: 1000
    nsample: 2000

  lsd2:
    model: null
"""


class MockDataBuilder:
    """
    测试数据构建器

    用于创建临时测试文件和目录
    """

    def __init__(self):
        self.temp_dir = Path(tempfile.mkdtemp(prefix="phylodater_test_"))
        self._files: List[Path] = []

    def create_tree_file(self, content: Optional[str] = None) -> Path:
        """创建树文件"""
        content = content or create_mock_tree_newick()
        file_path = self.temp_dir / "tree.nwk"
        file_path.write_text(content)
        self._files.append(file_path)
        return file_path

    def create_alignment_file(self, content: Optional[str] = None) -> Path:
        """创建比对文件"""
        content = content or create_mock_alignment_fasta()
        file_path = self.temp_dir / "alignment.fasta"
        file_path.write_text(content)
        self._files.append(file_path)
        return file_path

    def create_calibrations_file(self, content: Optional[str] = None) -> Path:
        """创建校准点配置文件"""
        content = content or create_mock_calibrations_yaml()
        file_path = self.temp_dir / "calibrations.yaml"
        file_path.write_text(content)
        self._files.append(file_path)
        return file_path

    def create_config_file(self, content: Optional[str] = None) -> Path:
        """创建配置文件"""
        content = content or create_mock_config_yaml()
        file_path = self.temp_dir / "config.yaml"
        file_path.write_text(content)
        self._files.append(file_path)
        return file_path

    def create_output_dir(self, name: str = "output") -> Path:
        """创建输出目录"""
        output_dir = self.temp_dir / name
        output_dir.mkdir(exist_ok=True)
        return output_dir

    def cleanup(self):
        """清理所有临时文件"""
        import shutil

        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
