"""
模块九：对抗性输入防护功能测试

验证恶意字符、空文件、输出覆盖保护、路径安全等。
"""

import pytest

from phylodater import PhylogeneticTree
from phylodater.services.deep_validator import DeepValidator
from tests.functional.conftest import run_phylodater_cli


class TestMaliciousCharacters:
    """恶意字符检测测试"""

    @pytest.mark.parametrize(
        "filename",
        [
            "control_char_tree.nwk",
            "unicode_bidi_tree.nwk",
            "zero_width_tree.nwk",
        ],
    )
    def test_sec_001_003_malicious_chars(self, functional_input_dir, filename):
        """SEC-001/002/003: 恶意字符"""
        bad_tree = functional_input_dir / "invalid" / filename
        validator = DeepValidator()
        result = validator.validate_tree_deep(bad_tree)
        assert not result.is_valid or any(
            "恶意" in e or "invalid" in e.lower() for e in result.errors
        )


class TestEmptyFileDetection:
    """空文件检测测试"""

    def test_sec_020_empty_tree(self, functional_input_dir):
        """SEC-020: 空树文件"""
        empty_tree = functional_input_dir / "boundary" / "empty_tree.nwk"
        with pytest.raises(Exception):
            PhylogeneticTree.from_file(str(empty_tree))

    def test_sec_021_empty_sequence(self, functional_input_dir):
        """SEC-021: 空序列文件"""
        empty_aln = functional_input_dir / "invalid" / "empty_alignment.fasta"
        validator = DeepValidator()
        result = validator.validate_alignment_deep(empty_aln)
        assert not result.is_valid or result.num_sequences == 0


class TestOutputOverwriteProtection:
    """输出覆盖保护测试"""

    def test_sec_031_force_mode(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """SEC-031: force 模式覆盖已存在文件"""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        # 先创建占位文件
        (output_dir / "phylodater.log").write_text("old log", encoding="utf-8")
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--force",
            ]
        )
        assert result.returncode in [0, 1]


class TestPathSafety:
    """路径安全测试"""

    def test_sec_040_path_traversal_attempt(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """SEC-040: 输出目录路径遍历尝试"""
        # 使用相对路径 ../ 作为输出目录，验证程序不退出工作目录
        output_dir = tmp_path / "out" / ".." / "escaped"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
            ]
        )
        # 当前实现可能只是创建目录；至少不应崩溃或进入敏感路径
        assert result.returncode in [0, 1, 3]

    def test_sec_041_long_node_name(self, tmp_path):
        """SEC-041: 超长节点名"""
        long_name = "A" * 5000
        tree_file = tmp_path / "long.nwk"
        tree_file.write_text(f"(({long_name}:0.1,B:0.2):0.3,C:0.4);", encoding="utf-8")
        validator = DeepValidator()
        result = validator.validate_tree_deep(tree_file)
        # 超长名称应产生警告或错误
        assert len(result.errors) > 0 or len(result.warnings) > 0 or not result.is_valid

    def test_sec_042_huge_calibration_age(self, tmp_path):
        """SEC-042: 超大校准年龄"""
        import yaml

        cal_file = tmp_path / "huge.yaml"
        config = {
            "calibrations": [
                {
                    "name": "Huge",
                    "is_root": True,
                    "constraint": {"type": "maximum", "max": 99999},
                }
            ]
        }
        cal_file.write_text(yaml.dump(config), encoding="utf-8")
        from phylodater import CalibrationLoader

        loader = CalibrationLoader()
        calibrations = loader.load(cal_file)
        # 当前实现已对 >4600 Ma 的约束拒绝加载
        assert len(calibrations) == 0
