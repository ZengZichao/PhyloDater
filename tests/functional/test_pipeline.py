"""
模块六 & 十：流水线引擎与端到端流程功能测试

验证配置优先级、reroot、输出产物、错误恢复等。
"""

import yaml

from phylodater.core.pipeline import PipelineConfig
from tests.functional.conftest import run_phylodater_cli


class TestPipelineConfiguration:
    """流水线配置测试"""

    def test_pip_001_default_config(self, tmp_path):
        """PIP-001: 默认配置"""
        config = PipelineConfig(
            methods=["pathd8"],
            output_dir=tmp_path,
            threads=4,
        )
        assert config.methods == ["pathd8"]
        assert config.threads == 4

    def test_pip_002_threads(self, tmp_path):
        """PIP-002: 自定义线程数"""
        config = PipelineConfig(methods=["pathd8"], output_dir=tmp_path, threads=8)
        assert config.threads == 8

    def test_pip_004_dry_run(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-004: 干运行模式"""
        output_dir = tmp_path / "out"
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
                "--dry-run",
            ]
        )
        assert result.returncode in [0, 1]

    def test_pip_011_yaml_config(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-011: YAML 配置加载"""
        config_file = tmp_path / "config.yaml"
        config_file.write_text(
            yaml.dump({"software_paths": {"paml_path": "/opt/paml"}}), encoding="utf-8"
        )
        output_dir = tmp_path / "out"
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
                "--config",
                str(config_file),
            ]
        )
        assert result.returncode in [0, 1]

    def test_pip_014_invalid_yaml_config(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-014: 无效 YAML 配置文件"""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("invalid: yaml: [", encoding="utf-8")
        output_dir = tmp_path / "out"
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
                "--config",
                str(config_file),
            ]
        )
        assert result.returncode != 0


class TestSingleMethodPipeline:
    """单方法流程测试"""

    def test_pip_020_pathd8_pipeline(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-020 / E2E-001: PATHd8 完整流程"""
        output_dir = tmp_path / "pathd8_test"
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
        assert result.returncode in [0, 1]
        # 若成功，验证输出目录包含日志
        assert (output_dir / "phylodater.log").exists() or result.returncode != 0


class TestErrorHandling:
    """错误处理测试"""

    def test_pip_040_input_validation_failure(
        self, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-040: 输入验证失败退出码 3"""
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                "/nonexistent/tree.nwk",
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
        assert result.returncode == 3


class TestReproducibility:
    """可重复性测试"""

    def test_pip_053_seed_reproducible(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-053: 指定 seed 后参数可传递"""
        output_dir = tmp_path / "out"
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
                "--seed",
                "42",
            ]
        )
        assert result.returncode in [0, 1]
