"""
配置管理模块的单元测试
"""

import pytest
import yaml

from phylodater.infrastructure.configuration import (
    CommonConfig,
    Configuration,
    LSD2Config,
    MCMCTreeConfig,
    SoftwarePaths,
    ToolConfig,
    derive_method_seed,
)
from phylodater.infrastructure.safe_io import safe_writer


class TestSoftwarePaths:
    """测试 SoftwarePaths 数据类"""

    def test_default_creation(self):
        """测试默认创建"""
        paths = SoftwarePaths()
        assert paths.paml_path is None
        assert paths.iqtree_bin is None

    def test_custom_creation(self):
        """测试自定义路径创建"""
        paths = SoftwarePaths(paml_path="/opt/paml", iqtree_bin="/usr/bin/iqtree2")
        assert paths.paml_path == "/opt/paml"
        assert paths.iqtree_bin == "/usr/bin/iqtree2"


class TestCommonConfig:
    """测试 CommonConfig 数据类"""

    def test_default_values(self):
        """测试默认值"""
        config = CommonConfig()
        assert config.nthreads == 4
        assert config.seed is None
        assert config.verbose is False

    def test_custom_values(self):
        """测试自定义值"""
        config = CommonConfig(nthreads=8, seed=42, verbose=True)
        assert config.nthreads == 8
        assert config.seed == 42
        assert config.verbose is True


class TestMCMCTreeConfig:
    """测试 MCMCTreeConfig 数据类"""

    def test_default_values(self):
        """测试默认值"""
        config = MCMCTreeConfig()
        assert config.clock == 2
        assert config.num_runs == 2
        assert config.burnin == 20000
        assert config.nsample == 50000

    def test_custom_values(self):
        """测试自定义值"""
        config = MCMCTreeConfig(clock=3, num_runs=4, burnin=50000)
        assert config.clock == 3
        assert config.num_runs == 4
        assert config.burnin == 50000


class TestLSD2Config:
    """测试 LSD2Config 数据类"""

    def test_default_values(self):
        """测试默认值：未指定 model 时由比对类型自动选择"""
        config = LSD2Config()
        assert config.model is None
        assert config.iqtree_bin == "iqtree2"


class TestToolConfig:
    """测试 ToolConfig 数据类"""

    def test_default_creation(self):
        """测试默认创建包含所有工具配置"""
        config = ToolConfig()
        assert isinstance(config.common, CommonConfig)
        assert isinstance(config.mcmctree, MCMCTreeConfig)
        assert isinstance(config.lsd2, LSD2Config)

    def test_common_config_default(self):
        """测试 common 配置的默认值"""
        config = ToolConfig()
        assert config.common.nthreads == 4
        assert config.common.seed is None
        assert config.common.verbose is False


class TestMethodSeedDerivation:
    """修复项 E-3：方法级子种子派生 + ``ToolConfig.clone_for_method``"""

    METHODS = [
        "mcmctree",
        "lsd2",
        "treepl",
        "pathd8",
        "r8s",
        "wlogdate",
        "mdcat",
    ]

    def test_derive_is_deterministic(self):
        assert derive_method_seed(42, "mcmctree") == derive_method_seed(42, "mcmctree")

    def test_derive_is_a_32bit_int(self):
        value = derive_method_seed(42, "mcmctree")
        assert isinstance(value, int)
        assert 0 <= value < 2**32

    def test_derive_gives_pairwise_distinct_seeds(self):
        values = [derive_method_seed(42, name) for name in self.METHODS]
        assert len(set(values)) == len(values)

    def test_derive_changes_with_base_seed(self):
        assert (
            derive_method_seed(42, "treepl") != derive_method_seed(43, "treepl")
        )

    def test_derive_ignores_method_name_case(self):
        assert derive_method_seed(7, "MCMCTree") == derive_method_seed(7, "mcmctree")

    def test_derive_rejects_negative_base(self):
        with pytest.raises(ValueError, match="non-negative"):
            derive_method_seed(-1, "mcmctree")

    def test_clone_for_method_replaces_common_seed(self):
        base = ToolConfig()
        base.common.seed = 42
        base.common.nthreads = 3

        clone = base.clone_for_method("mcmctree", 42)

        assert clone.common.seed == derive_method_seed(42, "mcmctree")
        # 其余字段照搬，且原对象一个字节都没被改
        assert clone.common.nthreads == 3
        assert base.common.seed == 42
        assert clone is not base
        assert clone.common is not base.common

    def test_clone_for_method_is_deep_copy(self):
        base = ToolConfig()
        base.common.seed = 42
        clone = base.clone_for_method("treepl", 42)

        clone.treepl.initial_smooth = 0.5
        clone.common.nthreads = 99

        assert base.treepl.initial_smooth == 100.0
        assert base.common.nthreads == 4

    def test_clone_for_method_without_base_seed_keeps_seed(self):
        """base_seed=None 时不假装能派生：副本种子保持原样（公开 API 兼容）"""
        base = ToolConfig()
        base.common.seed = 42

        clone = base.clone_for_method("lsd2", None)

        assert clone.common.seed == 42
        assert clone.common is not base.common

    def test_clone_for_method_is_stable_across_batch_sizes(self):
        """同一方法在"3 方法并行"与"7 方法并行"下拿到同一个子种子。"""
        base = ToolConfig()
        small = base.clone_for_method("mdcat", 42).common.seed
        large = ToolConfig().clone_for_method("mdcat", 42).common.seed
        assert small == large


class TestLegacyYamlKeyAliases:
    """修复项 E-6：YAML 里 ``software.common.threads`` 与规范键 ``nthreads`` 的统一"""

    def _write(self, tmp_path, common_block: dict):
        path = tmp_path / "config.yaml"
        with safe_writer(path) as f:
            yaml.dump({"software": {"common": common_block}}, f)
        return path

    def test_canonical_key_silently_accepted(self, tmp_path):
        cfg = Configuration().load_from_file(self._write(tmp_path, {"nthreads": 9}))
        assert cfg.get_config().common.nthreads == 9

    def test_legacy_key_still_works(self, tmp_path):
        """旧模板/旧文档写的 threads 不能"写了白写"：仍生效 + 弃用告警。"""
        cfg = Configuration().load_from_file(self._write(tmp_path, {"threads": 7}))
        assert cfg.get_config().common.nthreads == 7

    def test_legacy_key_emits_deprecation_warning(self, tmp_path, capsys):
        Configuration().load_from_file(self._write(tmp_path, {"threads": 7}))
        captured = capsys.readouterr()
        combined = captured.out + captured.err
        assert "threads" in combined
        assert "nthreads" in combined

    def test_canonical_key_wins_when_both_present(self, tmp_path):
        cfg = Configuration().load_from_file(
            self._write(tmp_path, {"threads": 2, "nthreads": 11})
        )
        assert cfg.get_config().common.nthreads == 11

    def test_cli_override_beats_yaml_legacy_key(self, tmp_path):
        cfg = Configuration().load_from_file(self._write(tmp_path, {"threads": 7}))
        cfg.set_cli_override("common", "nthreads", 3)
        assert cfg.get_config().common.nthreads == 3

    def test_method_args_still_reject_legacy_key(self):
        """CLI 侧不做静默别名：--method-args threads=... 必须报错并给出建议。"""
        cfg = Configuration()
        with pytest.raises(ValueError, match="Unknown parameter"):
            cfg.set_cli_override("common", "threads", 4)

    def test_example_config_uses_canonical_key(self):
        """仓库里的 config.example.yaml 必须已经写规范键（不再自相矛盾）。"""
        from pathlib import Path

        example = Path(__file__).resolve().parents[3] / "config.example.yaml"
        if not example.exists():  # pragma: no cover - 打包安装后可能不存在
            pytest.skip("config.example.yaml not present")
        text = example.read_text(encoding="utf-8")
        assert "nthreads:" in text
        # 不再有"software.common.threads"这一旧写法
        assert "\n    threads:" not in text


class TestConfiguration:
    """测试 Configuration 管理器"""

    def test_init(self):
        """测试初始化"""
        config = Configuration()
        assert config.software_paths is not None
        assert isinstance(config.software_paths, SoftwarePaths)

    def test_set_software_path(self):
        """测试设置软件路径"""
        config = Configuration()
        config.set_software_path("paml_path", "/opt/paml4.10.8")
        assert config.software_paths.paml_path == "/opt/paml4.10.8"

    def test_set_unknown_software_path(self):
        """测试设置未知软件路径（应被忽略）"""
        config = Configuration()
        config.set_software_path("unknown_tool", "/path/to/tool")
        # 不应抛出异常，但也不会设置

    def test_load_from_file_not_exists(self, tmp_path, capsys):
        """测试加载不存在的配置文件"""
        config = Configuration()
        non_existent = tmp_path / "non_existent.yaml"

        config.load_from_file(non_existent)

        captured = capsys.readouterr()
        assert "Config file not found" in captured.out

    def test_load_from_file_valid(self, tmp_path):
        """测试加载有效的配置文件"""
        config_file = tmp_path / "config.yaml"
        config_data = {
            "software_paths": {"paml_path": "/opt/paml"},
            "software": {
                "common": {"nthreads": 8, "seed": 42},
                "mcmctree": {"clock": 3, "num_runs": 4},
            },
        }
        with safe_writer(config_file) as f:
            yaml.dump(config_data, f)

        config = Configuration()
        config.load_from_file(config_file)

        assert config.software_paths.paml_path == "/opt/paml"
        final_config = config.get_config()
        assert final_config.common.nthreads == 8
        assert final_config.common.seed == 42
        assert final_config.mcmctree.clock == 3
        assert final_config.mcmctree.num_runs == 4

    def test_cli_override(self):
        """测试命令行参数覆盖"""
        config = Configuration()
        config.set_cli_override("mcmctree", "clock", 3)

        final_config = config.get_config()
        assert final_config.mcmctree.clock == 3

    def test_cli_override_common(self):
        """测试 common 配置的命令行覆盖"""
        config = Configuration()
        config.set_cli_override("common", "nthreads", 16)
        config.set_cli_override("common", "seed", 123)

        final_config = config.get_config()
        assert final_config.common.nthreads == 16
        assert final_config.common.seed == 123

    def test_save_runtime_config(self, tmp_path):
        """测试保存运行时配置"""
        config = Configuration()
        output_path = tmp_path / "runtime_config.yaml"

        config.save_runtime_config(output_path)

        assert output_path.exists()
        with open(output_path, "r") as f:
            data = yaml.safe_load(f)
        assert "software" in data
        assert "common" in data["software"]
        assert "mcmctree" in data["software"]

    def test_create_default_config_file(self, tmp_path):
        """测试创建默认配置文件"""
        output_path = tmp_path / "default_config.yaml"

        Configuration.create_default_config_file(output_path)

        assert output_path.exists()
        with open(output_path, "r") as f:
            data = yaml.safe_load(f)
        assert "software" in data
        assert "common" in data["software"]
