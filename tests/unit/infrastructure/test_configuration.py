"""
配置管理模块的单元测试
"""

import yaml

from phylodater.infrastructure.configuration import (
    CommonConfig,
    Configuration,
    LSD2Config,
    MCMCTreeConfig,
    SoftwarePaths,
    ToolConfig,
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
