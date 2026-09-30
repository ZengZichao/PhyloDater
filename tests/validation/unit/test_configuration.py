"""Unit tests for Configuration."""

from phylodater.infrastructure import Configuration


class TestConfiguration:
    def test_default_config_loads(self):
        config = Configuration()
        cfg = config.get_config()
        assert cfg is not None

    def test_load_from_file(self, default_config_path):
        from pathlib import Path

        config = Configuration()
        config.load_from_file(Path(default_config_path))
        cfg = config.get_config()
        # Verify that a configuration object is populated after loading.
        assert cfg is not None
        assert cfg.mcmctree is not None
        assert cfg.common is not None

    def test_set_software_path(self):
        config = Configuration()
        config.set_software_path("iqtree_bin", "/usr/bin/iqtree")
        assert config.software_paths.iqtree_bin == "/usr/bin/iqtree"

    def test_set_cli_override(self):
        config = Configuration()
        config.set_cli_override("mcmctree", "burnin", "5000")
        cfg = config.get_config()
        assert cfg.mcmctree.burnin == 5000

    def test_env_threads_override(self, monkeypatch):
        monkeypatch.setenv("PHYLODATER_THREADS", "8")
        from phylodater.core.pipeline import get_env_threads

        assert get_env_threads(4) == 8

    def test_env_threads_invalid_ignored(self, monkeypatch):
        monkeypatch.setenv("PHYLODATER_THREADS", "abc")
        from phylodater.core.pipeline import get_env_threads

        assert get_env_threads(4) == 4

    def test_save_runtime_config(self, tmp_test_dir):
        config = Configuration()
        out = tmp_test_dir / "runtime.yaml"
        config.save_runtime_config(out)
        assert out.exists()
