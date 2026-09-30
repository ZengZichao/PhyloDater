"""Unit tests for PipelineConfig."""

import pytest

from phylodater.core.pipeline import PipelineConfig


class TestPipelineConfig:
    def test_default_values(self):
        cfg = PipelineConfig()
        assert cfg.methods == ["mcmctree", "treepl", "pathd8"]
        assert cfg.threads >= 1
        assert cfg.burnin >= 0
        assert cfg.nsample >= 1
        assert cfg.timeout >= 1

    def test_invalid_threads(self):
        with pytest.raises(ValueError):
            PipelineConfig(threads=0)

    def test_invalid_seed(self):
        with pytest.raises(ValueError):
            PipelineConfig(seed=-1)

    def test_invalid_burnin(self):
        with pytest.raises(ValueError):
            PipelineConfig(burnin=-1)

    def test_invalid_reroot_strategy(self):
        with pytest.raises(ValueError):
            PipelineConfig(reroot_strategy="invalid")

    def test_empty_methods(self):
        with pytest.raises(ValueError):
            PipelineConfig(methods=[])
