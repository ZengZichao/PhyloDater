"""Unit tests for DatingMethodRegistry."""

import pytest

from phylodater.core import DatingMethodRegistry
from phylodater.core.exceptions import UnknownMethodError


class TestDatingMethodRegistry:
    def test_list_methods(self):
        methods = DatingMethodRegistry.list_methods()
        assert "mcmctree" in methods
        assert "treepl" in methods
        assert "pathd8" in methods
        assert "lsd2" in methods
        assert "r8s" in methods
        assert "wlogdate" in methods
        assert "mdcat" in methods

    def test_is_registered(self):
        assert DatingMethodRegistry.is_registered("mcmctree") is True
        assert DatingMethodRegistry.is_registered("unknown") is False

    def test_create_unknown_method(self, tmp_test_dir):
        from phylodater.infrastructure import Configuration

        config = Configuration()
        with pytest.raises(UnknownMethodError):
            DatingMethodRegistry.create(
                "nonexistent",
                config.get_config(),
                tmp_test_dir,
                config.software_paths,
            )
