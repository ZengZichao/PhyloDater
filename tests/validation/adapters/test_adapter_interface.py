"""Tests for the DatingMethod abstract interface and registry."""

import pytest

from phylodater.core import DatingMethodRegistry
from phylodater.infrastructure import Configuration


class TestAdapterInterface:
    def test_all_methods_registered(self):
        methods = DatingMethodRegistry.list_methods()
        expected = {"mcmctree", "treepl", "pathd8", "lsd2", "r8s", "wlogdate", "mdcat"}
        assert expected.issubset(set(methods))

    @pytest.mark.requires_mcmctree
    def test_mcmctree_method_name(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "mcmctree", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.method_name == "mcmctree"

    @pytest.mark.requires_treepl
    def test_treepl_method_name(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "treepl", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.method_name == "treepl"

    @pytest.mark.requires_pathd8
    def test_pathd8_method_name(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "pathd8", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.method_name == "pathd8"

    @pytest.mark.requires_lsd2
    def test_lsd2_method_name(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "lsd2", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.method_name == "lsd2"

    @pytest.mark.requires_r8s
    def test_r8s_method_name(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "r8s", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.method_name == "r8s"

    @pytest.mark.requires_wlogdate
    def test_wlogdate_method_name(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "wlogdate", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.method_name == "wlogdate"

    @pytest.mark.requires_mdcat
    def test_mdcat_method_name(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "mdcat", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.method_name == "mdcat"
