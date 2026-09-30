"""Tests for the MCMCTree adapter."""

import pytest

from phylodater.core import DatingMethodRegistry
from phylodater.infrastructure import Configuration
from phylodater.models import PhylogeneticTree
from phylodater.services import CalibrationLoader


@pytest.mark.requires_mcmctree
class TestMCMCTreeAdapter:
    def test_validate_environment(self, tmp_test_dir):
        config = Configuration()
        # Allow any detected PAML version rather than enforcing the default 4.10.8.
        cfg = config.get_config()
        cfg.mcmctree.paml_version = None
        adapter = DatingMethodRegistry.create(
            "mcmctree", cfg, tmp_test_dir, config.software_paths
        )
        assert adapter.validate_environment() is True

    def test_validate_calibrations_fixed(self, tmp_test_dir, small_tree_path):
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "mcmctree", config.get_config(), tmp_test_dir, config.software_paths
        )
        cal = CalibrationPoint(
            name="AB",
            age_constraint=FixedAgeConstraint(fixed_age=10.0),
            resolved_taxa=["A", "B"],
            mrca_leaf_pair=("A", "B"),
        )
        errors = adapter.validate_calibrations([cal])
        assert not errors

    def test_prepare_inputs(
        self,
        tmp_test_dir,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
    ):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "mcmctree", config.get_config(), tmp_test_dir, config.software_paths
        )
        tree = PhylogeneticTree.from_file(small_tree_path)
        calibrations = CalibrationLoader().load(small_tree_calibration_path)
        adapter.prepare_inputs(tree, calibrations, small_tree_alignment_path)
        # input files should be recorded
        assert adapter._input_files
        adapter.cleanup(preserve_intermediates=False)
