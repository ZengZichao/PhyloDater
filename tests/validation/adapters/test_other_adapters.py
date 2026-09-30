"""Basic tests for non-MCMCTree adapters."""

import pytest

from phylodater.core import DatingMethodRegistry
from phylodater.infrastructure import Configuration
from phylodater.models import PhylogeneticTree


@pytest.mark.requires_pathd8
class TestPathD8Adapter:
    def test_validate_calibrations_accepts_maximum(self, tmp_test_dir):
        """PATHd8 can express a maximum age as a hard ``maxage`` directive, so
        ``validate_calibrations`` must NOT reject a MaximumAgeConstraint."""
        from phylodater.models import CalibrationPoint, MaximumAgeConstraint

        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "pathd8", config.get_config(), tmp_test_dir, config.software_paths
        )
        cal = CalibrationPoint(
            name="Root",
            age_constraint=MaximumAgeConstraint(max_age=100.0),
            is_root_node=True,
        )
        assert adapter.validate_calibrations([cal]) == []

    def test_prepare_inputs_requires_fixed_age(
        self, tmp_test_dir, small_tree_path, small_tree_alignment_path
    ):
        """PATHd8 refuses to prepare inputs unless at least one FIXAGE
        calibration is present (its fixage stage runs on every invocation)."""
        from phylodater.core.exceptions import CalibrationError
        from phylodater.models import CalibrationPoint, MaximumAgeConstraint

        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "pathd8", config.get_config(), tmp_test_dir, config.software_paths
        )
        tree = PhylogeneticTree.from_file(small_tree_path)
        cal = CalibrationPoint(
            name="Root",
            age_constraint=MaximumAgeConstraint(max_age=100.0),
            is_root_node=True,
        )
        with pytest.raises(CalibrationError, match="fixed age"):
            adapter.prepare_inputs(tree, [cal], small_tree_alignment_path)

    def test_prepare_inputs_with_fixed_age(
        self, tmp_test_dir, small_tree_path, small_tree_alignment_path
    ):
        """With at least one FIXAGE calibration the PATHd8 infile is generated
        and contains the fixage directive."""
        from phylodater.models import CalibrationPoint, FixedAgeConstraint

        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "pathd8", config.get_config(), tmp_test_dir, config.software_paths
        )
        tree = PhylogeneticTree.from_file(small_tree_path)
        cal = CalibrationPoint(
            name="Anchor",
            age_constraint=FixedAgeConstraint(fixed_age=10.0),
            is_root_node=True,
            resolved_taxa=tree.tip_names[:2],
            mrca_leaf_pair=tuple(tree.tip_names[:2]),
        )
        inputs = adapter.prepare_inputs(tree, [cal], small_tree_alignment_path)
        infile = inputs["infile"]
        text = infile.read_text(encoding="utf-8")
        assert "fixage=" in text
        adapter.cleanup(preserve_intermediates=False)


@pytest.mark.requires_lsd2
class TestLSD2Adapter:
    def test_validate_environment(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "lsd2", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.validate_environment() is True


@pytest.mark.requires_treepl
class TestTreePLAdapter:
    def test_validate_environment(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "treepl", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.validate_environment() is True


@pytest.mark.requires_r8s
class TestR8SAdapter:
    def test_validate_environment(self, tmp_test_dir):
        config = Configuration()
        adapter = DatingMethodRegistry.create(
            "r8s", config.get_config(), tmp_test_dir, config.software_paths
        )
        assert adapter.validate_environment() is True
