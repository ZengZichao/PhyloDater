"""Functional tests for pipeline dry-run mode."""

from phylodater.core.pipeline import DatingPipeline, PipelineConfig
from phylodater.models import PhylogeneticTree
from phylodater.services import CalibrationLoader


class TestPipelineDryRun:
    def test_dry_run_mcmctree(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        tree = PhylogeneticTree.from_file(small_tree_path)
        calibrations = CalibrationLoader().load(small_tree_calibration_path)
        config = PipelineConfig(
            methods=["mcmctree"],
            output_dir=tmp_test_dir,
            dry_run=True,
        )
        pipeline = DatingPipeline(config)
        results = pipeline.run(tree, small_tree_alignment_path, calibrations)
        assert results == {}

    def test_dry_run_pathd8_without_fixed(
        self,
        small_tree_path,
        small_tree_alignment_path,
        small_tree_calibration_path,
        tmp_test_dir,
    ):
        tree = PhylogeneticTree.from_file(small_tree_path)
        calibrations = CalibrationLoader().load(small_tree_calibration_path)
        config = PipelineConfig(
            methods=["pathd8"],
            output_dir=tmp_test_dir,
            dry_run=True,
        )
        pipeline = DatingPipeline(config)
        results = pipeline.run(tree, small_tree_alignment_path, calibrations)
        assert results == {}
