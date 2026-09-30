"""Unit tests for PipelineCheckpointManager."""

from phylodater.infrastructure.checkpoint import PipelineCheckpointManager


class TestPipelineCheckpointManager:
    def _make_result(self):
        from phylodater.models import DatingResult, NodeAgeEstimate

        result = DatingResult(
            method_name="mcmctree",
            run_id="run_1",
            dated_tree_newick="(A:1,B:1);",
        )
        result.execution_seconds = 1.0
        result.node_ages = {"A": NodeAgeEstimate(mean_age=1.0)}
        return result

    def test_start_complete_method(self, tmp_test_dir):
        mgr = PipelineCheckpointManager(tmp_test_dir)
        mgr.start_method("mcmctree")
        mgr.complete_method("mcmctree", result=self._make_result())
        assert mgr.is_method_completed("mcmctree") is True

    def test_fail_method(self, tmp_test_dir):
        mgr = PipelineCheckpointManager(tmp_test_dir)
        mgr.start_method("treepl")
        mgr.fail_method("treepl", "Environment not available")
        assert mgr.is_method_failed("treepl") is True

    def test_get_methods_to_run(self, tmp_test_dir):
        mgr = PipelineCheckpointManager(tmp_test_dir)
        mgr.start_method("mcmctree")
        mgr.complete_method("mcmctree", result=self._make_result())
        methods = mgr.get_methods_to_run(["mcmctree", "treepl"], skip_completed=True)
        assert "mcmctree" not in methods
        assert "treepl" in methods

    def test_input_hash_roundtrip(self, tmp_test_dir):
        mgr = PipelineCheckpointManager(tmp_test_dir)
        mgr.record_input_hash("tree", "abc123")
        assert mgr.get_input_hash("tree") == "abc123"

    def test_input_matches(self, tmp_test_dir):
        mgr = PipelineCheckpointManager(tmp_test_dir)
        key = (
            mgr.RUN_FINGERPRINT_KEY
            if hasattr(mgr, "RUN_FINGERPRINT_KEY")
            else "phylodater:run-fingerprint"
        )
        mgr.record_input_hash(key, "fingerprint_1")
        assert mgr.input_matches("fingerprint_1") is True
        assert mgr.input_matches("fingerprint_2") is False
