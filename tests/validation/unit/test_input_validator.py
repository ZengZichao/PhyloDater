"""Unit tests for input validators."""

from phylodater.services import InputValidator
from phylodater.services.deep_validator import DeepValidator


class TestInputValidator:
    def test_validate_tree_file(self, small_tree_path):
        validator = InputValidator()
        result = validator.validate_tree_file(small_tree_path)
        assert result.is_valid is True
        assert "newick" in result.format_detected.lower()

    def test_validate_alignment_file(self, small_alignment_path):
        validator = InputValidator()
        result = validator.validate_alignment_file(small_alignment_path)
        assert result.is_valid is True
        assert "fasta" in result.format_detected.lower()

    def test_validate_calibration_file(self, valid_calibration_path):
        validator = InputValidator()
        result = validator.validate_calibration_file(valid_calibration_path)
        assert result.is_valid is True
        assert "yaml" in result.format_detected.lower()

    def test_validate_taxonomy_file(self, taxonomy_table_path):
        validator = InputValidator()
        result = validator.validate_taxonomy_file(taxonomy_table_path)
        assert result.is_valid is True


class TestDeepValidator:
    def test_validate_tree_deep(self, small_tree_path):
        validator = DeepValidator()
        result = validator.validate_tree_deep(small_tree_path)
        assert result.is_valid is True
        assert result.tree_count == 1

    def test_validate_alignment_deep(self, small_alignment_path):
        validator = DeepValidator()
        result = validator.validate_alignment_deep(small_alignment_path)
        assert result.is_valid is True
        assert result.num_sequences > 0

    def test_duplicate_ids(self, fixtures_dir):
        path = fixtures_dir / "alignments" / "duplicate_id.fasta"
        validator = DeepValidator()
        result = validator.validate_alignment_deep(path)
        assert result.is_valid is False
        assert len(result.duplicate_ids) > 0

    def test_malicious_tree(self, fixtures_dir):
        path = fixtures_dir / "trees" / "malicious_tree.nwk"
        validator = DeepValidator()
        result = validator.validate_tree_deep(path)
        assert result.is_valid is False
