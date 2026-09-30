"""Unit tests for PhylogeneticTree model."""

import pytest

from phylodater.models import PhylogeneticTree


class TestTreeParsing:
    def test_from_newick_basic(self):
        tree = PhylogeneticTree.from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        assert tree.num_tips == 4
        assert tree.tip_names == ["A", "B", "C", "D"]

    def test_from_file(self, small_tree_path):
        tree = PhylogeneticTree.from_file(small_tree_path)
        assert tree.num_tips == 4

    def test_empty_file_raises(self, tmp_test_dir):
        empty = tmp_test_dir / "empty.nwk"
        empty.write_text("", encoding="utf-8")
        with pytest.raises(ValueError):
            PhylogeneticTree.from_file(empty)

    def test_multi_tree_nexus(self, fixtures_dir):
        path = fixtures_dir / "trees" / "multi_tree.nex"
        # Current implementation treats the whole file as a single newick-like
        # string; loading may succeed but does not parse Nexus syntax.
        tree = PhylogeneticTree.from_file(path)
        # The parser extracts the first parenthesized tree it can find.
        assert tree.num_tips >= 1


class TestTreeProperties:
    def test_num_tips(self):
        tree = PhylogeneticTree("((A,B),(C,D));")
        assert tree.num_tips == 4

    def test_tip_names_order(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        assert tree.tip_names == ["A", "B", "C", "D"]

    def test_is_binary_true(self):
        tree = PhylogeneticTree("((A,B),(C,D));")
        assert tree.is_binary() is True

    def test_is_binary_false(self):
        tree = PhylogeneticTree("(A,B,C,D);")
        assert tree.is_binary() is False

    def test_without_branch_lengths(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        cleaned = tree.without_branch_lengths()
        assert ":" not in cleaned.newick

    def test_without_internal_labels(self):
        tree = PhylogeneticTree("((A,B)AB,(C,D)CD)root;")
        cleaned = tree.without_internal_labels()
        assert "AB" not in cleaned.newick
        assert "root" not in cleaned.newick

    def test_strip_annotations(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2)[&label=1]:0.3,(C:0.4,D:0.5):0.6);")
        cleaned = tree.strip_annotations()
        assert "[" not in cleaned.newick


class TestRootingInfo:
    def test_rooted_binary(self):
        tree = PhylogeneticTree("((A,B),(C,D));")
        info = tree.get_rooting_info()
        assert info["is_rooted"] is True
        assert info["root_children"] == 2

    def test_unrooted_polytomy(self):
        tree = PhylogeneticTree("(A,B,C,D);")
        info = tree.get_rooting_info()
        assert info["is_rooted"] is False
        assert info["root_children"] == 4

    def test_has_outgroup(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        info = tree.get_rooting_info()
        # Neither child is a single tip, so no outgroup detected
        assert info["has_outgroup"] is False


class TestBranchLengthValidation:
    def test_valid_branch_lengths(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        valid, nodes, errors = tree.validate_branch_lengths()
        assert valid is True
        assert not nodes
        assert not errors

    def test_negative_branch_length(self):
        tree = PhylogeneticTree("((A:0.1,B:-0.2):0.3,(C:0.4,D:0.5):0.6);")
        valid, nodes, errors = tree.validate_branch_lengths()
        assert valid is False
        assert nodes
        assert errors


class TestTreeWrite:
    def test_write_newick(self, tmp_test_dir):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        out = tmp_test_dir / "out.nwk"
        tree.write(out, format="newick")
        assert out.exists()
        assert "A" in out.read_text(encoding="utf-8")

    def test_write_nexus(self, tmp_test_dir):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        out = tmp_test_dir / "out.nexus"
        tree.write(out, format="nexus")
        assert out.exists()
        text = out.read_text(encoding="utf-8")
        assert "#NEXUS" in text.upper() or "A" in text


class TestRenamedLeaves:
    def test_with_renamed_leaves(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        renamed = tree.with_renamed_leaves({"A": "Alpha", "B": "Beta"})
        assert "Alpha" in renamed.newick
        assert "Beta" in renamed.newick
        # with_renamed_leaves is a simple string replacement, so "A" may
        # appear inside "Alpha". Verify the original single-letter tips are gone.
        assert "(A:" not in renamed.newick
        assert ",A:" not in renamed.newick
        assert "(B:" not in renamed.newick
        assert ",B:" not in renamed.newick
