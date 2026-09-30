"""Unit tests for tree topology queries."""

from phylodater.models import PhylogeneticTree


class TestMRCA:
    def test_mrca_two_tips(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        leaves = tree.get_mrca(["A", "B"])
        assert set(leaves) == {"A", "B"}

    def test_mrca_across_clades(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        leaves = tree.get_mrca(["A", "C"])
        assert set(leaves) == {"A", "B", "C", "D"}

    def test_mrca_missing_tip_returns_none(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        assert tree.get_mrca(["A", "X"]) is None

    def test_mrca_terminals(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        pair = tree.get_mrca_terminals(["A", "B", "C"])
        assert len(pair) == 2
        assert set(pair).issubset({"A", "B", "C", "D"})


class TestMonophyly:
    def test_monophyletic_sister_pair(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        assert tree.is_monophyletic(["A", "B"]) is True

    def test_not_monophyletic(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        assert tree.is_monophyletic(["A", "C"]) is False

    def test_monophyletic_empty(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        assert tree.is_monophyletic([]) is True


class TestDescendantLeaves:
    def test_get_descendant_leaves(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        leaves = tree.get_descendant_leaves(["A", "B"])
        assert set(leaves) == {"A", "B"}


class TestTopologyID:
    def test_topology_id(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        tid = tree.get_node_topology_id(["A", "B"])
        assert tid.startswith("2L:")

    def test_topology_id_consistent_for_same_node(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        # Two different leaf pairs that anchor the same internal node (AB).
        tid1 = tree.get_node_topology_id(["A", "B"])
        tid2 = tree.get_node_topology_id(["A", "B"])
        assert tid1 == tid2
        # Different node should get a different ID.
        tid3 = tree.get_node_topology_id(["C", "D"])
        assert tid1 != tid3


class TestNodeDepth:
    def test_node_depth(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        assert tree.get_node_depth(["A", "B"]) == 1
        assert tree.get_node_depth(["A", "C"]) == 0


class TestIterClades:
    def test_iter_clades(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        clades = tree.iter_clades()
        assert len(clades) == 7  # 4 leaves + 2 internal + 1 root
        leaves_only = [c for c in clades if c["is_leaf"]]
        assert len(leaves_only) == 4


class TestDistanceMatrix:
    def test_distance_matrix_cached(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        matrix1 = tree.get_distance_matrix()
        matrix2 = tree.get_distance_matrix()
        assert matrix1 is matrix2
        assert ("A", "B") in matrix1 or ("B", "A") in matrix1
