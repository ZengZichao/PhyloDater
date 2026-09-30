"""Unit tests for tree rooting strategies."""

from phylodater.models import PhylogeneticTree


class TestOutgroupRooting:
    def test_outgroup_rooting(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        rooted, outgroup_leaves = tree.rerooted_with_outgroup(["C"])
        assert rooted is not None
        assert set(rooted.tip_names) == {"A", "B", "C", "D"}
        assert outgroup_leaves == {"C"}

    def test_outgroup_rooting_preserved_leaves(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        rooted, _ = tree.rerooted_with_outgroup(["C"])
        assert sorted(rooted.tip_names) == ["A", "B", "C", "D"]

    def test_outgroup_not_found(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        rooted, leaves = tree.rerooted_with_outgroup(["X"])
        assert rooted is None
        assert leaves == frozenset()


class TestMidpointRooting:
    def test_midpoint_rooting(self):
        tree = PhylogeneticTree("(((A:1,B:1):1,C:2):1,D:3);")
        rooted = tree.rerooted_at_midpoint()
        assert rooted is not None
        assert sorted(rooted.tip_names) == ["A", "B", "C", "D"]

    def test_midpoint_rooting_preserved_leaves(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        rooted = tree.rerooted_at_midpoint()
        assert rooted is not None
        assert sorted(rooted.tip_names) == ["A", "B", "C", "D"]

    def test_midpoint_on_single_tip_tree(self):
        tree = PhylogeneticTree("(A:0.1);")
        rooted = tree.rerooted_at_midpoint()
        assert rooted is None


class TestMinVarianceRooting:
    def test_min_variance_rooting(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        result = tree.rerooted_at_min_variance()
        assert result is not None
        rooted, variance = result
        assert rooted is not None
        assert variance >= 0
        assert sorted(rooted.tip_names) == ["A", "B", "C", "D"]


class TestRootedDigest:
    def test_digest_changes_after_rooting(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        before = tree._rooted_digest()
        rooted, _ = tree.rerooted_with_outgroup(["C"])
        after = rooted._rooted_digest()
        assert before != after
