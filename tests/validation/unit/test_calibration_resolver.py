"""Unit tests for CalibrationResolver."""

import pytest

from phylodater.models import PhylogeneticTree
from phylodater.services import CalibrationResolver


class TestCalibrationResolver:
    def test_resolve_mrca_pair(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        resolver = CalibrationResolver(tree)
        result = resolver.resolve_mrca_pair(["A", "B"])
        assert result is not None
        assert set(result.resolved_taxa) == {"A", "B"}

    def test_resolve_luca(self):
        tree = PhylogeneticTree(
            "((Bact1_d_Bacteria:0.1,Bact2_d_Bacteria:0.1):0.2,(Arch1_d_Archaea:0.3,Arch2_d_Archaea:0.3):0.1);"
        )
        resolver = CalibrationResolver(tree)
        result = resolver.resolve("LUCA")
        assert result is not None
        assert len(result.resolved_taxa) == 4

    def test_resolve_root(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        resolver = CalibrationResolver(tree)
        result = resolver.resolve("ROOT")
        assert result is not None
        assert result.is_root_node is True

    def test_validate_ancestry_no_conflict(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        resolver = CalibrationResolver(tree)
        from phylodater.models import CalibrationPoint, UniformAgeConstraint

        cal1 = CalibrationPoint(
            name="AB",
            age_constraint=UniformAgeConstraint(min_age=1.0, max_age=2.0),
            resolved_taxa=["A", "B"],
            mrca_leaf_pair=("A", "B"),
        )
        cal2 = CalibrationPoint(
            name="CD",
            age_constraint=UniformAgeConstraint(min_age=1.0, max_age=2.0),
            resolved_taxa=["C", "D"],
            mrca_leaf_pair=("C", "D"),
        )
        resolver.validate_ancestry([cal1, cal2])

    def test_ancestor_conflict(self):
        tree = PhylogeneticTree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        resolver = CalibrationResolver(tree)
        from phylodater.models import CalibrationPoint, UniformAgeConstraint

        cal1 = CalibrationPoint(
            name="Root",
            age_constraint=UniformAgeConstraint(min_age=10.0, max_age=20.0),
            resolved_taxa=["A", "B", "C", "D"],
            mrca_leaf_pair=("A", "C"),
        )
        cal2 = CalibrationPoint(
            name="AB",
            age_constraint=UniformAgeConstraint(min_age=30.0, max_age=40.0),
            resolved_taxa=["A", "B"],
            mrca_leaf_pair=("A", "B"),
        )
        with pytest.raises(Exception):
            resolver.validate_ancestry([cal1, cal2])
