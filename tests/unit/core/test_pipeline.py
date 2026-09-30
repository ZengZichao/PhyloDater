"""Regression tests for parallel per-method sub-seed derivation.

Covers review-report item **P0-1** ("seed derivation degenerates to the same
value across the seven methods"): ``SeedSequence(base).spawn(n)`` children all
report ``.entropy == base``, so reading ``.entropy`` collapsed every parallel
method onto one shared seed. ``_derive_method_seeds`` must instead draw each
child's own stream so the seeds are distinct yet reproducible.
"""

from phylodater.core.pipeline import _derive_method_seeds

ALL_METHODS = [
    "mcmctree",
    "lsd2",
    "treepl",
    "pathd8",
    "r8s",
    "wlogdate",
    "mdcat",
]


def test_every_method_is_assigned():
    seeds = _derive_method_seeds(42, ALL_METHODS)
    assert set(seeds) == set(ALL_METHODS)


def test_sub_seeds_are_pairwise_distinct():
    """The core P0-1 assertion: no two of the seven methods share a seed."""
    seeds = _derive_method_seeds(42, ALL_METHODS)
    values = list(seeds.values())
    assert len(set(values)) == len(
        values
    ), "P0-1 regression: parallel methods collapsed onto a shared seed"


def test_sub_seeds_are_reproducible_for_fixed_base():
    a = _derive_method_seeds(42, ALL_METHODS)
    b = _derive_method_seeds(42, ALL_METHODS)
    assert a == b


def test_different_base_seed_changes_streams():
    a = _derive_method_seeds(42, ["mcmctree", "lsd2"])
    b = _derive_method_seeds(1234, ["mcmctree", "lsd2"])
    assert a != b


def test_single_and_empty_method_lists():
    one = _derive_method_seeds(7, ["mcmctree"])
    assert list(one) == ["mcmctree"]
    assert isinstance(one["mcmctree"], int)
    assert _derive_method_seeds(7, []) == {}
