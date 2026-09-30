"""Unit tests for calibration constraint models."""

import pytest

from phylodater.models import (
    FixedAgeConstraint,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)


class TestFixedConstraint:
    def test_fixed_age(self):
        c = FixedAgeConstraint(fixed_age=10.0)
        rendered = c.to_mcmctree_calib_string()
        assert rendered.startswith("B(")
        assert "10" in rendered


class TestUniformConstraint:
    def test_uniform(self):
        c = UniformAgeConstraint(min_age=8.0, max_age=12.0)
        rendered = c.to_mcmctree_calib_string()
        assert rendered.startswith("B(")

    def test_uniform_with_tail_probs(self):
        c = UniformAgeConstraint(
            min_age=8.0, max_age=12.0, tail_lower=0.01, tail_upper=0.01
        )
        rendered = c.to_mcmctree_calib_string()
        assert rendered.startswith("B(")


class TestSoftLowerConstraint:
    def test_soft_lower(self):
        c = SoftLowerBoundConstraint(min_age=8.0)
        rendered = c.to_mcmctree_calib_string()
        assert rendered.startswith("L(")


class TestMaximumConstraint:
    def test_maximum(self):
        c = MaximumAgeConstraint(max_age=12.0)
        rendered = c.to_mcmctree_calib_string()
        assert rendered.startswith("U(")


class TestSoftBoundsConstraint:
    def test_soft_bounds(self):
        c = SoftBoundsConstraint(min_age=8.0, max_age=12.0)
        rendered = c.to_mcmctree_calib_string()
        assert rendered.startswith("B(")


class TestGammaConstraint:
    def test_gamma(self):
        c = GammaPriorConstraint(alpha=2.0, beta=1.0)
        rendered = c.to_mcmctree_calib_string()
        assert rendered.startswith("G(")


class TestInvalidConstraints:
    def test_uniform_min_greater_than_max(self):
        with pytest.raises((ValueError, AssertionError)):
            UniformAgeConstraint(min_age=12.0, max_age=8.0)

    def test_negative_age_rejected(self):
        with pytest.raises((ValueError, AssertionError)):
            FixedAgeConstraint(fixed_age=-5.0)
