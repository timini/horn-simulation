"""Tests for pure-Python horn profile radius functions."""

import math

import numpy as np
import pytest

from horn_core.profiles import (
    conical_radius,
    exponential_radius,
    hyperbolic_radius,
    tractrix_radius,
    os_radius,
    lecleach_radius,
    cd_radius,
    get_radius_func,
)


THROAT_R = 0.025
MOUTH_R = 0.15
LENGTH = 0.3
ALL_PROFILES = ["conical", "exponential", "hyperbolic", "tractrix", "os", "lecleach", "cd"]


class TestEndpoints:
    """Each profile must satisfy r(0) == throat_radius and r(L) == mouth_radius."""

    @pytest.fixture(params=ALL_PROFILES)
    def profile(self, request):
        return request.param

    def test_throat_radius(self, profile):
        func = get_radius_func(profile, THROAT_R, MOUTH_R, LENGTH)
        # Le Cléac'h clips the tractrix curve, so throat radius is approximate
        tol = 0.02 if profile == "lecleach" else 1e-3
        assert func(0.0) == pytest.approx(THROAT_R, rel=tol)

    def test_mouth_radius(self, profile):
        func = get_radius_func(profile, THROAT_R, MOUTH_R, LENGTH)
        tol = 0.02 if profile == "lecleach" else 1e-3
        assert func(LENGTH) == pytest.approx(MOUTH_R, rel=tol)


class TestMonotonicity:
    """Horn profiles should expand monotonically (r(z2) >= r(z1) for z2 > z1)."""

    @pytest.fixture(params=ALL_PROFILES)
    def profile(self, request):
        return request.param

    def test_monotonic_expansion(self, profile):
        func = get_radius_func(profile, THROAT_R, MOUTH_R, LENGTH)
        z_vals = np.linspace(0, LENGTH, 200)
        r_vals = [func(z) for z in z_vals]
        for i in range(1, len(r_vals)):
            assert r_vals[i] >= r_vals[i - 1] - 1e-10, (
                f"{profile}: non-monotonic at z={z_vals[i]:.4f}, "
                f"r[{i-1}]={r_vals[i-1]:.6f} > r[{i}]={r_vals[i]:.6f}"
            )


class TestConicalAnalytical:
    """Conical profile has exact analytical form."""

    def test_midpoint(self):
        z_mid = LENGTH / 2
        expected = THROAT_R + (MOUTH_R - THROAT_R) / 2
        assert conical_radius(z_mid, THROAT_R, MOUTH_R, LENGTH) == pytest.approx(expected)

    def test_quarter_point(self):
        z = LENGTH / 4
        expected = THROAT_R + (MOUTH_R - THROAT_R) / 4
        assert conical_radius(z, THROAT_R, MOUTH_R, LENGTH) == pytest.approx(expected)


class TestExponentialAnalytical:
    """Exponential profile: r(z) = r_t * exp(m*z/L)."""

    def test_midpoint(self):
        m = math.log(MOUTH_R / THROAT_R)
        z_mid = LENGTH / 2
        expected = THROAT_R * math.exp(m * 0.5)
        assert exponential_radius(z_mid, THROAT_R, MOUTH_R, LENGTH) == pytest.approx(expected)


class TestHyperbolicAnalytical:
    """Hyperbolic profile: r(z) = r_t * cosh(m*z/L)."""

    def test_midpoint(self):
        m = np.arccosh(MOUTH_R / THROAT_R)
        z_mid = LENGTH / 2
        expected = float(THROAT_R * np.cosh(m * 0.5))
        assert hyperbolic_radius(z_mid, THROAT_R, MOUTH_R, LENGTH) == pytest.approx(expected, rel=1e-6)


class TestGetRadiusFunc:
    def test_invalid_profile_raises(self):
        with pytest.raises(ValueError, match="Unknown profile"):
            get_radius_func("parabolic", THROAT_R, MOUTH_R, LENGTH)

    def test_closure_returns_float(self):
        for profile in ALL_PROFILES:
            func = get_radius_func(profile, THROAT_R, MOUTH_R, LENGTH)
            result = func(LENGTH / 2)
            assert isinstance(result, float)

    def test_different_geometries(self):
        """Different throat/mouth should give different mid-section radii."""
        f1 = get_radius_func("exponential", 0.01, 0.1, 0.3)
        f2 = get_radius_func("exponential", 0.01, 0.2, 0.3)
        assert f1(0.15) < f2(0.15)


class TestCrossValidation:
    """Verify profiles.py matches horn-geometry/generator.py radius functions."""

    def test_conical_matches_generator(self):
        """Conical is trivial — exact match expected."""
        from horn_core.profiles import conical_radius
        for z in np.linspace(0, LENGTH, 50):
            expected = THROAT_R + (MOUTH_R - THROAT_R) * z / LENGTH
            assert conical_radius(z, THROAT_R, MOUTH_R, LENGTH) == pytest.approx(expected, abs=1e-12)

    def test_exponential_matches_generator(self):
        from horn_core.profiles import exponential_radius
        m = np.log(MOUTH_R / THROAT_R)
        for z in np.linspace(0, LENGTH, 50):
            expected = THROAT_R * np.exp(m * z / LENGTH)
            assert exponential_radius(z, THROAT_R, MOUTH_R, LENGTH) == pytest.approx(float(expected), rel=1e-10)

    def test_os_matches_generator(self):
        from horn_core.profiles import os_radius
        theta = math.atan2(math.sqrt(MOUTH_R**2 - THROAT_R**2), LENGTH)
        for z in np.linspace(0, LENGTH, 50):
            expected = math.sqrt(THROAT_R**2 + (z * math.tan(theta))**2)
            assert os_radius(z, THROAT_R, MOUTH_R, LENGTH) == pytest.approx(expected, rel=1e-10)

    def test_cd_matches_generator(self):
        from horn_core.profiles import cd_radius
        frac = 0.3
        z_t = frac * LENGTH
        r_trans = THROAT_R * (MOUTH_R / THROAT_R) ** frac
        for z in np.linspace(0, LENGTH, 50):
            if z <= z_t:
                expected = THROAT_R * math.exp(math.log(r_trans / THROAT_R) * z / z_t)
            else:
                expected = r_trans + (MOUTH_R - r_trans) * (z - z_t) / (LENGTH - z_t)
            assert cd_radius(z, THROAT_R, MOUTH_R, LENGTH) == pytest.approx(expected, rel=1e-10)
