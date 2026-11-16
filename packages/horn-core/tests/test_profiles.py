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
    scaled_tractrix_radius,
    cd_radius,
    hypex_radius,
    natural_tractrix_length,
    tractrix_constant,
    get_radius_func,
    ProfileGeometryError,
    CD_THROAT_FRACTION,
)


THROAT_R = 0.025
MOUTH_R = 0.15
LENGTH = 0.3
ALL_PROFILES = ["conical", "exponential", "hyperbolic", "tractrix", "os", "scaled_tractrix", "cd"]


class TestEndpoints:
    """Each profile must satisfy r(0) == throat_radius and r(L) == mouth_radius."""

    @pytest.fixture(params=ALL_PROFILES)
    def profile(self, request):
        return request.param

    def test_throat_radius(self, profile):
        func = get_radius_func(profile, THROAT_R, MOUTH_R, LENGTH)
        tol = 1e-3
        assert func(0.0) == pytest.approx(THROAT_R, rel=tol)

    def test_mouth_radius(self, profile):
        func = get_radius_func(profile, THROAT_R, MOUTH_R, LENGTH)
        tol = 1e-3
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


def _slope(func, z, h=1e-7):
    return (func(z + h) - func(z - h)) / (2 * h)


class TestTractrixGeometry:
    """The tractrix is the curve whose tangent, from wall to axis, has constant length a."""

    R_T, R_M = 0.0254, 0.182

    def _wall(self, L, n=400):
        func = get_radius_func("tractrix", self.R_T, self.R_M, L)
        z = np.linspace(0, L, n)
        return z, np.array([func(v) for v in z])

    @pytest.mark.parametrize("length_factor", [1.0, 1.3])
    def test_tangent_length_is_constant(self, length_factor):
        L = length_factor * natural_tractrix_length(self.R_T, self.R_M)
        z, r = self._wall(L)
        # Tangent length from finite differences of the sampled wall, away from
        # the mouth where dz/dr -> 0 makes differencing ill conditioned.
        keep = slice(2, -20)
        dz_dr = np.gradient(z, r)
        tangent = r * np.sqrt(1 + dz_dr**2)
        a = tangent[keep]
        assert np.ptp(a) / np.mean(a) < 0.01
        if length_factor == 1.0:
            assert np.mean(a) == pytest.approx(self.R_M, rel=0.01)
        else:
            assert np.mean(a) > self.R_M * 1.01

    def test_matches_closed_form_axial_position(self):
        """z(r) = x(r_t) - x(r), x(r) = a ln((a + sqrt(a^2 - r^2)) / r) - sqrt(a^2 - r^2)."""
        a = self.R_M
        L = natural_tractrix_length(self.R_T, self.R_M)

        def x(r):
            root = math.sqrt(a * a - r * r)
            return a * math.log((a + root) / r) - root

        func = get_radius_func("tractrix", self.R_T, self.R_M, L)
        for r in np.linspace(self.R_T * 1.01, self.R_M * 0.99, 25):
            z = x(self.R_T) - x(r)
            assert func(z) == pytest.approx(r, rel=1e-6)

    def test_natural_length_and_mouth_perpendicular(self):
        L = natural_tractrix_length(self.R_T, self.R_M)
        assert L == pytest.approx(0.30345, rel=1e-4)  # independent hand value for 25.4/182 mm
        func = get_radius_func("tractrix", self.R_T, self.R_M, L)
        assert func(0.0) == pytest.approx(self.R_T, rel=1e-9)
        assert func(L) == pytest.approx(self.R_M, rel=1e-9)
        # Wall turns perpendicular to the axis at the mouth: dr/dz grows without bound.
        near, nearer = _slope(func, L - 1e-5, 1e-6), _slope(func, L - 1e-7, 1e-8)
        assert near > 10 and nearer > 3 * near

    def test_flares_from_the_throat(self):
        """Regression: the old 'tractrix' stayed within 10% of the flare for 85% of its length."""
        z, r = self._wall(0.4)
        mid = r[len(r) // 2]
        assert (mid - self.R_T) / (self.R_M - self.R_T) > 0.15

    def test_volume_matches_closed_form(self):
        """V = (pi a^3 / 3) [(1 - r_t^2/a^2)^1.5 - (1 - r_m^2/a^2)^1.5]."""
        for L in (natural_tractrix_length(self.R_T, self.R_M), 0.45):
            a = tractrix_constant(self.R_T, self.R_M, L)
            expected = math.pi * a**3 / 3 * ((1 - (self.R_T / a) ** 2) ** 1.5 - (1 - (self.R_M / a) ** 2) ** 1.5)
            z, r = self._wall(L, n=20001)
            actual = np.trapezoid(np.pi * r**2, z)
            assert actual == pytest.approx(expected, rel=1e-4)

    def test_shorter_than_natural_length_is_rejected(self):
        with pytest.raises(ProfileGeometryError, match="natural length"):
            get_radius_func("tractrix", self.R_T, self.R_M, 0.2)


class TestOSGeometry:
    def test_throat_parallel_and_asymptotic_angle(self):
        func = get_radius_func("os", THROAT_R, MOUTH_R, LENGTH)
        assert _slope(func, 1e-6) == pytest.approx(0.0, abs=1e-4)
        tan_theta = math.sqrt(MOUTH_R**2 - THROAT_R**2) / LENGTH
        # Far from the throat the wall is a cone of half-angle theta.
        far = 1e4 * LENGTH
        assert _slope(func, far, 1e-3) == pytest.approx(tan_theta, rel=1e-6)


class TestHypex:
    @pytest.mark.parametrize("T", [0.0, 0.5, 0.7, 1.0, 2.0])
    def test_hits_mouth_radius(self, T):
        assert hypex_radius(LENGTH, THROAT_R, MOUTH_R, LENGTH, T) == pytest.approx(MOUTH_R, rel=1e-12)
        assert hypex_radius(0.0, THROAT_R, MOUTH_R, LENGTH, T) == pytest.approx(THROAT_R, rel=1e-12)

    def test_t1_is_exponential(self):
        for z in np.linspace(0, LENGTH, 31):
            assert hypex_radius(z, THROAT_R, MOUTH_R, LENGTH, 1.0) == pytest.approx(
                THROAT_R * (MOUTH_R / THROAT_R) ** (z / LENGTH), rel=1e-9)

    def test_t0_is_catenoidal_and_is_the_hyperbolic_default(self):
        m = math.acosh(MOUTH_R / THROAT_R)
        func = get_radius_func("hyperbolic", THROAT_R, MOUTH_R, LENGTH)
        for z in np.linspace(0, LENGTH, 31):
            expected = THROAT_R * math.cosh(m * z / LENGTH)
            assert hypex_radius(z, THROAT_R, MOUTH_R, LENGTH, 0.0) == pytest.approx(expected, rel=1e-12)
            assert func(z) == pytest.approx(expected, rel=1e-12)

    def test_throat_slope_follows_t(self):
        """dr/dz at the throat is T r_t / x0: zero for the catenoid, r_t/x0 for the exponential."""
        slopes = [_slope(get_radius_func("hyperbolic", THROAT_R, MOUTH_R, LENGTH, hypex_t=T), 1e-6)
                  for T in (0.0, 0.5, 1.0)]
        assert slopes[0] == pytest.approx(0.0, abs=1e-5)
        assert 0 < slopes[1] < slopes[2]
        assert slopes[2] == pytest.approx(THROAT_R * math.log(MOUTH_R / THROAT_R) / LENGTH, rel=1e-5)


class TestCD:
    @pytest.mark.parametrize("geometry", [(0.0254, 0.182, 0.4), (0.0125, 0.2, 0.5), (0.025, 0.15, 0.3)])
    def test_join_is_slope_continuous(self, geometry):
        r_t, r_m, L = geometry
        func = get_radius_func("cd", r_t, r_m, L)
        z_j = CD_THROAT_FRACTION * L
        h = 1e-6
        left = (func(z_j) - func(z_j - h)) / h
        right = (func(z_j + h) - func(z_j)) / h
        assert right == pytest.approx(left, rel=1e-4)
        assert func(0) == pytest.approx(r_t, rel=1e-12)
        assert func(L) == pytest.approx(r_m, rel=1e-9)

    def test_cone_is_straight_after_join(self):
        func = get_radius_func("cd", THROAT_R, MOUTH_R, LENGTH)
        z = np.linspace(CD_THROAT_FRACTION * LENGTH, LENGTH, 7)
        slopes = np.diff([func(v) for v in z]) / np.diff(z)
        assert np.ptp(slopes) < 1e-9


class TestScaledTractrixAndAlias:
    def test_lecleach_is_deprecated_alias(self):
        new = get_radius_func("scaled_tractrix", THROAT_R, MOUTH_R, LENGTH)
        with pytest.warns(DeprecationWarning, match="scaled_tractrix"):
            old = get_radius_func("lecleach", THROAT_R, MOUTH_R, LENGTH)
        for z in np.linspace(0, LENGTH, 11):
            assert old(z) == new(z)

    def test_equals_tractrix_only_at_natural_length(self):
        L_nat = natural_tractrix_length(THROAT_R, MOUTH_R)
        z = np.linspace(0, 0.98, 50)
        def max_diff(L):
            a = get_radius_func("scaled_tractrix", THROAT_R, MOUTH_R, L)
            b = get_radius_func("tractrix", THROAT_R, MOUTH_R, L)
            return max(abs(a(v * L) - b(v * L)) / b(v * L) for v in z)
        assert max_diff(L_nat) < 2e-3
        assert max_diff(1.5 * L_nat) > 0.05
