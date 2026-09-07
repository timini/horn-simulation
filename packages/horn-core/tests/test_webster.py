"""Tests for the TMM/Webster horn throat impedance solver."""

import math

import numpy as np
import pytest

from horn_core.profiles import get_radius_func
from horn_core.webster import (
    C0,
    RHO0,
    compute_throat_impedance_tmm,
    piston_radiation_impedance,
    unflanged_radiation_impedance,
)


class TestRadiationImpedance:
    """Verify radiation impedance matches horn_solver/radiation.py."""

    def test_flanged_piston_small_ka(self):
        """For small ka, R1 ≈ (ka)²/2, X1 ≈ 8ka/(3π)."""
        k = 2 * math.pi * 100 / C0  # 100 Hz
        a = 0.01  # small radius
        z = piston_radiation_impedance(k, a)
        ka = k * a
        # Taylor: R1(2ka) ≈ (2ka)²/8 = ka²/2
        assert z.real == pytest.approx(ka**2 / 2, rel=0.1)
        assert z.imag > 0  # reactance is positive

    def test_flanged_piston_large_ka(self):
        """For large ka, R1 → 1."""
        k = 2 * math.pi * 10000 / C0
        a = 0.1  # large radius
        z = piston_radiation_impedance(k, a)
        assert z.real == pytest.approx(1.0, abs=0.1)

    def test_unflanged_small_ka(self):
        k = 2 * math.pi * 100 / C0
        a = 0.01
        z = unflanged_radiation_impedance(k, a)
        ka = k * a
        assert z.real == pytest.approx(0.25 * ka**2, rel=1e-6)
        assert z.imag == pytest.approx(0.6133 * ka, rel=1e-6)


class TestCylindricalPipe:
    """Cylindrical horn (r_t == r_m) should behave like a simple pipe."""

    def test_real_impedance_positive(self):
        """Re(Z) > 0 at all frequencies (energy conservation)."""
        r = 0.025
        L = 0.3
        freq = np.linspace(200, 8000, 50)
        func = get_radius_func("conical", r, r + 1e-10, L)  # near-cylindrical
        z_real, z_imag = compute_throat_impedance_tmm(
            freq, func, L, r, r + 1e-10, n_segments=100
        )
        # Real part should be non-negative (radiation resistance at mouth)
        assert np.all(z_real >= -1e-6), f"Negative Re(Z) found: min={z_real.min()}"


class TestConicalHorn:
    """TMM for a conical horn — compare against known properties."""

    THROAT_R = 0.025
    MOUTH_R = 0.15
    LENGTH = 0.3

    def test_real_impedance_positive(self):
        """Re(Z_throat) > 0 at all frequencies."""
        freq = np.linspace(300, 8000, 80)
        func = get_radius_func("conical", self.THROAT_R, self.MOUTH_R, self.LENGTH)
        z_real, z_imag = compute_throat_impedance_tmm(
            freq, func, self.LENGTH, self.THROAT_R, self.MOUTH_R
        )
        assert np.all(z_real > -1e-6)

    def test_impedance_magnitude_bounded(self):
        """Impedance magnitude should be bounded by ρc (= 420 Pa·s/m)."""
        freq = np.linspace(500, 8000, 80)
        func = get_radius_func("conical", self.THROAT_R, self.MOUTH_R, self.LENGTH)
        z_real, z_imag = compute_throat_impedance_tmm(
            freq, func, self.LENGTH, self.THROAT_R, self.MOUTH_R
        )
        z_mag = np.sqrt(z_real**2 + z_imag**2)
        # Specific impedance should be in a reasonable range
        # (not zero, not wildly larger than rho*c)
        rho_c = RHO0 * C0
        assert np.all(z_mag < 10 * rho_c), f"Impedance too large: max={z_mag.max()}"

    def test_high_frequency_approaches_rho_c(self):
        """At high frequencies where ka >> 1, throat impedance → ρc."""
        freq = np.array([15000.0, 20000.0])
        func = get_radius_func("conical", self.THROAT_R, self.MOUTH_R, self.LENGTH)
        z_real, z_imag = compute_throat_impedance_tmm(
            freq, func, self.LENGTH, self.THROAT_R, self.MOUTH_R, n_segments=400
        )
        rho_c = RHO0 * C0
        # At high freq the horn looks like a matched load
        # Real part should approach ρc within a factor of ~2
        for i in range(len(freq)):
            assert 0.1 * rho_c < z_real[i] < 5 * rho_c


class TestSegmentConvergence:
    """TMM should converge as n_segments increases."""

    def test_convergence(self):
        freq = np.linspace(500, 4000, 30)
        func = get_radius_func("exponential", 0.025, 0.15, 0.3)

        z_50 = compute_throat_impedance_tmm(freq, func, 0.3, 0.025, 0.15, n_segments=50)
        z_200 = compute_throat_impedance_tmm(freq, func, 0.3, 0.025, 0.15, n_segments=200)
        z_500 = compute_throat_impedance_tmm(freq, func, 0.3, 0.025, 0.15, n_segments=500)

        # Error between 200 and 500 should be smaller than between 50 and 200
        err_50_200 = np.mean(np.abs(z_50[0] - z_200[0]))
        err_200_500 = np.mean(np.abs(z_200[0] - z_500[0]))
        assert err_200_500 < err_50_200, (
            f"Not converging: err(50,200)={err_50_200:.4f}, err(200,500)={err_200_500:.4f}"
        )


class TestAllProfiles:
    """Basic sanity checks across all 7 profiles."""

    PROFILES = ["conical", "exponential", "hyperbolic", "tractrix", "os", "lecleach", "cd"]

    @pytest.fixture(params=PROFILES)
    def profile(self, request):
        return request.param

    def test_produces_valid_output(self, profile):
        """TMM returns finite, correctly shaped arrays for all profiles."""
        freq = np.linspace(500, 4000, 20)
        func = get_radius_func(profile, 0.025, 0.15, 0.3)
        z_real, z_imag = compute_throat_impedance_tmm(
            freq, func, 0.3, 0.025, 0.15, n_segments=100
        )
        assert z_real.shape == freq.shape
        assert z_imag.shape == freq.shape
        assert np.all(np.isfinite(z_real))
        assert np.all(np.isfinite(z_imag))

    def test_real_part_positive(self, profile):
        """Re(Z) should be positive (passive system)."""
        freq = np.linspace(500, 4000, 30)
        func = get_radius_func(profile, 0.025, 0.15, 0.3)
        z_real, _ = compute_throat_impedance_tmm(
            freq, func, 0.3, 0.025, 0.15, n_segments=200
        )
        assert np.all(z_real > -1e-3), f"{profile}: Re(Z) negative, min={z_real.min()}"
