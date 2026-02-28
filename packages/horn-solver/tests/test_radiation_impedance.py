"""Unit tests for radiation impedance models (no dolfinx needed)."""

import numpy as np
import pytest
from scipy.special import j1, struve

from horn_solver.radiation import (
    piston_radiation_impedance,
    unflanged_radiation_impedance,
)


class TestPistonRadiationImpedance:
    """Tests for the flanged circular piston model."""

    def test_dc_limit(self):
        """At k=0, radiation impedance should be zero."""
        z = piston_radiation_impedance(k=0.0, a=0.05)
        assert z == 0.0 + 0.0j

    def test_small_argument_taylor(self):
        """For ka << 1, verify against Taylor expansion.

        R1(x) ≈ x²/8,  X1(x) ≈ 4x/(3π)  where x = 2ka.
        """
        k = 1.0  # small k
        a = 0.001  # small a => ka = 0.001
        x = 2.0 * k * a

        expected_real = x**2 / 8.0
        expected_imag = 4.0 * x / (3.0 * np.pi)

        z = piston_radiation_impedance(k, a)

        assert abs(z.real - expected_real) < 1e-8
        assert abs(z.imag - expected_imag) < 1e-8

    def test_high_frequency_limit(self):
        """For ka >> 1, Z should approach 1 + 0j (plane wave limit).

        At large x, J1(x)/x -> 0 and H1(x)/x -> 2/(πx) -> 0.
        """
        k = 1000.0
        a = 0.1  # ka = 100

        z = piston_radiation_impedance(k, a)

        assert abs(z.real - 1.0) < 0.05, f"Real part should approach 1.0, got {z.real}"
        assert abs(z.imag) < 0.05, f"Imag part should approach 0, got {z.imag}"

    def test_known_tabulated_value(self):
        """Cross-check against a known value at ka = 1.

        At x = 2ka = 2:
            R1(2) = 1 - 2*J1(2)/2 = 1 - J1(2) ≈ 1 - 0.5767 = 0.4233
            X1(2) = 2*H1(2)/2 = H1(2) ≈ 0.6468
        """
        k = 10.0
        a = 0.1  # ka = 1, x = 2

        expected_real = 1.0 - j1(2.0)
        expected_imag = struve(1, 2.0)

        z = piston_radiation_impedance(k, a)

        assert abs(z.real - expected_real) < 1e-10
        assert abs(z.imag - expected_imag) < 1e-10

    def test_impedance_is_passive(self):
        """Real part should always be non-negative (passive system)."""
        for ka in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
            k = ka / 0.05
            z = piston_radiation_impedance(k, 0.05)
            assert z.real >= 0.0, f"Real part negative at ka={ka}: {z.real}"

    def test_monotonic_real_part_trend(self):
        """Real part should generally increase from 0 towards 1 with ka."""
        a = 0.05
        ka_values = [0.01, 0.1, 0.5, 1.0, 5.0, 20.0]
        reals = [piston_radiation_impedance(ka / a, a).real for ka in ka_values]

        # Not strictly monotonic (oscillations around 1), but should increase
        # from near 0 to near 1
        assert reals[0] < 0.01
        assert reals[-1] > 0.9


class TestUnflangedRadiationImpedance:
    """Tests for the Levine-Schwinger unflanged pipe model."""

    def test_dc_limit(self):
        """At k=0, impedance should be zero."""
        z = unflanged_radiation_impedance(k=0.0, a=0.05)
        assert z == 0.0 + 0.0j

    def test_formula_verification(self):
        """Verify the formula Z/(rho*c) = (ka)²/4 + j*0.6133*ka."""
        k = 20.0
        a = 0.05
        ka = k * a  # = 1.0

        z = unflanged_radiation_impedance(k, a)

        expected_real = 0.25 * ka**2
        expected_imag = 0.6133 * ka

        assert abs(z.real - expected_real) < 1e-12
        assert abs(z.imag - expected_imag) < 1e-12

    def test_real_part_quadratic(self):
        """Real part should scale as (ka)²."""
        a = 0.05
        z1 = unflanged_radiation_impedance(10.0, a)  # ka = 0.5
        z2 = unflanged_radiation_impedance(20.0, a)  # ka = 1.0

        # ratio of real parts should be (1.0/0.5)² = 4
        ratio = z2.real / z1.real
        assert abs(ratio - 4.0) < 1e-10

    def test_imag_part_linear(self):
        """Imaginary part should scale linearly with ka."""
        a = 0.05
        z1 = unflanged_radiation_impedance(10.0, a)  # ka = 0.5
        z2 = unflanged_radiation_impedance(20.0, a)  # ka = 1.0

        ratio = z2.imag / z1.imag
        assert abs(ratio - 2.0) < 1e-10

    def test_impedance_is_passive(self):
        """Real part should always be non-negative."""
        for ka in [0.01, 0.1, 0.5, 1.0, 1.5]:
            k = ka / 0.05
            z = unflanged_radiation_impedance(k, 0.05)
            assert z.real >= 0.0


class TestModelComparison:
    """Sanity checks comparing the two models."""

    def test_unflanged_less_than_flanged_at_low_freq(self):
        """Unflanged pipe radiates less efficiently than flanged piston."""
        k = 20.0
        a = 0.05  # ka = 1.0

        z_flanged = piston_radiation_impedance(k, a)
        z_unflanged = unflanged_radiation_impedance(k, a)

        assert z_unflanged.real < z_flanged.real, (
            f"Unflanged real ({z_unflanged.real:.4f}) should be less than "
            f"flanged real ({z_flanged.real:.4f}) at ka=1"
        )
