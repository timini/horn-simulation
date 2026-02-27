"""V5: Hyperbolic (cosh) horn — Webster equation analytical validation.

Physics
-------
The Webster horn equation for a cosh horn has an analytical solution:

    p(x) = [A * exp(-j*k_h*x) + B * exp(j*k_h*x)] / cosh(beta*x)

where beta = acosh(r_mouth/r_throat)/L and
k_h = sqrt(k² - beta²) is the effective wavenumber (complex below cutoff).

The cutoff frequency is fc = beta*c0/(2*pi). Below cutoff, waves are evanescent.

BCs: Dirichlet p=1 at throat, Robin dp/dx = -jkp at mouth.
The Robin BC includes a tanh(beta*L) term from the cosh envelope derivative.

Geometry: throat=12.5mm, mouth=100mm, length=300mm (8:1 radius ratio,
matching V2 but with hyperbolic flare to isolate profile shape effects).

Tolerance: 3 dB (Webster is 1D approximation of 3D).
"""

import numpy as np
import pytest

from .conftest import assert_spl_within_tolerance


def webster_hyperbolic_spl(
    frequencies: np.ndarray,
    throat_radius: float,
    mouth_radius: float,
    length: float,
    c0: float = 343.0,
    p_ref: float = 20e-6,
) -> np.ndarray:
    """Compute SPL at the mouth of a hyperbolic (cosh) horn using the Webster equation.

    For a cosh horn, the area function is S(x) = S_t * cosh²(beta*x),
    giving radius r(x) = r_t * cosh(beta*x) where beta = acosh(r_m/r_t)/L.

    The general solution is:
        p(x) = [A * exp(-j*k_h*x) + B * exp(j*k_h*x)] / cosh(beta*x)

    where k_h = sqrt(k² - beta²).

    Boundary conditions:
        - Throat (x = 0): p = 1 (Dirichlet, unit pressure)
        - Mouth  (x = L): dp/dx = -jkp (radiation / Sommerfeld)

    Parameters
    ----------
    frequencies : frequency array (Hz)
    throat_radius, mouth_radius, length : geometry (m)
    c0 : speed of sound (m/s)
    p_ref : reference pressure (Pa)

    Returns
    -------
    spl : SPL at the mouth for each frequency (dB)
    """
    beta = np.arccosh(mouth_radius / throat_radius) / length

    spl = np.zeros_like(frequencies)
    for i, freq in enumerate(frequencies):
        k = 2 * np.pi * freq / c0

        # Effective wavenumber (complex below cutoff)
        k_h = np.sqrt(k**2 - beta**2 + 0j)

        # p(x) = [A*exp(-j*k_h*x) + B*exp(j*k_h*x)] / cosh(beta*x)
        #
        # dp/dx = [A*(-j*k_h)*exp(-j*k_h*x) + B*(j*k_h)*exp(j*k_h*x)] / cosh(beta*x)
        #       - beta*tanh(beta*x) * [A*exp(-j*k_h*x) + B*exp(j*k_h*x)] / cosh(beta*x)
        #
        # At throat (x=0): cosh(0)=1, so p(0) = A + B = 1
        #
        # Robin BC at x=L: dp/dx(L) = -j*k*p(L)
        # [A*(-j*k_h)*e_neg + B*(j*k_h)*e_pos] / cosh(beta*L)
        #   - beta*tanh(beta*L) * [A*e_neg + B*e_pos] / cosh(beta*L)
        #   = -j*k * [A*e_neg + B*e_pos] / cosh(beta*L)
        #
        # Multiply through by cosh(beta*L):
        # A*(-j*k_h - beta*tanh(beta*L))*e_neg + B*(j*k_h - beta*tanh(beta*L))*e_pos
        #   = -j*k * [A*e_neg + B*e_pos]

        e_neg = np.exp(-1j * k_h * length)
        e_pos = np.exp(1j * k_h * length)
        tanh_bL = np.tanh(beta * length)

        # Row 1 (throat): A + B = 1
        # Row 2 (mouth Robin, rearranged):
        #   A*(-j*k_h - beta*tanh_bL + j*k)*e_neg + B*(j*k_h - beta*tanh_bL + j*k)*e_pos = 0
        c1 = (-1j * k_h - beta * tanh_bL + 1j * k)
        c2 = (1j * k_h - beta * tanh_bL + 1j * k)

        M = np.array([
            [1.0, 1.0],
            [c1 * e_neg, c2 * e_pos],
        ])
        rhs = np.array([1.0, 0.0])

        AB = np.linalg.solve(M, rhs)
        A, B = AB[0], AB[1]

        # Pressure at mouth (x=L)
        cosh_bL = np.cosh(beta * length)
        p_mouth = (A * e_neg + B * e_pos) / cosh_bL
        p_rms = np.abs(p_mouth)

        spl[i] = 20 * np.log10(p_rms / p_ref + 1e-12)

    return spl


@pytest.mark.validation
class TestHyperbolicHornWebster:
    """V5 validation: hyperbolic horn FEM vs Webster equation."""

    def test_spl_within_webster_tolerance(self, hyperbolic_horn_results):
        """FEM SPL should be within 3 dB of Webster prediction at each frequency."""
        frequencies, fem_spl, ref = hyperbolic_horn_results
        geom = ref["geometry"]
        tolerance = ref["expected"]["tolerance_db"]

        webster_spl = webster_hyperbolic_spl(
            frequencies=frequencies,
            throat_radius=geom["throat_radius_m"],
            mouth_radius=geom["mouth_radius_m"],
            length=geom["length_m"],
        )

        assert_spl_within_tolerance(
            computed=fem_spl,
            reference=webster_spl,
            tolerance_db=tolerance,
            frequencies=frequencies,
            label="V5 hyperbolic horn (FEM vs Webster)",
        )

    def test_highpass_behavior(self, hyperbolic_horn_results):
        """SPL should be higher above cutoff than below (high-pass characteristic)."""
        frequencies, spl, ref = hyperbolic_horn_results
        cutoff = ref["expected"]["cutoff_hz"]

        below = spl[frequencies < cutoff]
        above = spl[frequencies > cutoff]

        if len(below) > 0 and len(above) > 0:
            assert np.mean(above) > np.mean(below), (
                f"Hyperbolic horn should show high-pass behavior around {cutoff} Hz.\n"
                f"Mean SPL below cutoff: {np.mean(below):.1f} dB\n"
                f"Mean SPL above cutoff: {np.mean(above):.1f} dB"
            )

    def test_webster_function_sanity(self):
        """Sanity check: Webster hyperbolic solution gives reasonable SPL."""
        frequencies = np.array([1000.0])
        spl = webster_hyperbolic_spl(
            frequencies=frequencies,
            throat_radius=0.0125,
            mouth_radius=0.1,
            length=0.3,
        )
        assert np.isfinite(spl[0]), f"Webster SPL is not finite: {spl[0]}"
        assert 70 < spl[0] < 120, f"Webster SPL outside reasonable range: {spl[0]:.1f} dB"
