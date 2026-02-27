"""V4: Exponential horn — Webster equation analytical validation.

Physics
-------
The Webster horn equation for an exponential horn has an analytical solution:

    p(x) = exp(-m*x) * [A * exp(-j*k_e*x) + B * exp(j*k_e*x)]

where m = ln(r_mouth/r_throat)/L is the flare rate and
k_e = sqrt(k² - m²) is the effective wavenumber (complex below cutoff).

The cutoff frequency is fc = m*c0/(2*pi). Below cutoff, waves are evanescent.

BCs: Dirichlet p=1 at throat, Robin dp/dx = -jkp at mouth.
This gives a 2x2 system for coefficients A and B.

Geometry: Rasetshwane & Neely (JASA 2012, PMC3316681) — axisymmetric
exponential horn. throat=5mm, mouth=22.95mm, length=120mm.

Tolerance: 3 dB (Webster is 1D approximation of 3D).
"""

import numpy as np
import pytest

from .conftest import assert_spl_within_tolerance


def webster_exponential_spl(
    frequencies: np.ndarray,
    throat_radius: float,
    mouth_radius: float,
    length: float,
    c0: float = 343.0,
    p_ref: float = 20e-6,
) -> np.ndarray:
    """Compute SPL at the mouth of an exponential horn using the Webster equation.

    For an exponential horn, the area function is S(x) = S_t * exp(2*m*x/L),
    giving radius r(x) = r_t * exp(m*x/L) where m = ln(r_m/r_t).

    The general solution is:
        p(x) = exp(-m'*x) * [A * exp(-j*k_e*x) + B * exp(j*k_e*x)]

    where m' = m/L is the flare rate per unit length (applied to area: 2m'/1),
    and k_e = sqrt(k² - m'²).

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
    m = np.log(mouth_radius / throat_radius)  # ln(r_m/r_t)
    alpha = m / length  # flare rate (1/m), applied to radius; area flare = 2*alpha

    spl = np.zeros_like(frequencies)
    for i, freq in enumerate(frequencies):
        k = 2 * np.pi * freq / c0

        # Effective wavenumber (complex below cutoff)
        k_e = np.sqrt(k**2 - alpha**2 + 0j)

        # At throat (x=0): p(0) = A + B = 1
        # At mouth (x=L): dp/dx = -jkp
        #
        # p(x) = exp(-alpha*x) * [A*exp(-j*k_e*x) + B*exp(j*k_e*x)]
        # dp/dx = exp(-alpha*x) * [A*(-alpha - j*k_e)*exp(-j*k_e*x)
        #                         + B*(-alpha + j*k_e)*exp(j*k_e*x)]
        #
        # Robin BC at x=L: dp/dx(L) = -j*k*p(L)
        # A*(-alpha - j*k_e)*exp(-j*k_e*L) + B*(-alpha + j*k_e)*exp(j*k_e*L)
        #   = -j*k * [A*exp(-j*k_e*L) + B*exp(j*k_e*L)]

        e_neg = np.exp(-1j * k_e * length)
        e_pos = np.exp(1j * k_e * length)

        # Row 1 (throat): A + B = 1
        # Row 2 (mouth Robin): A*c1*e_neg + B*c2*e_pos = 0
        #   where c1 = (-alpha - j*k_e) + j*k = j*(k - k_e) - alpha
        #         c2 = (-alpha + j*k_e) + j*k = j*(k + k_e) - alpha
        c1 = -alpha - 1j * k_e + 1j * k
        c2 = -alpha + 1j * k_e + 1j * k

        M = np.array([
            [1.0, 1.0],
            [c1 * e_neg, c2 * e_pos],
        ])
        rhs = np.array([1.0, 0.0])

        AB = np.linalg.solve(M, rhs)
        A, B = AB[0], AB[1]

        # Pressure at mouth (x=L)
        p_mouth = np.exp(-alpha * length) * (A * e_neg + B * e_pos)
        p_rms = np.abs(p_mouth)

        spl[i] = 20 * np.log10(p_rms / p_ref + 1e-12)

    return spl


@pytest.mark.validation
class TestExponentialHornWebster:
    """V4 validation: exponential horn FEM vs Webster equation."""

    def test_spl_within_webster_tolerance(self, exponential_horn_results):
        """FEM SPL should be within 3 dB of Webster prediction at each frequency."""
        frequencies, fem_spl, ref = exponential_horn_results
        geom = ref["geometry"]
        tolerance = ref["expected"]["tolerance_db"]

        webster_spl = webster_exponential_spl(
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
            label="V4 exponential horn (FEM vs Webster)",
        )

    def test_highpass_behavior(self, exponential_horn_results):
        """SPL should be higher above cutoff than below (high-pass characteristic)."""
        frequencies, spl, ref = exponential_horn_results
        cutoff = ref["expected"]["cutoff_hz"]

        below = spl[frequencies < cutoff]
        above = spl[frequencies > cutoff]

        if len(below) > 0 and len(above) > 0:
            assert np.mean(above) > np.mean(below), (
                f"Exponential horn should show high-pass behavior around {cutoff} Hz.\n"
                f"Mean SPL below cutoff: {np.mean(below):.1f} dB\n"
                f"Mean SPL above cutoff: {np.mean(above):.1f} dB"
            )

    def test_webster_function_sanity(self):
        """Sanity check: Webster exponential solution gives reasonable SPL."""
        frequencies = np.array([1000.0])
        spl = webster_exponential_spl(
            frequencies=frequencies,
            throat_radius=0.005,
            mouth_radius=0.02295,
            length=0.12,
        )
        assert np.isfinite(spl[0]), f"Webster SPL is not finite: {spl[0]}"
        assert 70 < spl[0] < 120, f"Webster SPL outside reasonable range: {spl[0]:.1f} dB"
