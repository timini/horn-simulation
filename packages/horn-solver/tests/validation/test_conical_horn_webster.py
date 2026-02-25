"""V2: Conical horn — Webster equation analytical validation.

Physics
-------
The Webster horn equation for a conical frustum has an analytical solution
in terms of spherical waves:

    p(x) = (x_t / x) * [A * exp(-jkx) + B * exp(jkx)]

where x is distance from the virtual apex. A and B are set by boundary
conditions: Dirichlet p=1 at the throat and Robin dp/dn = -jkp at the mouth.

The FEM solves the full 3D Helmholtz equation; Webster is a 1D plane-wave
approximation. Agreement should be qualitative (within ~3 dB), not exact.

Parameters: throat_radius=25mm, mouth_radius=100mm, length=500mm (4:1 expansion).

Tolerance: 3 dB (Webster is 1D approximation of 3D).
"""

import numpy as np
import pytest

from .conftest import assert_spl_within_tolerance


def webster_conical_spl(
    frequencies: np.ndarray,
    throat_radius: float,
    mouth_radius: float,
    length: float,
    c0: float = 343.0,
    p_ref: float = 20e-6,
) -> np.ndarray:
    """Compute SPL at the mouth of a conical horn using the Webster equation.

    For a conical frustum, the area function is S(z) = S_t * ((x_t + z) / x_t)^2,
    where x_t = L * r_t / (r_m - r_t) is the distance from virtual apex to throat.

    The general solution in conical coordinates is:
        p(x) = (1/x) * [A * exp(-jkx) + B * exp(jkx)]

    Boundary conditions:
        - Throat (x = x_t): p = 1 (Dirichlet, unit pressure)
        - Mouth  (x = x_m): dp/dx = -jkp (radiation / Sommerfeld)

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
    # Distance from virtual apex to throat and mouth
    x_t = length * throat_radius / (mouth_radius - throat_radius)
    x_m = x_t + length

    spl = np.zeros_like(frequencies)
    for i, freq in enumerate(frequencies):
        k = 2 * np.pi * freq / c0

        e_t = np.exp(-1j * k * x_t)
        f_t = np.exp(1j * k * x_t)
        e_m = np.exp(-1j * k * x_m)
        f_m = np.exp(1j * k * x_m)

        # 2x2 system from throat Dirichlet + mouth Robin BCs:
        #   Row 1 (throat): A*e_t + B*f_t = x_t
        #   Row 2 (mouth):  -A*e_m + B*f_m*(-1 + 2jk*x_m) = 0
        M = np.array([
            [e_t, f_t],
            [-e_m, f_m * (-1 + 2j * k * x_m)],
        ])
        rhs = np.array([x_t, 0.0])

        AB = np.linalg.solve(M, rhs)
        A, B = AB[0], AB[1]

        # Pressure at mouth
        p_mouth = (1.0 / x_m) * (A * e_m + B * f_m)
        p_rms = np.abs(p_mouth)

        spl[i] = 20 * np.log10(p_rms / p_ref + 1e-12)

    return spl


@pytest.mark.validation
class TestConicalHornWebster:
    """V2 validation: conical horn FEM vs Webster equation."""

    def test_spl_within_webster_tolerance(self, conical_horn_results):
        """FEM SPL should be within 3 dB of Webster prediction at each frequency."""
        frequencies, fem_spl, ref = conical_horn_results
        geom = ref["geometry"]
        tolerance = ref["expected"]["tolerance_db"]

        # Compute Webster analytical reference
        webster_spl = webster_conical_spl(
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
            label="V2 conical horn (FEM vs Webster)",
        )

    def test_spl_increases_with_frequency(self, conical_horn_results):
        """SPL should generally increase with frequency in the horn passband."""
        frequencies, spl, ref = conical_horn_results

        # Compare first quarter and last quarter averages
        n = len(spl)
        quarter = n // 4
        low_freq_avg = np.mean(spl[:quarter])
        high_freq_avg = np.mean(spl[-quarter:])

        assert high_freq_avg > low_freq_avg, (
            f"SPL should increase with frequency in conical horn passband.\n"
            f"Low-freq avg ({frequencies[0]:.0f}-{frequencies[quarter-1]:.0f} Hz): "
            f"{low_freq_avg:.1f} dB\n"
            f"High-freq avg ({frequencies[-quarter]:.0f}-{frequencies[-1]:.0f} Hz): "
            f"{high_freq_avg:.1f} dB"
        )

    def test_webster_function_sanity(self):
        """Sanity check: Webster solution for a mild cone gives reasonable SPL."""
        frequencies = np.array([1000.0])
        spl = webster_conical_spl(
            frequencies=frequencies,
            throat_radius=0.05,
            mouth_radius=0.1,
            length=0.5,
        )
        # Should produce a finite, reasonable SPL value
        assert np.isfinite(spl[0]), f"Webster SPL is not finite: {spl[0]}"
        assert 70 < spl[0] < 120, f"Webster SPL outside reasonable range: {spl[0]:.1f} dB"
