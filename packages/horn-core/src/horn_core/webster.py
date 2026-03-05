"""Transfer Matrix Method (TMM) solver for horn throat impedance.

Slices the horn into N thin cylindrical segments, builds a 2×2 transfer
matrix per segment using plane-wave propagation at the local cross-section
area, and cascades from throat → mouth. A radiation impedance boundary
condition is applied at the mouth.

With fine segmentation (N ≥ 200) this converges to the Webster horn
equation solution, giving millisecond-speed approximations of the FEM
solver output suitable for prescreening large candidate grids.
"""

import numpy as np
from scipy.special import j1, struve
from typing import Callable, Tuple

# Speed of sound and air density at ~20 °C
C0 = 343.0
RHO0 = 1.225


def piston_radiation_impedance(k: float, a: float) -> complex:
    """Flanged circular piston radiation impedance normalised by ρc.

    Z_rad / (ρc) = R1(2ka) + j·X1(2ka)
    where R1(x) = 1 - 2·J1(x)/x,  X1(x) = 2·H1(x)/x.

    Duplicated from horn_solver/radiation.py to avoid dolfinx dependency.
    """
    ka = k * a
    x = 2.0 * ka

    if x < 1e-12:
        r1 = x**2 / 8.0
        x1 = 4.0 * x / (3.0 * np.pi)
    else:
        r1 = 1.0 - 2.0 * j1(x) / x
        x1 = 2.0 * struve(1, x) / x

    return complex(r1, x1)


def unflanged_radiation_impedance(k: float, a: float) -> complex:
    """Unflanged pipe radiation impedance (Levine-Schwinger), normalised by ρc."""
    ka = k * a
    return complex(0.25 * ka**2, 0.6133 * ka)


def _segment_matrix(k: float, S: float, dz: float) -> np.ndarray:
    """2×2 transfer matrix for a cylindrical segment (plane-wave).

    Convention: [p_in, U_in] = T · [p_out, U_out]
    where p is pressure and U is volume velocity.

    Args:
        k: Wavenumber (rad/m).
        S: Cross-section area (m²).
        dz: Segment length (m).
    """
    cos_kl = np.cos(k * dz)
    sin_kl = np.sin(k * dz)
    Zc = RHO0 * C0 / S  # characteristic acoustic impedance

    return np.array([
        [cos_kl, 1j * Zc * sin_kl],
        [1j * sin_kl / Zc, cos_kl],
    ], dtype=complex)


def compute_throat_impedance_tmm(
    frequencies: np.ndarray,
    radius_func: Callable[[float], float],
    length: float,
    throat_radius: float,
    mouth_radius: float,
    n_segments: int = 200,
    radiation_model: str = "flanged_piston",
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute specific acoustic impedance at the horn throat using TMM.

    Slices the horn into cylindrical segments at their midpoint radius,
    builds the cascaded transfer matrix from throat to mouth, and applies
    a radiation impedance boundary condition at the mouth.

    Args:
        frequencies: Array of frequencies in Hz.
        radius_func: Callable mapping z ∈ [0, length] to radius (m).
        length: Horn length (m).
        throat_radius: Throat radius (m).
        mouth_radius: Mouth radius (m).
        n_segments: Number of segments (≥50 recommended).
        radiation_model: "flanged_piston" or "unflanged".

    Returns:
        Tuple of (z_real, z_imag) arrays — specific acoustic impedance
        at the throat in Pa·s/m, matching the FEM solver CSV format.
    """
    # Segment boundaries and midpoint radii
    z_edges = np.linspace(0, length, n_segments + 1)
    dz = length / n_segments
    z_mid = 0.5 * (z_edges[:-1] + z_edges[1:])
    r_mid = np.array([radius_func(z) for z in z_mid])
    S_mid = np.pi * r_mid**2

    z_real_out = np.zeros(len(frequencies))
    z_imag_out = np.zeros(len(frequencies))

    a_mouth = radius_func(length)
    S_throat = np.pi * throat_radius**2

    for i_f, freq in enumerate(frequencies):
        k = 2.0 * np.pi * freq / C0

        # Radiation impedance at the mouth (normalised → acoustic)
        if radiation_model == "flanged_piston":
            z_rad_norm = piston_radiation_impedance(k, a_mouth)
        else:
            z_rad_norm = unflanged_radiation_impedance(k, a_mouth)

        S_mouth = np.pi * a_mouth**2
        Z_load = z_rad_norm * RHO0 * C0 / S_mouth  # acoustic impedance

        # Cascade transfer matrices: T_total = T_0 · T_1 · ... · T_{N-1}
        # Maps (p, U) at throat to (p, U) at mouth
        T = np.eye(2, dtype=complex)
        for seg in range(n_segments):
            T_seg = _segment_matrix(k, S_mid[seg], dz)
            T = T @ T_seg

        # Throat acoustic impedance from transfer matrix:
        # p_throat = T11*p_mouth + T12*U_mouth
        # U_throat = T21*p_mouth + T22*U_mouth
        # At mouth: p_mouth = Z_load * U_mouth
        # So: Z_throat = p_throat/U_throat
        #   = (T11*Z_load + T12) / (T21*Z_load + T22)
        Z_throat_acoustic = (T[0, 0] * Z_load + T[0, 1]) / (T[1, 0] * Z_load + T[1, 1])

        # Convert acoustic → specific acoustic impedance
        Z_throat_specific = Z_throat_acoustic * S_throat

        z_real_out[i_f] = np.real(Z_throat_specific)
        z_imag_out[i_f] = np.imag(Z_throat_specific)

    return z_real_out, z_imag_out
