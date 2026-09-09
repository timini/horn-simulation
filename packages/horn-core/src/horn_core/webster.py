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
from horn_core.duct import DEFAULT_AIR, circular_duct_properties, circular_pipe_radiation

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


def compute_horn_transfer_tmm(
    frequencies: np.ndarray,
    radius_func: Callable[[float], float],
    length: float,
    throat_radius: float,
    mouth_radius: float,
    n_segments: int = 200,
    radiation_model: str = "flanged_piston",
    loss_model: str = "lossless",
    flange_width: float = 0.,
    air=DEFAULT_AIR,
) -> dict:
    """Compute unit-inlet-pressure impedance and mouth transfer using TMM.

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
        Dictionary of specific throat impedance (Pa·s/m), mouth pressure
        (Pa/Pa), mouth volume velocity (m³/s/Pa), and physical areas (m²).
    """
    from horn_core.acoustics import validate_response
    frequencies = validate_response(frequencies)
    if not np.isfinite([length, throat_radius, mouth_radius]).all() or min(length, throat_radius, mouth_radius) <= 0:
        raise ValueError("Horn dimensions must be finite and positive")
    if n_segments < 2:
        raise ValueError("At least two TMM segments are required")
    if radiation_model not in {"plane_wave", "flanged_piston", "modal_baffled", "unflanged", "unflanged_piston", "finite_flange", "closed"}:
        raise ValueError("TMM supports local radiation models only")
    z_mid = (np.arange(n_segments) + 0.5) * length / n_segments
    areas = np.pi * np.array([radius_func(z) for z in z_mid])**2
    if not np.isfinite(areas).all() or np.any(areas <= 0):
        raise ValueError("Invalid horn section area")
    k = 2 * np.pi * frequencies / air.c
    a_mouth = mouth_radius
    if radiation_model == "closed":
        z_norm = np.zeros(len(k), dtype=complex)
    elif radiation_model == "finite_flange":
        z_norm = circular_pipe_radiation(k, a_mouth, flange_width)
    elif radiation_model == "plane_wave":
        z_norm = np.ones(len(k), dtype=complex)
    elif radiation_model in {"flanged_piston", "modal_baffled"}:
        # Webster is a plane-mode prescreen. The subsequent modal FEM and
        # nonuniform aperture observer determine the final recommendation.
        z_norm = np.array([piston_radiation_impedance(ki, a_mouth) for ki in k])
    else:
        if np.any(k * a_mouth >= 1.5):
            raise ValueError("Unflanged low-ka approximation requires ka < 1.5")
        z_norm = np.array([unflanged_radiation_impedance(ki, a_mouth) for ki in k])
    mouth_area = np.pi * mouth_radius**2
    load = z_norm * air.rho * air.c / mouth_area
    matrices = np.tile(np.eye(2, dtype=complex), (len(k), 1, 1))
    for area in areas:
        propagation_k, specific_zc = circular_duct_properties(frequencies, np.sqrt(area/np.pi), air=air, loss_model=loss_model)
        cos_kl, sin_kl = np.cos(propagation_k*length/n_segments), np.sin(propagation_k*length/n_segments)
        zc = specific_zc / area
        segment = np.empty_like(matrices)
        segment[:, 0, 0] = segment[:, 1, 1] = cos_kl
        segment[:, 0, 1] = 1j * zc * sin_kl
        segment[:, 1, 0] = 1j * sin_kl / zc
        matrices = matrices @ segment
    pressure_denominator = matrices[:, 0, 0]*load + matrices[:, 0, 1]
    flow_denominator = matrices[:, 1, 0]*load + matrices[:, 1, 1]
    if radiation_model == "closed":
        pressure_denominator, flow_denominator = matrices[:, 0, 0], matrices[:, 1, 0]
    impedance = pressure_denominator / flow_denominator * np.pi*throat_radius**2
    return {
        "z_real": impedance.real, "z_imag": impedance.imag,
        "mouth_pressure_transfer": (1 if radiation_model == "closed" else load) / pressure_denominator,
        "mouth_volume_velocity_transfer": np.zeros(len(k), dtype=complex) if radiation_model == "closed" else 1 / pressure_denominator,
        "inlet_area_m2": np.pi*throat_radius**2,
        "mouth_area_m2": mouth_area,
    }


def compute_throat_impedance_tmm(*args, **kwargs):
    """Compatibility wrapper returning specific throat impedance components."""
    result = compute_horn_transfer_tmm(*args, **kwargs)
    return result["z_real"], result["z_imag"]
