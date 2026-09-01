"""Pure-Python radius functions for all 7 horn flare profiles.

Extracted from horn-geometry/generator.py without the gmsh dependency.
Only requires numpy. Each function maps axial position z ∈ [0, L] to
the local horn radius r(z).
"""

import math
from typing import Callable

import numpy as np


def conical_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Linear flare: r(z) = r_t + (r_m - r_t) * z / L."""
    return r_t + (r_m - r_t) * z / L


def exponential_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Exponential flare: r(z) = r_t * exp(m * z / L), m = ln(r_m / r_t)."""
    m = math.log(r_m / r_t)
    return r_t * math.exp(m * z / L)


def hyperbolic_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Hyperbolic (hypex) flare: r(z) = r_t * cosh(m * z / L), m = acosh(r_m / r_t)."""
    m = np.arccosh(r_m / r_t)
    return float(r_t * np.cosh(m * z / L))


def _build_tractrix_interp(r_t: float, r_m: float, L: float):
    """Pre-compute the tractrix interpolation arrays."""
    t = np.linspace(np.pi - 1e-6, np.pi / 2, 500)
    y, x = np.sin(t), np.log(np.tan(t / 2)) + np.cos(t)
    x -= x[0]
    x_n, y_n = x / x[-1], y / y[-1]
    return x_n * L, y_n


def tractrix_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Tractrix flare: rapid initial expansion, decelerating toward mouth."""
    x_scaled, y_n = _build_tractrix_interp(r_t, r_m, L)
    return r_t + (r_m - r_t) * float(np.interp(z, x_scaled, y_n))


def os_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Oblate spheroidal (OS / Geddes) flare: r(z) = sqrt(r_t² + (z·tanθ)²)."""
    theta = math.atan2(math.sqrt(r_m**2 - r_t**2), L)
    return math.sqrt(r_t**2 + (z * math.tan(theta)) ** 2)


def _build_lecleach_interp(r_t: float, r_m: float, L: float):
    """Pre-compute the Le Cléac'h interpolation arrays."""
    t = np.linspace(np.pi - 1e-6, np.pi / 2, 500)
    y, x = np.sin(t), np.log(np.tan(t / 2)) + np.cos(t)
    x -= x[0]
    idx = np.searchsorted(y, r_t / r_m)
    x_c, y_c = x[idx:] - x[idx], y[idx:]
    return x_c / x_c[-1] * L, y_c / y_c[-1] * r_m


def lecleach_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Le Cléac'h flare: tractrix with a = r_m, clipped at r >= r_t."""
    x_scaled, y_scaled = _build_lecleach_interp(r_t, r_m, L)
    return float(np.interp(z, x_scaled, y_scaled))


def cd_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Constant Directivity (CD) horn: exponential throat (30%) + conical body."""
    frac = 0.3
    z_t = frac * L
    r_trans = r_t * (r_m / r_t) ** frac

    if z <= z_t:
        return r_t * math.exp(math.log(r_trans / r_t) * z / z_t)
    else:
        return r_trans + (r_m - r_trans) * (z - z_t) / (L - z_t)


_PROFILE_DISPATCH = {
    "conical": conical_radius,
    "exponential": exponential_radius,
    "hyperbolic": hyperbolic_radius,
    "tractrix": tractrix_radius,
    "os": os_radius,
    "lecleach": lecleach_radius,
    "cd": cd_radius,
}


def get_radius_func(
    profile: str, r_t: float, r_m: float, L: float
) -> Callable[[float], float]:
    """Return a closure r(z) for the given profile and geometry.

    Args:
        profile: One of "conical", "exponential", "hyperbolic", "tractrix",
                 "os", "lecleach", "cd".
        r_t: Throat radius (m).
        r_m: Mouth radius (m).
        L: Horn length (m).

    Returns:
        Callable that maps z ∈ [0, L] → radius (m).
    """
    if profile not in _PROFILE_DISPATCH:
        raise ValueError(
            f"Unknown profile '{profile}'. Choose from: {list(_PROFILE_DISPATCH.keys())}"
        )

    raw_func = _PROFILE_DISPATCH[profile]

    # For profiles with expensive interpolation setup, pre-compute once
    if profile == "tractrix":
        x_scaled, y_n = _build_tractrix_interp(r_t, r_m, L)

        def _tractrix(z: float) -> float:
            return r_t + (r_m - r_t) * float(np.interp(z, x_scaled, y_n))

        return _tractrix

    if profile == "lecleach":
        x_scaled, y_scaled = _build_lecleach_interp(r_t, r_m, L)

        def _lecleach(z: float) -> float:
            return float(np.interp(z, x_scaled, y_scaled))

        return _lecleach

    def _radius(z: float) -> float:
        return raw_func(z, r_t, r_m, L)

    return _radius
