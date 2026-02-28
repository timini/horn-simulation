"""Analytical radiation impedance models for circular apertures.

No dolfinx dependency — uses only numpy and scipy.special.
"""

import numpy as np
from scipy.special import j1, struve


def piston_radiation_impedance(k: float, a: float) -> complex:
    """Flanged circular piston radiation impedance (normalised by rho*c).

    Z_rad / (rho*c) = R1(2ka) + j * X1(2ka)

    where:
        R1(x) = 1 - 2*J1(x)/x
        X1(x) = 2*H1(x)/x    (H1 = Struve function of order 1)

    Parameters
    ----------
    k : wavenumber (rad/m)
    a : piston radius (m)

    Returns
    -------
    Complex specific impedance ratio Z_rad / (rho*c).
    """
    ka = k * a
    x = 2.0 * ka

    if x < 1e-12:
        # Taylor expansion for small x:
        # R1(x) ≈ x²/8,  X1(x) ≈ 4x/(3π)
        r1 = x**2 / 8.0
        x1 = 4.0 * x / (3.0 * np.pi)
    else:
        r1 = 1.0 - 2.0 * j1(x) / x
        x1 = 2.0 * struve(1, x) / x

    return complex(r1, x1)


def unflanged_radiation_impedance(k: float, a: float) -> complex:
    """Unflanged pipe radiation impedance (Levine-Schwinger approximation).

    Z_rad / (rho*c) ≈ (ka)²/4 + j * 0.6133 * ka

    Accurate for ka < 1.5 approximately.

    Parameters
    ----------
    k : wavenumber (rad/m)
    a : pipe radius (m)

    Returns
    -------
    Complex specific impedance ratio Z_rad / (rho*c).
    """
    ka = k * a
    return complex(0.25 * ka**2, 0.6133 * ka)
