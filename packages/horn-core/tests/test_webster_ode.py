"""Check the Webster TMM against an independent ODE solution of the horn equation.

The reference integrates the lossless Webster equation written as a first-order
system in pressure p and volume velocity U,

    dp/dz = -j w rho U / S(z),    dU/dz = -j w S(z) p / (rho c^2),

with scipy's adaptive DOP853 integrator, from the mouth (p = Z_L U) back to the
throat. It shares only the profile function and the mouth radiation load with
the TMM, so it checks the segment cascade and its throat/mouth bookkeeping.

The tolerances are frozen from the measured error of the current TMM
(midpoint stepped cascade, second order in the segment length). Do not widen
them to make a change pass; a larger error is a regression.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from horn_core.duct import DEFAULT_AIR as AIR
from horn_core.profiles import get_radius_func
from horn_core.webster import compute_horn_transfer_tmm, piston_radiation_impedance

THROAT_RADIUS, MOUTH_RADIUS, LENGTH = 0.0125, 0.15, 0.4
FREQUENCIES = np.geomspace(100.0, 8000.0, 40)
# Measured max relative throat impedance error at n=200: exponential 0.22 %,
# conical 0.56 %. At n=2000: 2.2e-5 and 5.6e-5.
DEFAULT_SEGMENT_TOLERANCE = 0.01
FINE_SEGMENT_TOLERANCE = 1e-4


def webster_ode(radius_func):
    """Return specific throat impedance and mouth/throat pressure ratio."""
    impedance, transfer = [], []
    mouth_area = np.pi * MOUTH_RADIUS**2
    for frequency in FREQUENCIES:
        omega = 2 * np.pi * frequency
        load = piston_radiation_impedance(omega / AIR.c, MOUTH_RADIUS) * AIR.rho * AIR.c / mouth_area

        def rhs(z, y):
            p, u = y[0] + 1j * y[1], y[2] + 1j * y[3]
            area = np.pi * radius_func(z) ** 2
            dp = -1j * omega * AIR.rho * u / area
            du = -1j * omega * area * p / (AIR.rho * AIR.c**2)
            return [dp.real, dp.imag, du.real, du.imag]

        solution = solve_ivp(rhs, (LENGTH, 0.0), [load.real, load.imag, 1.0, 0.0],
                             method="DOP853", rtol=1e-10, atol=1e-14)
        assert solution.success, solution.message
        p0 = solution.y[0, -1] + 1j * solution.y[1, -1]
        u0 = solution.y[2, -1] + 1j * solution.y[3, -1]
        impedance.append(p0 / u0 * np.pi * THROAT_RADIUS**2)
        transfer.append(load / p0)
    return np.array(impedance), np.array(transfer)


def tmm(radius_func, n_segments):
    result = compute_horn_transfer_tmm(FREQUENCIES, radius_func, LENGTH, THROAT_RADIUS,
                                       MOUTH_RADIUS, n_segments=n_segments)
    return result["z_real"] + 1j * result["z_imag"], result["mouth_pressure_transfer"]


@pytest.fixture(scope="module", params=["exponential", "conical"])
def case(request):
    radius_func = get_radius_func(request.param, THROAT_RADIUS, MOUTH_RADIUS, LENGTH)
    return radius_func, webster_ode(radius_func)


def relative_error(value, reference):
    return np.max(np.abs(value - reference) / np.abs(reference))


@pytest.mark.parametrize("n_segments,tolerance", [(200, DEFAULT_SEGMENT_TOLERANCE),
                                                  (2000, FINE_SEGMENT_TOLERANCE)])
def test_tmm_matches_webster_ode(case, n_segments, tolerance):
    radius_func, (reference_z, reference_transfer) = case
    impedance, transfer = tmm(radius_func, n_segments)
    assert relative_error(impedance, reference_z) < tolerance
    assert relative_error(transfer, reference_transfer) < tolerance


def test_tmm_converges_at_second_order(case):
    radius_func, (reference_z, _) = case
    coarse = relative_error(tmm(radius_func, 200)[0], reference_z)
    fine = relative_error(tmm(radius_func, 400)[0], reference_z)
    assert 3.5 < coarse / fine < 4.5
