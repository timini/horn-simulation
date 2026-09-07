"""Shared acoustic conventions: RMS phasors with time dependence exp(+i wt).

Pressure is Pa, particle velocity m/s, volume velocity m³/s, specific
impedance Pa s/m. Observation distances are measured from the mouth plane.
"""
import numpy as np

C0 = 343.0
RHO0 = 1.225
P_REF = 20e-6
SCHEMA_VERSION = 2


def validate_band(low, high):
    if not np.isfinite([low, high]).all() or not 0 < low < high:
        raise ValueError("Frequency band must be finite and satisfy 0 < low < high")


def validate_response(frequencies, *arrays):
    f = np.asarray(frequencies, dtype=float)
    if f.ndim != 1 or len(f) < 2 or not np.isfinite(f).all():
        raise ValueError("Response needs at least two finite frequency samples")
    if np.any(f <= 0) or np.any(np.diff(f) <= 0):
        raise ValueError("Frequencies must be positive, unique and increasing")
    for array in arrays:
        a = np.asarray(array)
        if a.shape != f.shape or not np.isfinite(a).all():
            raise ValueError("Response arrays must match frequencies and be finite")
    return f


def pressure_level(pressure):
    return 20 * np.log10(np.maximum(np.abs(pressure), 1e-30) / P_REF)


def baffled_piston_on_axis(frequencies, volume_velocity, area_m2, distance_m=1.0):
    """Exact on-axis Rayleigh result for a UNIFORM circular baffled piston.

    Applying this to an arbitrary horn aperture is an approximation, not a
    full exterior solve. No claim is made about an unbaffled enclosure.
    Formula: rho*c*v*(exp(-ikr)-exp(-ik*sqrt(r²+a²))).
    """
    if not np.isfinite([area_m2, distance_m]).all() or min(area_m2, distance_m) <= 0:
        raise ValueError("Area and observation distance must be finite and positive")
    f = np.asarray(frequencies)
    u = np.asarray(volume_velocity)
    k = 2 * np.pi * f / C0
    a2 = area_m2 / np.pi
    # Rationalised difference and expm1 avoid cancellation for small apertures.
    delta = a2 / (np.sqrt(distance_m**2 + a2) + distance_m)
    return RHO0 * C0 * u / area_m2 * np.exp(-1j*k*distance_m) * (-np.expm1(-1j*k*delta))


def inlet_area_from_frame(frame, fallback_radius=None):
    """Use tagged physical boundary area; legacy CSVs need a circular assumption."""
    if "inlet_area_m2" not in frame:
        if fallback_radius is None or not np.isfinite(fallback_radius) or fallback_radius <= 0:
            raise ValueError("Missing mesh inlet area; supply an explicit legacy circular radius")
        return float(np.pi * fallback_radius**2)
    area = frame["inlet_area_m2"].to_numpy()
    if not np.isfinite(area).all() or np.any(area <= 0) or not np.allclose(area, area[0], rtol=1e-6, atol=1e-12):
        raise ValueError("Inconsistent or invalid mesh inlet areas")
    return float(area[0])
