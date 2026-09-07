"""Driver–horn transfer function for Phase A screening.

Couples a compression driver model (from T-S parameters) with the horn's
throat impedance (from the FEM solver CSV) to compute the actual SPL at
the horn mouth without re-running the FEM solver.

Physics chain (impedance domains documented at every step)::

    V_g  (input voltage, default 2.83 V RMS)
      → Z_e = Re + jωLe                                    (electrical, Ω)
      → Z_mech_driver = Rms + jωMms + 1/(jωCms)            (mechanical, kg/s)
      → Z_horn = z_real + j·z_imag                          (specific acoustic, Pa·s/m)
      → Z_mech_load = Z_horn · Sd² / S_throat               (mechanical, kg/s)
      → Z_mech_total = Z_mech_driver + Z_mech_load
      → Z_mot = BL² / Z_mech_total                          (motional impedance, Ω)
      → I = V_g / (Z_e + Z_mot)                             (current, A)
      → v = BL·I / Z_mech_total                             (diaphragm velocity, m/s)
      → p_throat = Z_horn · v · (Sd / S_throat)             (throat pressure, Pa)
"""
from horn_core.acoustics import inlet_area_from_frame

import argparse
import json
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_drivers.loader import load_drivers


def compute_driver_operating_point(
    driver: DriverParameters,
    frequencies: np.ndarray,
    z_horn_real: np.ndarray,
    z_horn_imag: np.ndarray,
    throat_area: float,
    v_g: float = 2.83,
) -> dict:
    """Compute the linear motor operating point for a driver coupled to a horn.

    Args:
        driver: Driver T-S parameters (SI units).
        frequencies: Array of frequencies in Hz.
        z_horn_real: Real part of specific acoustic throat impedance (Pa·s/m).
        z_horn_imag: Imaginary part of specific acoustic throat impedance (Pa·s/m).
        throat_area: Physical throat cross-section area (m²).
        v_g: Input voltage (V). Default 2.83 V RMS.

    Returns:
        Complex RMS pressure, velocity and current arrays, peak displacement,
        electrical impedance, and real power terms at each frequency.
    """
    if not np.isfinite(throat_area) or throat_area <= 0 or not np.isfinite(v_g) or v_g <= 0:
        raise ValueError("Throat area and RMS voltage must be positive and finite")
    frequencies = np.asarray(frequencies, dtype=float)
    z_horn_real = np.asarray(z_horn_real, dtype=float)
    z_horn_imag = np.asarray(z_horn_imag, dtype=float)
    if (frequencies.ndim != 1 or frequencies.size == 0
            or np.any(frequencies <= 0) or not np.isfinite(frequencies).all()
            or z_horn_real.shape != frequencies.shape or z_horn_imag.shape != frequencies.shape
            or not np.isfinite(z_horn_real).all() or not np.isfinite(z_horn_imag).all()):
        raise ValueError("Positive finite frequencies and matching finite impedance arrays required")
    # Keep legacy Mms estimates explicitly provisional. When known, use Mmd
    # plus separately specified rear air load; never subtract an invented mass.
    moving_mass = driver.coupled_moving_mass_kg
    if moving_mass <= 0 or driver.cms_m_per_n <= 0:
        raise ValueError("Invalid driver mass or compliance")
    omega = 2.0 * np.pi * frequencies

    # Electrical impedance
    z_e = driver.re_ohm + 1j * omega * driver.le_h

    # Mechanical impedance of the driver suspension
    rms = driver.rms_kg_per_s if driver.rms_kg_per_s is not None else 0.0
    z_mech_driver = rms + 1j * omega * moving_mass + 1.0 / (1j * omega * driver.cms_m_per_n)

    # Horn throat impedance (specific acoustic → mechanical via area ratio)
    z_horn = z_horn_real + 1j * z_horn_imag
    z_mech_load = z_horn * (driver.sd_m2 ** 2) / throat_area

    # Total mechanical impedance seen by the voice coil
    z_mech_total = z_mech_driver + z_mech_load

    # Solve the two motor equations without dividing by mechanical impedance.
    # The undamped, unloaded resonance has zero current and finite velocity.
    denominator = z_e * z_mech_total + driver.bl_tm**2
    current = v_g * z_mech_total / denominator
    velocity = driver.bl_tm * v_g / denominator
    with np.errstate(divide="ignore", invalid="ignore"):
        electrical_impedance = np.divide(v_g, current)
    electrical_impedance = np.where(current == 0, complex(np.inf, 0), electrical_impedance)

    # Throat pressure: p = Z_horn * U / S_throat, where U = v * Sd
    p_throat = z_horn * velocity * (driver.sd_m2 / throat_area)

    return {"throat_pressure": p_throat, "velocity_rms": velocity,
            "current_rms": current, "electrical_impedance": electrical_impedance,
            "displacement_peak_m": np.sqrt(2)*np.abs(velocity)/omega,
            "input_power_w": np.real(v_g*np.conjugate(current)),
            "copper_power_w": np.abs(current)**2*driver.re_ohm,
            "mechanical_loss_w": np.abs(velocity)**2*rms,
            "horn_power_w": np.abs(velocity)**2*np.real(z_mech_load)}


def compute_driver_response(driver, frequencies, z_horn_real, z_horn_imag, throat_area, v_g=2.83):
    """Compatibility API returning complex throat pressure at RMS drive voltage."""
    return compute_driver_operating_point(driver,frequencies,z_horn_real,z_horn_imag,throat_area,v_g)["throat_pressure"]


def scale_solver_spl(
    solver_spl: np.ndarray,
    p_throat: np.ndarray,
) -> np.ndarray:
    """Scale normalised solver SPL by the actual throat pressure.

    The FEM solver runs with p=1 at the inlet (Dirichlet).  The actual
    SPL is obtained by adding the dB-level of the real throat pressure::

        SPL_actual = SPL_solver + 20·log10(|p_throat|)

    Args:
        solver_spl: SPL array from the Dirichlet-BC solver (dB re 20 µPa).
        p_throat: Complex throat pressure array from ``compute_driver_response``.

    Returns:
        Scaled SPL array (dB re 20 µPa).
    """
    p_magnitude = np.abs(p_throat)
    # Avoid log(0)
    p_magnitude = np.maximum(p_magnitude, 1e-30)
    return solver_spl + 20.0 * np.log10(p_magnitude)


def screen_all_drivers(
    solver_csv: str,
    drivers: List[DriverParameters],
    throat_radius: float,
    v_g: float = 2.83,
) -> pd.DataFrame:
    """Batch-screen all drivers against a single horn geometry.

    Reads the solver CSV (frequency, spl, z_real, z_imag), computes the
    coupled SPL for every driver, and returns a DataFrame with per-driver
    results and basic KPIs.

    Args:
        solver_csv: Path to CSV with columns: frequency, spl, z_real, z_imag.
        drivers: List of DriverParameters to screen.
        throat_radius: Horn throat radius in metres.
        v_g: Input voltage (V).

    Returns:
        DataFrame with columns: driver_id, manufacturer, model_name,
        avg_spl_db, peak_spl_db, min_spl_db, ripple_db.
    """
    df = pd.read_csv(solver_csv)
    freq = df["frequency"].values
    solver_spl = df["spl"].values
    z_real = df["z_real"].values
    z_imag = df["z_imag"].values
    throat_area = inlet_area_from_frame(df, throat_radius)

    rows = []
    for drv in drivers:
        p_throat = compute_driver_response(drv, freq, z_real, z_imag, throat_area, v_g)
        scaled_spl = scale_solver_spl(solver_spl, p_throat)

        rows.append({
            "driver_id": drv.driver_id,
            "manufacturer": drv.manufacturer,
            "model_name": drv.model_name,
            "avg_spl_db": float(np.mean(scaled_spl)),
            "peak_spl_db": float(np.max(scaled_spl)),
            "min_spl_db": float(np.min(scaled_spl)),
            "ripple_db": float(np.max(scaled_spl) - np.min(scaled_spl)),
        })

    return pd.DataFrame(rows)


def main():
    """CLI for driver screening against a solver result."""
    parser = argparse.ArgumentParser(
        description="Screen drivers against horn solver results (Phase A).",
    )
    parser.add_argument("--solver-csv", required=True, help="Solver output CSV.")
    parser.add_argument("--drivers-db", required=True, help="Driver database JSON.")
    parser.add_argument("--throat-radius", type=float, required=True, help="Throat radius (m).")
    parser.add_argument("--voltage", type=float, default=2.83, help="Input voltage (V).")
    parser.add_argument("--output-dir", type=str, default="screening_results", help="Output directory.")
    args = parser.parse_args()

    drivers = load_drivers(args.drivers_db)
    print(f"Loaded {len(drivers)} drivers from {args.drivers_db}")

    results = screen_all_drivers(
        args.solver_csv, drivers, args.throat_radius, args.voltage,
    )

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    results.to_csv(out / "screening_summary.csv", index=False)
    print(f"\nScreening results written to {out / 'screening_summary.csv'}")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
