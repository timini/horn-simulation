"""Driver-horn transfer function coupling.

Couples a compression driver model (from T-S parameters) with the horn's
throat impedance (from the FEM solver CSV) to compute the actual SPL at
the horn mouth without re-running the FEM solver.

Physics chain (impedance domains documented at every step)::

    V_g  (input voltage, default 2.83 V = 1 W into 8 Ohm)
      -> Z_e = Re + j*omega*Le                              (electrical, Ohm)
      -> Z_mech_driver = Rms + j*omega*Mms + 1/(j*omega*Cms) (mechanical, kg/s)
      -> Z_horn = z_real + j*z_imag                          (specific acoustic, Pa*s/m)
      -> Z_mech_load = Z_horn * Sd^2 / S_throat              (mechanical, kg/s)
      -> Z_mech_total = Z_mech_driver + Z_mech_load
      -> Z_mot = BL^2 / Z_mech_total                         (motional impedance, Ohm)
      -> I = V_g / (Z_e + Z_mot)                             (current, A)
      -> v = BL*I / Z_mech_total                              (diaphragm velocity, m/s)
      -> p_throat = Z_horn * v * (Sd / S_throat)              (throat pressure, Pa)
"""

import argparse
import json
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters


def compute_driver_response(
    driver: DriverParameters,
    frequencies: np.ndarray,
    z_horn_real: np.ndarray,
    z_horn_imag: np.ndarray,
    throat_area: float,
    v_g: float = 2.83,
) -> np.ndarray:
    """Compute complex throat pressure for a driver coupled to a horn.

    Args:
        driver: Driver T-S parameters (SI units).
        frequencies: Array of frequencies in Hz.
        z_horn_real: Real part of specific acoustic throat impedance (Pa*s/m).
        z_horn_imag: Imaginary part of specific acoustic throat impedance (Pa*s/m).
        throat_area: Physical throat cross-section area (m^2).
        v_g: Input voltage (V). Default 2.83 V (1 W into 8 Ohm).

    Returns:
        Complex array of throat pressures (Pa) at each frequency.
    """
    omega = 2.0 * np.pi * frequencies

    # Electrical impedance
    z_e = driver.re_ohm + 1j * omega * driver.le_h

    # Mechanical impedance of the driver suspension
    rms = driver.rms_kg_per_s if driver.rms_kg_per_s is not None else 0.0
    z_mech_driver = rms + 1j * omega * driver.mms_kg + 1.0 / (1j * omega * driver.cms_m_per_n)

    # Horn throat impedance (specific acoustic -> mechanical via area ratio)
    z_horn = z_horn_real + 1j * z_horn_imag
    z_mech_load = z_horn * (driver.sd_m2 ** 2) / throat_area

    # Total mechanical impedance seen by the voice coil
    z_mech_total = z_mech_driver + z_mech_load

    # Motional impedance reflected into electrical domain
    z_mot = (driver.bl_tm ** 2) / z_mech_total

    # Voice coil current
    current = v_g / (z_e + z_mot)

    # Diaphragm velocity
    velocity = driver.bl_tm * current / z_mech_total

    # Throat pressure: p = Z_horn * U / S_throat, where U = v * Sd
    p_throat = z_horn * velocity * (driver.sd_m2 / throat_area)

    return p_throat


def scale_solver_spl(
    solver_spl: np.ndarray,
    p_throat: np.ndarray,
) -> np.ndarray:
    """Scale normalised solver SPL by the actual throat pressure.

    The FEM solver runs with p=1 at the inlet (Dirichlet).  The actual
    SPL is obtained by adding the dB-level of the real throat pressure::

        SPL_actual = SPL_solver + 20*log10(|p_throat|)

    Args:
        solver_spl: SPL array from the Dirichlet-BC solver (dB re 20 uPa).
        p_throat: Complex throat pressure array from ``compute_driver_response``.

    Returns:
        Scaled SPL array (dB re 20 uPa).
    """
    p_magnitude = np.abs(p_throat)
    p_magnitude = np.maximum(p_magnitude, 1e-30)
    return solver_spl + 20.0 * np.log10(p_magnitude)


def main():
    """CLI for driver screening against a solver result."""
    parser = argparse.ArgumentParser(
        description="Screen drivers against horn solver results.",
    )
    parser.add_argument("--solver-csv", required=True, help="Solver output CSV.")
    parser.add_argument("--drivers-db", required=True, help="Driver database JSON.")
    parser.add_argument("--throat-radius", type=float, required=True, help="Throat radius (m).")
    parser.add_argument("--voltage", type=float, default=2.83, help="Input voltage (V).")
    parser.add_argument("--output-dir", type=str, default="screening_results", help="Output directory.")
    args = parser.parse_args()

    from horn_drivers.loader import load_drivers

    drivers = load_drivers(args.drivers_db)
    print(f"Loaded {len(drivers)} drivers from {args.drivers_db}")

    df = pd.read_csv(args.solver_csv)
    freq = df["frequency"].values
    solver_spl = df["spl"].values
    z_real = df["z_real"].values
    z_imag = df["z_imag"].values
    throat_area = np.pi * args.throat_radius ** 2

    rows = []
    for drv in drivers:
        p_throat = compute_driver_response(drv, freq, z_real, z_imag, throat_area, args.voltage)
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

    results = pd.DataFrame(rows)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    results.to_csv(out / "screening_summary.csv", index=False)
    print(f"\nScreening results written to {out / 'screening_summary.csv'}")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
