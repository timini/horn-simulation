"""Couple a single driver with a single horn FEM result.

Reads the merged solver CSV (frequency, spl, z_real, z_imag) produced by
single-mode FEM, plus a driver from the drivers DB, and produces:

  * coupled_spl.csv        — frequency, coupled_spl, horn_only_spl
  * coupled_spl.png        — plot comparing horn-only vs driver+horn
  * driver_horn_kpis.json  — KPIs computed on the coupled SPL

This lets `nextflow run main.nf --mode single --driver_id …` give an
end-to-end driver+horn frequency response without going through auto mode.
"""
from horn_core.acoustics import inlet_area_from_frame

from horn_analysis.kpi import format_kpi

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from horn_drivers.loader import load_driver
from horn_analysis.kpi import extract_kpis_from_arrays
from horn_analysis.transfer_function import compute_driver_response, scale_solver_spl


def couple(
    solver_csv: str,
    drivers_db: str,
    driver_id: str,
    throat_radius: float | None = None,
    profile: str = "horn",
    voltage: float = 2.83,
    output_csv: str = "coupled_spl.csv",
    output_png: str = "coupled_spl.png",
    output_kpis: str = "driver_horn_kpis.json",
    imported_geometry: bool = False,
):
    """Compute and write the coupled driver+horn response."""
    drv = load_driver(drivers_db, driver_id)

    df = pd.read_csv(solver_csv)
    if "schema_version" in df and (
        not (df.schema_version == 2).all() or "bc_mode" not in df or not (df.bc_mode == "dirichlet").all()
    ):
        raise ValueError("Single-driver coupling requires a unit-pressure Dirichlet transfer run")
    freq = df["frequency"].values
    solver_spl = df["spl"].values
    z_real = df["z_real"].values
    z_imag = df["z_imag"].values

    throat_area = inlet_area_from_frame(df, None if imported_geometry else throat_radius)

    p_throat = compute_driver_response(drv, freq, z_real, z_imag, throat_area, voltage)
    coupled_spl = scale_solver_spl(solver_spl, p_throat)

    out = pd.DataFrame({
        "frequency": freq,
        "coupled_spl": coupled_spl,
        "horn_only_spl": solver_spl,
    })
    out.to_csv(output_csv, index=False)

    kpi = extract_kpis_from_arrays(freq, coupled_spl)
    kpi_payload = {
        **kpi.to_dict(),
        "driver_id": drv.driver_id,
        "manufacturer": drv.manufacturer,
        "model_name": drv.model_name,
        "fs_hz": drv.fs_hz,
        "qts": drv.qts,
        "bl_tm": drv.bl_tm,
        "le_h": drv.le_h,
        "mms_kg": drv.mms_kg,
        "sd_m2": drv.sd_m2,
        "throat_radius_m": None if imported_geometry else throat_radius,
        "inlet_area_m2": throat_area,
        "drive_voltage_v": voltage,
    }
    Path(output_kpis).write_text(json.dumps(kpi_payload, indent=2))

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 6.5))
    label = f"{drv.manufacturer} {drv.model_name} + {profile} horn"
    ax.semilogx(freq, coupled_spl, label=label, linewidth=2.0, color="#1f77b4")
    ax.semilogx(freq, solver_spl, label="Horn at unit inlet pressure (1 Pa RMS)",
                linewidth=1.0, color="#888888", linestyle="--")

    # Mark KPI points
    if kpi.f3_low_hz and not kpi.f3_low_is_bound:
        ax.axvline(kpi.f3_low_hz, color="green", alpha=0.4, linestyle=":", label=f"f3 low {kpi.f3_low_hz:.0f} Hz")
    if kpi.f3_high_hz and not kpi.f3_high_is_bound:
        ax.axvline(kpi.f3_high_hz, color="red", alpha=0.4, linestyle=":", label=f"f3 high {kpi.f3_high_hz:.0f} Hz")

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Mouth-plane pressure level (dB re 20 µPa)")
    interface_label = (f"Inlet area {throat_area:.8g} m²" if imported_geometry or throat_radius is None
                       else f"Throat radius {throat_radius * 1000:.1f} mm")
    ax.set_title(
        f"Driver + horn coupled response — {drv.manufacturer} {drv.model_name}\n"
        f"{interface_label}; driver drive {voltage:g} V RMS"
    )
    ax.legend(loc="lower center", fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(50, 20000)
    plt.tight_layout()
    plt.savefig(output_png, dpi=120)
    plt.close()

    f3lo = format_kpi(kpi, "f3_low_hz")
    f3hi = format_kpi(kpi, "f3_high_hz")
    rip = f"{kpi.passband_ripple_db:.2f}" if kpi.passband_ripple_db else "n/a"
    print(f"Wrote {output_csv}, {output_png}, {output_kpis}")
    print(f"Driver: {drv.manufacturer} {drv.model_name}  Fs={drv.fs_hz} Hz  Le={drv.le_h*1000:.2f} mH  Bl={drv.bl_tm}")
    print(f"KPI: f3 {f3lo}-{f3hi} Hz  peak {kpi.peak_spl_db:.1f} dB @ {kpi.peak_frequency_hz:.0f} Hz  ripple {rip} dB")


def main():
    parser = argparse.ArgumentParser(description="Couple a driver with a horn FEM result.")
    parser.add_argument("--solver-csv", required=True, help="Merged FEM solver CSV (frequency, spl, z_real, z_imag).")
    parser.add_argument("--drivers-db", required=True, help="Driver database directory or JSON file.")
    parser.add_argument("--driver-id", required=True, help="Driver ID to load from the database.")
    parser.add_argument("--throat-radius", type=float, default=None, help="Known circular throat radius (m); omitted for imported CAD.")
    parser.add_argument("--imported-geometry", action="store_true")
    parser.add_argument("--profile", default="horn", help="Horn profile label (cosmetic).")
    parser.add_argument("--voltage", type=float, default=2.83, help="Drive voltage (V).")
    parser.add_argument("--output-csv", default="coupled_spl.csv")
    parser.add_argument("--output-png", default="coupled_spl.png")
    parser.add_argument("--output-kpis", default="driver_horn_kpis.json")
    args = parser.parse_args()

    couple(
        solver_csv=args.solver_csv,
        drivers_db=args.drivers_db,
        driver_id=args.driver_id,
        throat_radius=args.throat_radius,
        profile=args.profile,
        voltage=args.voltage,
        output_csv=args.output_csv,
        output_png=args.output_png,
        output_kpis=args.output_kpis,
        imported_geometry=args.imported_geometry,
    )


if __name__ == "__main__":
    main()
