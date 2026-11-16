"""Couple a single driver with a single horn FEM result.

Reads the merged solver CSV produced by single-mode FEM, plus a driver from
the drivers DB, and produces:

  * coupled_spl.csv        : frequency, coupled_spl, horn_only_spl and the
                             electrical columns of electrical_output
  * coupled_spl.png        : plot comparing horn-only vs driver+horn
  * driver_horn_kpis.json  : evaluate_response KPIs on the coupled level

It uses the same observer path (evaluation.coupled_observer) and the same
assessment (evaluate_response through assess_with_chamber_policy) as auto
mode, so the same STEP response and driver give the same levels and KPIs.
"""
from horn_analysis.kpi import format_kpi

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from horn_drivers.loader import load_driver
from horn_analysis.kpi import extract_kpis_from_arrays
from horn_analysis.scoring import TargetSpec


def single_assessment(df, drv, target, legacy_radius=None, length=None):
    """Coupled levels, electrical table and auto-mode KPIs for one driver on one horn.

    The target band is clipped to the simulated sweep when the sweep does not
    bracket it; ``target_band_clipped_to_sweep`` records that. When the sweep
    misses the band entirely the whole sweep is assessed instead and
    ``target_band_outside_sweep`` records that.
    """
    from dataclasses import replace
    from horn_analysis.electrical_output import electrical_table
    from horn_analysis.evaluation import assess_with_chamber_policy, coupled_observer, evaluate_response
    from horn_analysis.transfer_function import compute_driver_operating_point

    f = df["frequency"].to_numpy(dtype=float)
    requested = [target.f_low_hz, target.f_high_hz]
    low, high = max(target.f_low_hz, f.min()), min(target.f_high_hz, f.max())
    outside = not low < high
    if outside:
        low, high = float(f.min()), float(f.max())
    clipped = not outside and (low, high) != (target.f_low_hz, target.f_high_hz)
    if outside or clipped:
        target = replace(target, f_low_hz=low, f_high_hz=high)
    observed = coupled_observer(df, drv, target, legacy_radius)
    mouth_radius = float(np.sqrt(df.mouth_area_m2.iloc[0] / np.pi)) if "mouth_area_m2" in df else None

    def evaluate(model):
        levels = observed["levels"] if model is drv else coupled_observer(df, model, target, legacy_radius)["levels"]
        point = compute_driver_operating_point(model, f, df["z_real"].to_numpy(), df["z_imag"].to_numpy(),
                                               observed["area"], target.voltage_rms)
        return evaluate_response(f, levels, target, model, observed["area"], length=length,
                                 mouth_radius=mouth_radius, operating_point=point)

    assessment = assess_with_chamber_policy(drv, evaluate, target)
    assessment.update(target_band_hz=[target.f_low_hz, target.f_high_hz], target_band_clipped_to_sweep=clipped,
                      target_band_outside_sweep=outside, target_band_requested_hz=requested,
                      output_metric=observed["metric"], inlet_area_m2=observed["area"])
    table = electrical_table(df, drv, target, legacy_radius, length)
    return observed["levels"], table, assessment


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
    target: TargetSpec | None = None,
    length: float | None = None,
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
    if target is None:
        target = TargetSpec(float(freq.min()), float(freq.max()), voltage_rms=voltage)
    elif target.voltage_rms != voltage:
        from dataclasses import replace
        target = replace(target, voltage_rms=voltage)

    coupled_spl, table, assessment = single_assessment(
        df, drv, target, None if imported_geometry else throat_radius, length)
    throat_area = assessment["inlet_area_m2"]

    out = pd.DataFrame({"frequency": freq, "coupled_spl": coupled_spl, "horn_only_spl": solver_spl})
    out = pd.concat([out, table.drop(columns=["frequency", "spl"])], axis=1)
    out.to_csv(output_csv, index=False)

    kpi = extract_kpis_from_arrays(freq, coupled_spl)
    kpi_payload = {
        **kpi.to_dict(),
        # Band-referenced KPIs, f3 and ripple from evaluate_response, as in auto mode.
        **assessment,
        "average_sensitivity_db": assessment["avg_sensitivity_db"],
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
    if not assessment["f3_low_is_bound"]:
        ax.axvline(assessment["f3_low_hz"], color="green", alpha=0.4, linestyle=":", label=f"f3 low {assessment['f3_low_hz']:.0f} Hz")
    if not assessment["f3_high_is_bound"]:
        ax.axvline(assessment["f3_high_hz"], color="red", alpha=0.4, linestyle=":", label=f"f3 high {assessment['f3_high_hz']:.0f} Hz")
    ax.axvspan(*assessment["target_band_hz"], color="#94a3b8", alpha=0.12, label="Target band")

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Mouth-plane pressure level (dB re 20 µPa)" if assessment["output_metric"] == "legacy_mouth_pressure"
                  else f"On-axis level at {target.observation_distance_m:g} m (dB re 20 µPa)")
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

    f3lo = format_kpi(assessment, "f3_low_hz")
    f3hi = format_kpi(assessment, "f3_high_hz")
    rip = f"{assessment['passband_ripple_db']:.2f}"
    print(f"Wrote {output_csv}, {output_png}, {output_kpis}")
    print(f"Driver: {drv.manufacturer} {drv.model_name}  Fs={drv.fs_hz} Hz  Le={drv.le_h*1000:.2f} mH  Bl={drv.bl_tm}")
    print(f"KPI: f3 {f3lo}-{f3hi} Hz  peak {kpi.peak_spl_db:.1f} dB @ {kpi.peak_frequency_hz:.0f} Hz  ripple {rip} dB  "
          f"status {assessment['eligibility_status']}" + (" (estimated front chamber)" if assessment["front_chamber_estimated"] else ""))


def main():
    parser = argparse.ArgumentParser(description="Couple a driver with a horn FEM result.")
    parser.add_argument("--solver-csv", required=True, help="Merged FEM solver CSV (frequency, spl, z_real, z_imag).")
    parser.add_argument("--drivers-db", required=True, help="Driver database directory or JSON file.")
    parser.add_argument("--driver-id", required=True, help="Driver ID to load from the database.")
    parser.add_argument("--throat-radius", type=float, default=None, help="Known circular throat radius (m); omitted for imported CAD.")
    parser.add_argument("--imported-geometry", action="store_true")
    parser.add_argument("--profile", default="horn", help="Horn profile label (cosmetic).")
    parser.add_argument("--voltage", type=float, default=2.83, help="Drive voltage (V).")
    parser.add_argument("--target-f-low", type=float, default=None, help="Target band low edge (Hz); default sweep start.")
    parser.add_argument("--target-f-high", type=float, default=None, help="Target band high edge (Hz); default sweep end.")
    parser.add_argument("--distance", type=float, default=1.0, help="On-axis observer distance (m).")
    parser.add_argument("--max-ripple", type=float, default=6.0)
    parser.add_argument("--max-compression", type=float, default=10.0)
    parser.add_argument("--test-power", type=float, default=None)
    parser.add_argument("--min-output", type=float, default=None, help="Minimum mean level in the band (dB), as in auto mode.")
    parser.add_argument("--score-weights", default=None, help="Score weights, e.g. ripple=1,extension=0.5, as in auto mode.")
    parser.add_argument("--length", type=float, default=None, help="Horn length (m), for time of flight and the size term.")
    parser.add_argument("--output-csv", default="coupled_spl.csv")
    parser.add_argument("--output-png", default="coupled_spl.png")
    parser.add_argument("--output-kpis", default="driver_horn_kpis.json")
    args = parser.parse_args()

    sweep = pd.read_csv(args.solver_csv, usecols=["frequency"])["frequency"]
    target = TargetSpec(args.target_f_low if args.target_f_low is not None else float(sweep.min()),
                        args.target_f_high if args.target_f_high is not None else float(sweep.max()),
                        voltage_rms=args.voltage, observation_distance_m=args.distance,
                        max_ripple_db=args.max_ripple, max_compression_ratio=args.max_compression,
                        min_output_db=args.min_output, score_weights=args.score_weights,
                        test_power_w=args.test_power)
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
        target=target,
        length=args.length,
    )


if __name__ == "__main__":
    main()
