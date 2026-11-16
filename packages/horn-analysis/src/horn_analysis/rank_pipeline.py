"""Integrated ranking pipeline for driver-horn combinations.

Takes a solver CSV, pre-screened drivers, and a target spec, then scores
and ranks all combinations using coupled SPL from the transfer function.
"""
from horn_core.acoustics import inlet_area_from_frame

import argparse
import json
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_analysis.kpi import extract_kpis_from_arrays
from horn_analysis.scoring import TargetSpec


def rank_horn_drivers(
    solver_csv: str,
    horn_label: str,
    throat_radius: float,
    drivers: List[DriverParameters],
    target: TargetSpec,
    top_n: int = 10,
    v_g: float = None,
    length: float = None,
) -> List[dict]:
    """Score all drivers against one horn geometry using coupled SPL.

    Per driver:
    1. coupled_output() -> observer level from the solver transfer
    2. compute_driver_operating_point() -> excursion and input power
    3. evaluate_response() -> constraints, f3, extension and composite score
    4. extract_kpis_from_arrays() -> peak-lobe diagnostics

    Args:
        solver_csv: Path to solver output CSV (frequency, spl, z_real, z_imag).
        horn_label: Label for this horn geometry (e.g. "exponential").
        throat_radius: Horn throat radius in metres.
        drivers: List of DriverParameters to evaluate.
        target: Target frequency specification.
        top_n: Number of top results to return.
        v_g: Input voltage (V).
        length: Horn length in metres, for the length cap and size term.

    Returns:
        List of dicts with score + KPI info for top N drivers, sorted by
        composite_score descending.
    """
    from horn_analysis.evaluation import coupled_output
    if top_n < 1:
        raise ValueError("top_n must be positive")
    if v_g is not None:
        from dataclasses import replace
        target = replace(target, voltage_rms=v_g)
    df = pd.read_csv(solver_csv)
    results = []
    from horn_analysis.evaluation import assess_driver
    from horn_analysis.transfer_function import compute_driver_operating_point
    for drv in drivers:
        coupled_spl, area, metric = coupled_output(df, drv, target, throat_radius)

        def response(driver):
            levels = coupled_spl if driver is drv else coupled_output(df, driver, target, throat_radius)[0]
            point = compute_driver_operating_point(driver,df.frequency.to_numpy(),df.z_real.to_numpy(),df.z_imag.to_numpy(),area,target.voltage_rms)
            return levels, point
        mouth_radius = float(np.sqrt(df.mouth_area_m2.iloc[0] / np.pi)) if "mouth_area_m2" in df else None
        assessment = assess_driver(df.frequency.to_numpy(), target, drv, area, response, length=length, mouth_radius=mouth_radius)
        kpi = extract_kpis_from_arrays(df.frequency.to_numpy(), coupled_spl)
        results.append({
            **assessment, "driver_id": drv.driver_id, "horn_label": horn_label,
            "kpi": kpi.to_dict(), "manufacturer": drv.manufacturer, "model_name": drv.model_name,
            "output_metric": metric, "observation_distance_m": target.observation_distance_m,
            "drive_voltage_rms": target.voltage_rms, "inlet_area_m2": area,
            "loss_model": str(df.loss_model.iloc[0]) if "loss_model" in df else "lossless",
            "radiation_model": str(df.radiation_model.iloc[0]) if "radiation_model" in df else "legacy_unknown",
            "element_degree": int(df.element_degree.iloc[0]) if "element_degree" in df else 1,
        })
    return sorted(results, key=lambda r: (-r["composite_score"], r["driver_id"]))[:top_n]


def rejection_rows(record, candidate, driver_ids):
    """Ranking rows for a candidate the solver rejected before producing a response.

    ``record`` is a merged rejection (horn_analysis.merge.merge_rejections) and
    ``candidate`` its candidates.csv row. One row per driver keeps the rejected
    list in the same shape as model rejections.
    """
    geometry = {key: float(candidate[key]) for key in ("throat_radius", "mouth_radius", "length")
                if candidate.get(key) not in (None, "")}
    return [{"horn_label": record["candidate_id"], "driver_id": driver_id,
             "profile": candidate.get("profile", ""), **geometry,
             "model_feasible": False, "simulation_eligible": False, "composite_score": 0.,
             "extension_octaves": None, "passband_ripple_db": None, "avg_sensitivity_db": None,
             "rejection_reasons": list(record["rejection_reasons"]),
             "rejection_detail": record.get("rejection_detail", ""), "evidence_gaps": []}
            for driver_id in sorted(driver_ids)]


def main():
    """CLI for ranking drivers against a horn solver result."""
    parser = argparse.ArgumentParser(
        description="Rank drivers against a horn solver result.",
    )
    parser.add_argument("--solver-csv", required=True, help="Solver output CSV.")
    parser.add_argument("--drivers-db", required=True, help="Driver database JSON.")
    parser.add_argument("--throat-radius", type=float, required=True, help="Throat radius (m).")
    parser.add_argument("--target-f-low", type=float, required=True, help="Target low freq (Hz).")
    parser.add_argument("--target-f-high", type=float, required=True, help="Target high freq (Hz).")
    parser.add_argument("--top-n", type=int, default=10, help="Number of top results.")
    parser.add_argument("--horn-label", type=str, default="horn", help="Horn geometry label.")
    parser.add_argument("--voltage", type=float, default=2.83, help="Input voltage (V).")
    parser.add_argument("--min-output", type=float, default=None, help="Minimum target-band level (dB).")
    parser.add_argument("--score-weights", default=None, help="e.g. sensitivity=0.625,ripple=0.25,extension=0.125,size=0")
    parser.add_argument("--output", type=str, default="ranked_results.json", help="Output JSON file.")
    args = parser.parse_args()

    from horn_drivers.loader import load_drivers

    drivers = load_drivers(args.drivers_db)
    print(f"Loaded {len(drivers)} drivers from {args.drivers_db}")

    target = TargetSpec(f_low_hz=args.target_f_low, f_high_hz=args.target_f_high,
                        min_output_db=args.min_output, score_weights=args.score_weights)

    results = rank_horn_drivers(
        solver_csv=args.solver_csv,
        horn_label=args.horn_label,
        throat_radius=args.throat_radius,
        drivers=drivers,
        target=target,
        top_n=args.top_n,
        v_g=args.voltage,
    )

    output = json.dumps(results, indent=2)
    Path(args.output).write_text(output)
    print(f"\nTop {len(results)} results written to {args.output}")
    for i, r in enumerate(results, 1):
        print(f"  {i}. {r['driver_id']} ({r['manufacturer']} {r['model_name']}) "
              f"score={r['composite_score']:.3f}")


if __name__ == "__main__":
    main()
