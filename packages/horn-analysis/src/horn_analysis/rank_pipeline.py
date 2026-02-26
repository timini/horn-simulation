"""Integrated ranking pipeline for driver-horn combinations.

Takes a solver CSV, pre-screened drivers, and a target spec, then scores
and ranks all combinations using coupled SPL from the transfer function.
"""

import argparse
import json
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_analysis.kpi import extract_kpis_from_arrays
from horn_analysis.scoring import TargetSpec, SelectionScore, compute_selection_score, rank_candidates
from horn_analysis.transfer_function import compute_driver_response, scale_solver_spl


def rank_horn_drivers(
    solver_csv: str,
    horn_label: str,
    throat_radius: float,
    drivers: List[DriverParameters],
    target: TargetSpec,
    top_n: int = 10,
    v_g: float = 2.83,
) -> List[dict]:
    """Score all drivers against one horn geometry using coupled SPL.

    Per driver:
    1. compute_driver_response() -> p_throat
    2. scale_solver_spl() -> coupled SPL array
    3. extract_kpis_from_arrays() -> HornKPI
    4. compute_selection_score() -> SelectionScore

    Args:
        solver_csv: Path to solver output CSV (frequency, spl, z_real, z_imag).
        horn_label: Label for this horn geometry (e.g. "exponential").
        throat_radius: Horn throat radius in metres.
        drivers: List of DriverParameters to evaluate.
        target: Target frequency specification.
        top_n: Number of top results to return.
        v_g: Input voltage (V).

    Returns:
        List of dicts with score + KPI info for top N drivers, sorted by
        composite_score descending.
    """
    df = pd.read_csv(solver_csv)
    freq = df["frequency"].values
    solver_spl = df["spl"].values
    z_real = df["z_real"].values
    z_imag = df["z_imag"].values
    throat_area = np.pi * throat_radius ** 2

    scores: List[SelectionScore] = []

    for drv in drivers:
        p_throat = compute_driver_response(
            drv, freq, z_real, z_imag, throat_area, v_g,
        )
        coupled_spl = scale_solver_spl(solver_spl, p_throat)
        kpi = extract_kpis_from_arrays(freq, coupled_spl)
        score = compute_selection_score(
            kpi, target,
            driver_id=drv.driver_id,
            horn_label=horn_label,
        )
        scores.append(score)

    ranked = rank_candidates(scores, top_n=top_n)

    results = []
    for s in ranked:
        # Re-compute KPI for the result dict
        drv = next(d for d in drivers if d.driver_id == s.driver_id)
        p_throat = compute_driver_response(
            drv, freq, z_real, z_imag, throat_area, v_g,
        )
        coupled_spl = scale_solver_spl(solver_spl, p_throat)
        kpi = extract_kpis_from_arrays(freq, coupled_spl)

        results.append({
            **s.to_dict(),
            "kpi": kpi.to_dict(),
            "manufacturer": drv.manufacturer,
            "model_name": drv.model_name,
        })

    return results


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
    parser.add_argument("--output", type=str, default="ranked_results.json", help="Output JSON file.")
    args = parser.parse_args()

    from horn_drivers.loader import load_drivers

    drivers = load_drivers(args.drivers_db)
    print(f"Loaded {len(drivers)} drivers from {args.drivers_db}")

    target = TargetSpec(f_low_hz=args.target_f_low, f_high_hz=args.target_f_high)

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
