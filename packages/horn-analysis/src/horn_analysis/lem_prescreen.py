"""LEM (Lumped Element Model) prescreening for auto pipeline.

Combines Webster/TMM throat impedance with driver coupling and scoring
to rapidly rank all geometry candidates analytically. Only the top-N
candidates proceed to expensive FEM simulation.

Lives in horn-analysis (not horn-core) because it imports transfer_function,
kpi, and scoring from horn-analysis.
"""

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from horn_core.candidates import CandidateGeometry
from horn_core.parameters import DriverParameters
from horn_core.profiles import get_radius_func
from horn_core.webster import compute_throat_impedance_tmm
from horn_analysis.kpi import extract_kpis_from_arrays
from horn_analysis.scoring import TargetSpec, compute_selection_score
from horn_analysis.transfer_function import compute_driver_response


P_REF = 20e-6  # reference pressure for dB SPL (Pa)


def _spl_from_pressure(p_throat: np.ndarray) -> np.ndarray:
    """Convert complex throat pressure to SPL (dB re 20 µPa)."""
    p_mag = np.abs(p_throat)
    p_mag = np.maximum(p_mag, 1e-30)
    return 20.0 * np.log10(p_mag / P_REF)


def lem_prescreen_candidates(
    candidates: List[CandidateGeometry],
    drivers: List[DriverParameters],
    target_f_low: float,
    target_f_high: float,
    sim_freq_range: tuple,
    num_frequencies: int = 100,
    top_n: int = 10,
    n_segments: int = 200,
) -> dict:
    """Score all (candidate, driver) pairs using TMM + driver coupling.

    For each pair:
    1. get_radius_func() → compute_throat_impedance_tmm() → z_real, z_imag
    2. compute_driver_response(driver, freq, z_real, z_imag, throat_area) → p_throat
    3. SPL = 20·log₁₀(|p_throat| / 20µPa)
    4. extract_kpis_from_arrays() → compute_selection_score() → composite score

    Ranks all scores, collects unique top-N candidate_ids.

    Args:
        candidates: List of CandidateGeometry from the candidates CSV.
        drivers: Pre-screened drivers.
        target_f_low: Target low frequency (Hz).
        target_f_high: Target high frequency (Hz).
        sim_freq_range: (min_freq, max_freq) for simulation.
        num_frequencies: Number of frequency points for TMM evaluation.
        top_n: Number of top candidates to pass to FEM.
        n_segments: Number of TMM segments per horn.

    Returns:
        Dict with keys: total_evaluated, total_pairs, top_n, rankings (list),
        filtered_candidate_ids (list of unique candidate IDs).
    """
    frequencies = np.linspace(sim_freq_range[0], sim_freq_range[1], num_frequencies)
    target = TargetSpec(f_low_hz=target_f_low, f_high_hz=target_f_high)

    all_scores = []

    # Cache TMM results per candidate (independent of driver)
    tmm_cache: Dict[str, tuple] = {}

    for cand in candidates:
        cache_key = cand.candidate_id
        if cache_key not in tmm_cache:
            radius_func = get_radius_func(
                cand.profile, cand.throat_radius, cand.mouth_radius, cand.length
            )
            z_real, z_imag = compute_throat_impedance_tmm(
                frequencies=frequencies,
                radius_func=radius_func,
                length=cand.length,
                throat_radius=cand.throat_radius,
                mouth_radius=cand.mouth_radius,
                n_segments=n_segments,
            )
            tmm_cache[cache_key] = (z_real, z_imag)

        z_real, z_imag = tmm_cache[cache_key]
        throat_area = np.pi * cand.throat_radius**2

        for drv in drivers:
            p_throat = compute_driver_response(
                drv, frequencies, z_real, z_imag, throat_area
            )
            coupled_spl = _spl_from_pressure(p_throat)
            kpi = extract_kpis_from_arrays(frequencies, coupled_spl)
            score = compute_selection_score(
                kpi, target,
                driver_id=drv.driver_id,
                horn_label=cand.candidate_id,
            )

            all_scores.append({
                "candidate_id": cand.candidate_id,
                "profile": cand.profile,
                "throat_radius": cand.throat_radius,
                "mouth_radius": cand.mouth_radius,
                "length": cand.length,
                "driver_id": drv.driver_id,
                "composite_score": score.composite_score,
                "bandwidth_coverage": score.bandwidth_coverage,
                "passband_ripple_db": score.passband_ripple_db,
                "avg_sensitivity_db": score.avg_sensitivity_db,
            })

    # Sort by composite score descending
    all_scores.sort(key=lambda x: x["composite_score"], reverse=True)

    # Collect unique top-N candidate IDs (by best score per candidate)
    seen_ids = set()
    filtered_ids = []
    for entry in all_scores:
        cid = entry["candidate_id"]
        if cid not in seen_ids:
            seen_ids.add(cid)
            filtered_ids.append(cid)
            if len(filtered_ids) >= top_n:
                break

    return {
        "total_evaluated": len(candidates),
        "total_pairs": len(all_scores),
        "top_n": top_n,
        "filtered_candidate_ids": filtered_ids,
        "rankings": all_scores,
    }


def _load_candidates_csv(csv_path: str) -> List[CandidateGeometry]:
    """Load candidates from a CSV file."""
    candidates = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            candidates.append(CandidateGeometry(
                candidate_id=row["candidate_id"],
                profile=row["profile"],
                throat_radius=float(row["throat_radius"]),
                mouth_radius=float(row["mouth_radius"]),
                length=float(row["length"]),
            ))
    return candidates


def _write_filtered_csv(
    candidates: List[CandidateGeometry],
    filtered_ids: List[str],
    output_path: str,
) -> None:
    """Write filtered candidates CSV containing only the top-N."""
    id_set = set(filtered_ids)
    filtered = [c for c in candidates if c.candidate_id in id_set]

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["candidate_id", "profile", "throat_radius", "mouth_radius", "length"])
        for c in filtered:
            writer.writerow([c.candidate_id, c.profile, c.throat_radius, c.mouth_radius, c.length])


def main():
    """CLI entry point for LEM prescreening."""
    parser = argparse.ArgumentParser(
        description="LEM/Webster prescreening for auto pipeline candidates.",
    )
    parser.add_argument("--candidates-csv", required=True, help="Input candidates CSV.")
    parser.add_argument("--prescreen-json", required=True, help="Driver prescreen result JSON.")
    parser.add_argument("--drivers-db", required=True, help="Driver database path.")
    parser.add_argument("--design-json", required=True, help="Design JSON with sim_freq_range.")
    parser.add_argument("--target-f-low", type=float, required=True, help="Target low freq (Hz).")
    parser.add_argument("--target-f-high", type=float, required=True, help="Target high freq (Hz).")
    parser.add_argument("--top-n", type=int, default=10, help="Number of top candidates for FEM.")
    parser.add_argument("--num-frequencies", type=int, default=100, help="Frequency points for TMM.")
    parser.add_argument("--output", required=True, help="Output LEM results JSON.")
    parser.add_argument("--filtered-csv", required=True, help="Output filtered candidates CSV.")
    args = parser.parse_args()

    from horn_drivers.loader import load_drivers

    # Load candidates
    candidates = _load_candidates_csv(args.candidates_csv)
    print(f"Loaded {len(candidates)} candidates from {args.candidates_csv}")

    # Load pre-screened driver IDs and full driver data
    prescreen = json.loads(Path(args.prescreen_json).read_text())
    driver_ids = set(prescreen["drivers"])
    all_drivers = load_drivers(args.drivers_db)
    drivers = [d for d in all_drivers if d.driver_id in driver_ids]
    print(f"Using {len(drivers)} pre-screened drivers")

    # Load simulation frequency range from design JSON
    design = json.loads(Path(args.design_json).read_text())
    sim_freq_range = tuple(design["sim_freq_range"])

    # Run LEM prescreening
    results = lem_prescreen_candidates(
        candidates=candidates,
        drivers=drivers,
        target_f_low=args.target_f_low,
        target_f_high=args.target_f_high,
        sim_freq_range=sim_freq_range,
        num_frequencies=args.num_frequencies,
        top_n=args.top_n,
    )

    # Write outputs
    Path(args.output).write_text(json.dumps(results, indent=2))
    _write_filtered_csv(candidates, results["filtered_candidate_ids"], args.filtered_csv)

    print(f"LEM prescreening complete:")
    print(f"  {results['total_evaluated']} candidates × {len(drivers)} drivers = {results['total_pairs']} pairs")
    print(f"  Top {len(results['filtered_candidate_ids'])} candidates passed to FEM")
    print(f"  Results: {args.output}")
    print(f"  Filtered CSV: {args.filtered_csv}")


if __name__ == "__main__":
    main()
