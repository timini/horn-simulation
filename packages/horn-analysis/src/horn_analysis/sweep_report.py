"""Sweep report generation for driver-horn ranking results.

Produces ranking JSON, comparison plots, individual CSVs, and a
human-readable summary from ranked driver-horn combinations.
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_analysis.compare import plot_multi_comparison
from horn_analysis.scoring import TargetSpec
from horn_analysis.transfer_function import compute_driver_response, scale_solver_spl


def generate_sweep_report(
    all_ranked: List[dict],
    solver_csvs: Dict[str, str],
    drivers: Dict[str, DriverParameters],
    throat_radius: float,
    target: TargetSpec,
    output_dir: str,
    top_n: int = 5,
) -> Path:
    """Generate the sweep report with rankings, plots, and CSVs.

    Args:
        all_ranked: Combined ranked results from all profiles.
        solver_csvs: Mapping of horn_label -> solver CSV path.
        drivers: Mapping of driver_id -> DriverParameters.
        throat_radius: Horn throat radius in metres.
        target: Target frequency specification.
        output_dir: Directory for output files.
        top_n: Number of top candidates to include in detailed output.

    Returns:
        Path to the output directory.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Sort all results by composite score
    all_ranked.sort(key=lambda r: r["composite_score"], reverse=True)
    top_results = all_ranked[:top_n]

    # 1. Full ranking JSON
    (out / "sweep_ranking.json").write_text(json.dumps(all_ranked, indent=2))

    # 2. Generate individual coupled-SPL CSVs for top candidates
    csv_pairs = []  # (csv_path, label) for comparison plot
    for rank, result in enumerate(top_results, 1):
        driver_id = result["driver_id"]
        horn_label = result["horn_label"]

        if driver_id not in drivers or horn_label not in solver_csvs:
            continue

        drv = drivers[driver_id]
        solver_csv = solver_csvs[horn_label]
        df = pd.read_csv(solver_csv)
        freq = df["frequency"].values
        solver_spl = df["spl"].values
        z_real = df["z_real"].values
        z_imag = df["z_imag"].values
        throat_area = np.pi * throat_radius ** 2

        p_throat = compute_driver_response(drv, freq, z_real, z_imag, throat_area)
        coupled_spl = scale_solver_spl(solver_spl, p_throat)

        csv_name = f"coupled_{rank:02d}_{driver_id}_{horn_label}.csv"
        csv_path = out / csv_name
        pd.DataFrame({
            "frequency": freq,
            "spl": coupled_spl,
        }).to_csv(csv_path, index=False)

        label = f"#{rank} {drv.manufacturer} {drv.model_name} ({horn_label})"
        csv_pairs.append((str(csv_path), label))

    # 3. Comparison plot
    if csv_pairs:
        plot_multi_comparison(
            csv_pairs,
            str(out / "sweep_comparison.png"),
            kpi_table=True,
        )

    # 4. Human-readable summary
    lines = [
        "Horn Driver Sweep Results",
        "=" * 40,
        f"Target: {target.f_low_hz:.0f} Hz - {target.f_high_hz:.0f} Hz",
        f"Throat radius: {throat_radius:.4f} m",
        f"Profiles evaluated: {', '.join(solver_csvs.keys())}",
        f"Total candidates scored: {len(all_ranked)}",
        "",
        f"Top {len(top_results)} Results:",
        "-" * 40,
    ]

    for rank, result in enumerate(top_results, 1):
        lines.append(
            f"  {rank}. {result.get('manufacturer', '')} {result.get('model_name', '')} "
            f"+ {result['horn_label']}"
        )
        lines.append(f"     Score: {result['composite_score']:.3f}  "
                      f"BW coverage: {result['bandwidth_coverage']:.1%}  "
                      f"Ripple: {result['passband_ripple_db']:.1f} dB  "
                      f"Sensitivity: {result['avg_sensitivity_db']:.1f} dB")
        if "kpi" in result:
            kpi = result["kpi"]
            f3l = f"{kpi['f3_low_hz']:.0f}" if kpi.get("f3_low_hz") else "N/A"
            f3h = f"{kpi['f3_high_hz']:.0f}" if kpi.get("f3_high_hz") else "N/A"
            lines.append(f"     f3: {f3l} - {f3h} Hz  "
                          f"Peak: {kpi['peak_spl_db']:.1f} dB @ {kpi['peak_frequency_hz']:.0f} Hz")
        lines.append("")

    (out / "sweep_summary.txt").write_text("\n".join(lines))

    print(f"Sweep report generated in {out}/")
    print(f"  - sweep_ranking.json ({len(all_ranked)} entries)")
    print(f"  - sweep_comparison.png (top {len(csv_pairs)} overlaid)")
    print(f"  - sweep_summary.txt")
    print(f"  - {len(csv_pairs)} individual coupled-SPL CSVs")

    return out


def main():
    """CLI for generating a sweep report from pre-computed rankings."""
    parser = argparse.ArgumentParser(
        description="Generate sweep report from ranked results.",
    )
    parser.add_argument("--ranked-json", required=True, help="Ranked results JSON file.")
    parser.add_argument("--solver-csvs", nargs="+", required=True,
                        help="Solver CSVs as profile:path pairs (e.g. conical:results.csv).")
    parser.add_argument("--drivers-db", required=True, help="Driver database JSON.")
    parser.add_argument("--throat-radius", type=float, required=True, help="Throat radius (m).")
    parser.add_argument("--target-f-low", type=float, required=True, help="Target low freq (Hz).")
    parser.add_argument("--target-f-high", type=float, required=True, help="Target high freq (Hz).")
    parser.add_argument("--top-n", type=int, default=5, help="Number of top results for detailed output.")
    parser.add_argument("--output-dir", type=str, default="sweep_report", help="Output directory.")
    args = parser.parse_args()

    from horn_drivers.loader import load_drivers

    all_ranked = json.loads(Path(args.ranked_json).read_text())

    solver_csvs = {}
    for pair in args.solver_csvs:
        profile, path = pair.split(":", 1)
        solver_csvs[profile] = path

    driver_list = load_drivers(args.drivers_db)
    drivers = {d.driver_id: d for d in driver_list}

    target = TargetSpec(f_low_hz=args.target_f_low, f_high_hz=args.target_f_high)

    generate_sweep_report(
        all_ranked=all_ranked,
        solver_csvs=solver_csvs,
        drivers=drivers,
        throat_radius=args.throat_radius,
        target=target,
        output_dir=args.output_dir,
        top_n=args.top_n,
    )


if __name__ == "__main__":
    main()
