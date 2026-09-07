"""Reports for experimental driver-horn predictions and their evidence gaps."""
from horn_analysis.evaluation import coupled_output

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_analysis.compare import plot_multi_comparison
from horn_analysis.html_report import generate_html_report
from horn_analysis.scoring import TargetSpec
from horn_analysis.transfer_function import compute_driver_response, scale_solver_spl


def generate_auto_report(
    all_ranked: List[dict],
    solver_csvs: Dict[str, str],
    drivers: Dict[str, DriverParameters],
    throat_radius: float | None,
    target: TargetSpec,
    output_dir: str,
    top_n: int = 5,
    mouth_radius: float | None = None,
    horn_length: float | None = None,
    derived_geometry: Optional[dict] = None,
    total_candidates: int | None = None,
    total_scored: int | None = None,
    lem_results: Optional[dict] = None,
    no_feasible_reason: Optional[str] = None,
) -> Path:
    """Generate the auto-select report with rankings, plots, and CSVs.

    Args:
        all_ranked: Combined ranked results from all profiles.
        solver_csvs: Mapping of horn_label -> solver CSV path.
        drivers: Mapping of driver_id -> DriverParameters.
        throat_radius: Horn throat radius in metres.
        target: Target frequency specification.
        output_dir: Directory for output files.
        top_n: Number of top candidates to include in detailed output.
        mouth_radius: Horn mouth radius in metres (for report display).
        horn_length: Horn length in metres (for report display).
        derived_geometry: Optional dict from geometry_designer (fullauto mode).
        total_candidates: Total number of geometry candidates simulated.
        total_scored: Total number of driver-horn combinations scored.

    Returns:
        Path to the output directory.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Sort all results by composite score
    all_ranked.sort(key=lambda r: r["composite_score"], reverse=True)
    from horn_analysis.search import annotate_comparable_candidates
    all_ranked = annotate_comparable_candidates(all_ranked)
    top_results = all_ranked[:top_n]

    # 1. Full ranking JSON
    (out / "auto_ranking.json").write_text(json.dumps(all_ranked, indent=2))

    # 2. Generate individual coupled-SPL CSVs for top candidates
    csv_pairs = []  # (csv_path, label) for comparison plot
    for rank, result in enumerate(top_results, 1):
        driver_id = result["driver_id"]
        horn_label = result["horn_label"]

        if driver_id not in drivers or horn_label not in solver_csvs:
            raise ValueError(f"Missing driver or solver response for {driver_id}/{horn_label}")

        drv = drivers[driver_id]
        solver_csv = solver_csvs[horn_label]
        df = pd.read_csv(solver_csv)
        freq = df["frequency"].values
        solver_spl = df["spl"].values
        z_real = df["z_real"].values
        z_imag = df["z_imag"].values
        cand_throat = result.get("throat_radius", throat_radius)

        coupled_spl, _, metric = coupled_output(df, drv, target, cand_throat)

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
            str(out / "auto_comparison.png"),
            kpi_table=True,
        )

    else:
        from horn_analysis import plot_theme
        fig, ax = plot_theme.create_figure()
        ax.text(.5, .5, "No feasible design in evaluated candidates", ha="center", va="center", transform=ax.transAxes)
        fig.savefig(out / "auto_comparison.png")
        import matplotlib.pyplot as plt
        plt.close(fig)

    # 4. Human-readable summary
    scored_display = total_scored if total_scored is not None else len(all_ranked)
    lines = [
        "Horn Driver Auto-Select Results — experimental predictions",
        "Status: experimental_candidates" if all_ranked else "Status: no_feasible_design",
        "No recommendation is physically validated. See evidence gaps in ranking JSON.",
        "Scores within 0.02 are near-ties for comparison; physical uncertainty is not quantified.",
        "=" * 40,
        f"Target: {target.f_low_hz:.0f} Hz - {target.f_high_hz:.0f} Hz",
        "Throat radii: see each candidate below" if derived_geometry else (f"Throat radius: {throat_radius:.4f} m" if throat_radius is not None else "Throat radius: not available"),
        f"Profiles evaluated: {', '.join(solver_csvs.keys())}",
        f"Total candidates scored: {scored_display}",
    ]

    if lem_results:
        lines.extend([
            "",
            "LEM Prescreening:",
            f"  Candidates evaluated by LEM: {lem_results.get('total_evaluated', 'N/A')}",
            f"  Driver-horn pairs scored: {lem_results.get('total_pairs', 'N/A')}",
            f"  Passed to FEM: {len(lem_results.get('filtered_candidate_ids', []))}",
        ])

    lines.extend([
        "",
        f"Top {len(top_results)} Results:",
        "-" * 40,
    ])

    for rank, result in enumerate(top_results, 1):
        lines.append(
            f"  {rank}. {result.get('manufacturer', '')} {result.get('model_name', '')} "
            f"+ {result['horn_label']}"
        )
        lines.append(f"     Score: {result['composite_score']:.3f}  "
                      f"BW coverage: {result['bandwidth_coverage']:.1%}  "
                      f"Ripple: {result['passband_ripple_db']:.1f} dB  "
                      f"Mean output level: {result['avg_sensitivity_db']:.1f} dB")
        if all(result.get(key) is not None for key in ("throat_radius", "mouth_radius", "length")):
            lines.append(f"     Throat radius: {result['throat_radius']:.6f} m; mouth radius: {result['mouth_radius']:.6f} m; length: {result['length']:.6f} m")
        if "kpi" in result:
            kpi = result["kpi"]
            f3l = f"{kpi['f3_low_hz']:.0f}" if kpi.get("f3_low_hz") else "N/A"
            f3h = f"{kpi['f3_high_hz']:.0f}" if kpi.get("f3_high_hz") else "N/A"
            lines.append(f"     f3: {f3l} - {f3h} Hz  "
                          f"Peak: {kpi['peak_spl_db']:.1f} dB @ {kpi['peak_frequency_hz']:.0f} Hz")
        lines.append("")

    if no_feasible_reason and not all_ranked:
        lines.append("Reason: " + no_feasible_reason)
    (out / "auto_summary.txt").write_text("\n".join(lines))

    # 5. Self-contained HTML report
    html_report = generate_html_report(
        all_ranked=all_ranked,
        solver_csvs=solver_csvs,
        drivers=drivers,
        throat_radius=throat_radius,
        target=target,
        csv_pairs=csv_pairs,
        top_n=top_n,
        mouth_radius=mouth_radius,
        length=horn_length,
        derived_geometry=derived_geometry,
        total_candidates=total_candidates,
        total_scored=total_scored,
        lem_results=lem_results,
        no_feasible_reason=no_feasible_reason,
    )
    (out / "auto_report.html").write_text(html_report)

    print(f"Auto-select report generated in {out}/")
    print(f"  - auto_ranking.json ({len(all_ranked)} entries)")
    print(f"  - auto_comparison.png (top {len(csv_pairs)} overlaid)")
    print(f"  - auto_summary.txt")
    print(f"  - auto_report.html")
    print(f"  - {len(csv_pairs)} individual coupled-SPL CSVs")

    return out


def main():
    """CLI for generating an auto-select report from pre-computed rankings."""
    parser = argparse.ArgumentParser(
        description="Generate auto-select report from ranked results.",
    )
    parser.add_argument("--ranked-json", required=True, help="Ranked results JSON file.")
    parser.add_argument("--solver-csvs", nargs="+", required=True,
                        help="Solver CSVs as profile:path pairs (e.g. conical:results.csv).")
    parser.add_argument("--drivers-db", required=True, help="Driver database JSON.")
    parser.add_argument("--throat-radius", type=float, required=True, help="Throat radius (m).")
    parser.add_argument("--target-f-low", type=float, required=True, help="Target low freq (Hz).")
    parser.add_argument("--target-f-high", type=float, required=True, help="Target high freq (Hz).")
    parser.add_argument("--top-n", type=int, default=5, help="Number of top results for detailed output.")
    parser.add_argument("--output-dir", type=str, default="auto_report", help="Output directory.")
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

    generate_auto_report(
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
