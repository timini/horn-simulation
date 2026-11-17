"""Reports for experimental driver-horn predictions and their evidence gaps."""

from horn_analysis.kpi import format_kpi

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_analysis.compare import plot_multi_comparison
from horn_analysis.html_report import generate_html_report
from horn_analysis.scoring import TargetSpec, size_compromise_note
from horn_analysis.electrical_output import electrical_table, plot_electrical


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
    scoring: Optional[dict] = None,
    rejected: Optional[List[dict]] = None,
    rejection_stage: Optional[str] = None,
    step_files: Optional[Dict[str, str]] = None,
    pricing: Optional[dict] = None,
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
        scoring: Score formula, weights and bound origins (``describe_scoring``).
            Derived from ``target`` when not given.
        rejected: Rows that failed a hard limit. Infeasible rows in
            ``all_ranked`` are moved here. They are summarised by reason with
            margins in the report, auto_summary.txt and auto_rejections.json.
        rejection_stage: Where those rows were rejected, e.g. "FEM" or
            "analytical prescreen".
        step_files: horn_label -> STEP path to show for each ranked horn.
        pricing: Run pricing summary from ranked_results.json; rows carry
            their own ``price``. Shown in the summary and HTML, never scored.

    Returns:
        Path to the output directory.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if scoring is None:
        from horn_analysis.scoring import describe_scoring
        scoring = describe_scoring(target, derived_geometry)

    from horn_analysis.report_tables import (candidate_label, horn_legend, chamber_assumption, has_parametric_geometry,
                                             profile_table, rejection_summary)
    rejected = list(rejected or []) + [r for r in all_ranked if r.get("model_feasible") is False]
    all_ranked = [r for r in all_ranked if r.get("model_feasible") is not False]
    step_files = step_files or {}
    if not all_ranked and rejected and not no_feasible_reason:
        no_feasible_reason = ("Every evaluated driver-horn pair failed at least one hard limit; "
                              "see the rejection table for the reasons and how close they came.")

    # Sort all results by composite score
    all_ranked.sort(key=lambda r: r["composite_score"], reverse=True)
    from horn_analysis.search import annotate_comparable_candidates
    all_ranked = annotate_comparable_candidates(all_ranked)
    top_results = all_ranked[:top_n]

    # 1. Full ranking JSON
    (out / "auto_ranking.json").write_text(json.dumps(all_ranked, indent=2))

    # 2. Generate individual coupled-SPL CSVs for top candidates
    csv_pairs = []  # (csv_path, label) for comparison plot
    electrical_tables = []
    builder_files = {}
    for rank, result in enumerate(top_results, 1):
        driver_id = result["driver_id"]
        horn_label = result["horn_label"]

        if driver_id not in drivers or horn_label not in solver_csvs:
            raise ValueError(f"Missing driver or solver response for {driver_id}/{horn_label}")

        drv = drivers[driver_id]
        solver_csv = solver_csvs[horn_label]
        df = pd.read_csv(solver_csv)
        cand_throat = result.get("throat_radius", throat_radius)

        # frequency, spl, Re/Im/|Ze| (ohm), excursion (mm), input power and
        # system group delay; see horn_analysis.electrical_output.
        table = electrical_table(df, drv, target, cand_throat, result.get("length") or None)

        csv_name = f"coupled_{rank:02d}_{driver_id}_{horn_label}.csv"
        csv_path = out / csv_name
        table.to_csv(csv_path, index=False)

        label = candidate_label(rank, result)
        csv_pairs.append((str(csv_path), label))
        electrical_tables.append((label, table, drv.xmax_m))

        if horn_label not in builder_files:
            files = {"step": step_files.get(horn_label)}
            if has_parametric_geometry(result):
                profile_name = f"profile_{rank:02d}_{horn_label}.csv"
                profile_table(result).to_csv(out / profile_name, index=False)
                files["profile_csv"] = profile_name
            builder_files[horn_label] = files

    summary = rejection_summary(rejected, target, drivers)
    (out / "auto_rejections.json").write_text(json.dumps(
        {"stage": rejection_stage, "rejected_pairs": len(rejected), "reasons": summary}, indent=2))

    if electrical_tables:
        plot_electrical(electrical_tables, str(out / "auto_electrical.png"))

    # Experimental infinite-baffle directivity for each modal horn shown.
    # The pattern is relative to the axis and so independent of the driver.
    from horn_analysis.directivity import directivity_table, beamwidth_figure, DIRECTIVITY_LABEL
    directivity_tables = []
    seen_horns = set()
    for rank, result in enumerate(top_results, 1):
        horn_label = result["horn_label"]
        if horn_label in seen_horns:
            continue
        seen_horns.add(horn_label)
        table = directivity_table(pd.read_csv(solver_csvs[horn_label]))
        if table is not None:
            table.to_csv(out / f"directivity_{horn_label}.csv", index=False)
            # Same '#rank horn' legend as the other horn-only plots.
            directivity_tables.append((horn_legend(rank, result), table))
    if directivity_tables:
        from horn_analysis import plot_theme
        plot_theme.save_figure(beamwidth_figure(directivity_tables, target), str(out / "auto_beamwidth.png"))

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
        "Horn Driver Auto-Select Results: EXPERIMENTAL PREDICTION",
        "Status: experimental_candidates" if all_ranked else "Status: no_feasible_design",
        "No recommendation is physically validated. See evidence gaps in ranking JSON.",
        "Candidates are named '#rank horn_label + driver_id'; plot legends and CSV names use the same IDs.",
        "Dimensions are in mm. Profile CSVs give the air-passage wall radius r_mm against axial position z_mm (z = 0 at the throat).",
        "f3 is 3 dB below the target-band mean level, searched outward from the band edges over the full sweep; ≤ / ≥ values are sweep bounds.",
        "Scores within 0.02 are near-ties for comparison; physical uncertainty is not quantified.",
        "=" * 40,
        f"Target: {target.f_low_hz:.0f} Hz - {target.f_high_hz:.0f} Hz",
        "Throat radii: see each candidate below" if derived_geometry else (f"Throat radius: {throat_radius:.4f} m" if throat_radius is not None else "Throat radius: not available"),
        f"Horns evaluated: {', '.join(solver_csvs.keys())}",
        f"Total candidates scored: {scored_display}",
    ]

    weights = scoring["weights"]
    lines.extend(["", "Score: " + scoring["formula"],
                  "Weights: " + ", ".join(f"{k} {v:.3f}" for k, v in weights.items()) + f" ({scoring['weights_source']})"])
    lines.extend(f"  {k.replace('_', ' ')}: {v:.2f}" for k, v in scoring["exchange_rates"].items() if isinstance(v, float))
    lines.extend(f"Bound {k}: {v['source']}" + (f" {v['searched_range']}" if v.get("searched_range") else "")
                 for k, v in scoring["geometry_bounds"].items())
    lines.extend(f"Constraint {k}: {v['value']} ({v['source']})" for k, v in scoring["constraints"].items())
    compromise = scoring.get("size_compromise")
    if compromise:
        lines.append("Size compromise: " + size_compromise_note(compromise, all_ranked))
    refinement = scoring.get("refinement")
    if refinement:
        lines.append(f"Refinement: {refinement['new_evaluations']} new geometries from "
                     f"{refinement['seed_candidate']} ({refinement['seed_source'].replace('_', ' ')})")
        lines.extend(f"  Widened heuristic {w['field']} {w['side']} bound: {w['from']} to {w['to']}"
                     for w in refinement.get("widened_bounds", []))
    if scoring.get("front_chamber_policy"):
        lines.append("Front chamber: " + scoring["front_chamber_policy"])
    from horn_analysis.pricing import describe_pricing, format_price
    lines.extend(describe_pricing(pricing))

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
            f"  {candidate_label(rank, result)} ({result.get('manufacturer', '')} {result.get('model_name', '')}); "
            f"evidence status: {result.get('eligibility_status', 'unknown')}"
        )
        lines.append(f"     Score: {result['composite_score']:.3f}  "
                      f"Extension: {result['extension_octaves']:+.2f} oct  "
                      f"Ripple: {result['passband_ripple_db']:.1f} dB  "
                      f"Mean output level: {result['avg_sensitivity_db']:.1f} dB")
        drv = drivers.get(result["driver_id"])
        if drv is not None:
            lines.append("     Driver rating: "
                         + (f"{drv.power_w:g} W nominal" if drv.power_w is not None else "nominal power not recorded")
                         + (f", {drv.peak_power_w:g} W peak" if drv.peak_power_w is not None else "")
                         + (f"; sensitivity {drv.sensitivity_db_1w1m:.1f} dB at 1 W, 1 m (datasheet)"
                            if drv.sensitivity_db_1w1m is not None else "; sensitivity not recorded"))
        if pricing:
            price = result.get("price")
            lines.append("     Price: " + format_price(price)
                         + (f" ({price['source_id']}, {price['retrieved_on']}, {price['url']})" if price else "")
                         + ("; cheapest feasible" if result.get("cheapest_feasible") else ""))
        if all(result.get(key) is not None for key in ("throat_radius", "mouth_radius", "length")):
            lines.append(f"     Throat radius: {result['throat_radius']*1e3:.2f} mm; mouth radius: {result['mouth_radius']*1e3:.2f} mm; length: {result['length']*1e3:.2f} mm")
        files = builder_files.get(result["horn_label"], {})
        if files.get("profile_csv"):
            lines.append(f"     Profile CSV: {files['profile_csv']}")
        lines.append(f"     STEP (relative to this folder): {files.get('step') or 'not available'}")
        interface = result.get("driver_interface")
        if interface:
            corner = interface.get("front_chamber_corner_hz")
            lines.append(f"     Front chamber: {chamber_assumption(result)}; "
                         f"Sd/St {interface['compression_ratio']:.2f}"
                         + (f"; lumped corner on a resistive load {corner:.0f} Hz" if corner else ""))
        if result.get("front_chamber_estimated"):
            only = result.get("estimated_chamber_reasons") or []
            ideal = result.get("without_front_chamber") or {}
            lines.append("     front_chamber_estimated: yes"
                         + (f"; fails only with the estimate: {', '.join(only)}" if only else "")
                         + (f"; ripple without a chamber {ideal['passband_ripple_db']:.1f} dB" if ideal.get("passband_ripple_db") is not None else ""))
        if result.get("test_power_w") is not None:
            margin = result.get("test_power_xmax_margin_mm")
            lines.append(f"     At {result['test_power_w']:g} W: max in-band excursion {result['test_power_band_max_excursion_mm']:.2f} mm"
                         + (f", Xmax margin {margin:.2f} mm" if margin is not None else ", Xmax unknown")
                         + (f"; max in-band input power {result['test_power_band_max_input_power_w']:.1f} W"
                            if result.get("test_power_band_max_input_power_w") is not None else ""))
        if "kpi" in result:
            kpi = result["kpi"]
            f3l = format_kpi(result, "f3_low_hz")
            f3h = format_kpi(result, "f3_high_hz")
            lines.append(f"     f3: {f3l} - {f3h} Hz  "
                          f"Peak: {kpi['peak_spl_db']:.1f} dB @ {kpi['peak_frequency_hz']:.0f} Hz")
        lines.append("")

    if any(r.get("radiation_model") == "flanged_piston" for r in top_results):
        from horn_analysis.html_report import FLANGED_PISTON_CAVEAT
        lines.append(FLANGED_PISTON_CAVEAT)
    if directivity_tables:
        lines.append(f"Directivity ({DIRECTIVITY_LABEL}): directivity_<horn>.csv and auto_beamwidth.png, "
                     "relative to the axis, far field, 0-90 degrees.")
    elif top_results:
        lines.append("Directivity: not available; it needs the modal_baffled radiation model.")
    if no_feasible_reason and not all_ranked:
        lines.append("Reason: " + no_feasible_reason)
    if summary:
        lines.extend(["", f"Rejected pairs: {len(rejected)}" + (f" ({rejection_stage})" if rejection_stage else ""),
                      "  reason: count, closest miss, largest excess"])
        for entry in summary:
            fmt = lambda v: "n/a" if v is None else f"{v:.2f} {entry['unit']}".strip()
            lines.append(f"  {entry['reason']}: {entry['count']}, {fmt(entry['closest_excess'])}, {fmt(entry['largest_excess'])}")
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
        scoring=scoring,
        directivity_tables=directivity_tables,
        electrical_tables=electrical_tables,
        electrical_png=str(out / "auto_electrical.png") if electrical_tables else None,
        rejected=rejected,
        rejection_stage=rejection_stage,
        builder_files=builder_files,
        pricing=pricing,
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
