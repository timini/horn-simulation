#!/usr/bin/env python3
"""Render a completed reference audit without rerunning or changing the models."""
import argparse
import base64
import hashlib
from html import escape
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from validate_reference_audit import interpolate

CASE_NAMES = {"Brass_O": "Brass cylinder, open", "Brass_C": "Brass cylinder, closed",
              "Wood_O": "Wood cylinder, open", "Wood_C": "Wood cylinder, closed",
              "3D_O": "ABS cylinder, open", "3D_C": "ABS cylinder, closed",
              "Cone_O": "ABS cone, open", "Cone_C": "ABS cone, closed"}


def table(rows):
    return pd.DataFrame(rows).to_html(index=False, border=0, escape=True, float_format=lambda x: f"{x:.3g}")


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_prediction_hashes(root, audit):
    if file_sha(root/"prediction_hashes.json") != audit["prediction_hashes_sha256"]:
        raise ValueError("Changed prediction hash manifest")
    for filename, digest in json.loads((root/"prediction_hashes.json").read_text()).items():
        if file_sha(root/filename) != digest:
            raise ValueError(f"Changed prediction: {filename}")


def archived_fem_impedance(frame, air):
    rho, c = air["rho"], air["c"]
    if not np.isfinite([rho, c]).all() or min(rho, c) <= 0:
        raise ValueError("Invalid frozen air constants")
    if not (np.allclose(frame.air_rho_kg_m3, rho, rtol=1e-12, atol=0)
            and np.allclose(frame.air_c_m_s, c, rtol=1e-12, atol=0)):
        raise ValueError("Solver CSV disagrees with frozen air constants")
    return (frame.z_real.to_numpy()+1j*frame.z_imag.to_numpy())/(rho*c)


def verify_import_inventory(path):
    root = path.parent
    inventory = json.loads(path.read_text())
    sources, dataset_ids = set(), set()
    curves = diy_curves = 0
    for dataset in inventory["datasets"]:
        reference = dataset["reference"]
        if reference["id"] in dataset_ids:
            raise ValueError("Duplicate dataset in inventory")
        dataset_ids.add(reference["id"])
        if file_sha(root/reference["archive_name"]) != reference["sha256"]:
            raise ValueError("Changed imported archive")
        directory = root/reference["id"]
        if json.loads((directory/"manifest.json").read_text()) != dataset:
            raise ValueError("Dataset manifest disagrees with import inventory")
        for curve in dataset["curves"]:
            if file_sha(directory/curve["csv"]) != curve["csv_sha256"]:
                raise ValueError("Changed imported curve")
            curves += 1
            diy_curves += reference["format"] == "zip_frd_zma"
            sources.add(curve["source_sha256"])
    if len(sources) != inventory["unique_curve_files"]:
        raise ValueError("Incorrect unique-curve inventory count")
    return {"archives": len(dataset_ids), "curves": curves, "unique_curve_files": len(sources),
            "diy_curves": diy_curves, "inventory_sha256": file_sha(path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--inventory-json", type=Path, required=True)
    parser.add_argument("--search-json", type=Path, required=True)
    parser.add_argument("--resonance-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.audit_dir, args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    audit = json.loads((root/"audit.json").read_text())
    search = json.loads(args.search_json.read_text())
    resonance = json.loads(args.resonance_json.read_text())
    manifest = json.loads((args.reference_dir/"manifest.json").read_text())
    protocol = json.loads((root/"protocol.json").read_text())
    imports = verify_import_inventory(args.inventory_json)
    verify_prediction_hashes(root, audit)
    if hashlib.sha256((root/"protocol.json").read_bytes()).hexdigest() != audit["protocol_sha256"]:
        raise ValueError("Changed protocol")
    if hashlib.sha256((args.reference_dir/"manifest.json").read_bytes()).hexdigest() != protocol["manifest_sha256"]:
        raise ValueError("Changed manifest")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(4, 2, figsize=(13, 15), sharex=True)
    cases_in_order = [c for c in protocol["cases"] if c in audit["by_case"]]
    if len(cases_in_order) != 8:
        raise ValueError("This report layout requires eight measured configurations")
    for ax, case in zip(axes.flat, cases_in_order):
        p = pd.read_csv(root/f"{case}-prediction.csv")
        frequency = p.frequency.to_numpy()
        measurements = []
        for c in manifest["curves"]:
            if c["reference_kind"] != "measurement" or c["configuration"] != case:
                continue
            path = args.reference_dir/c["csv"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != c["csv_sha256"]:
                raise ValueError("Changed measurement")
            z = interpolate(pd.read_csv(path), frequency, "z_real_normalized", "z_imag_normalized")
            measurements.append(20*np.log10(np.maximum(np.abs(z), 1e-15)))
        low, median, high = np.percentile(measurements, [10, 50, 90], axis=0)
        ax.fill_between(frequency, low, high, color="#8996a5", alpha=.25, label="Measured 10–90% range")
        ax.semilogx(frequency, median, color="#344457", lw=1.4, label="Measured median")
        z = p.z_real_normalized.to_numpy()+1j*p.z_imag_normalized.to_numpy()
        ax.semilogx(frequency, 20*np.log10(np.abs(z)), color="#007e91", lw=1.2, label="TMM: 20,001 points")
        fem = pd.read_csv(root/f"{protocol['model_aliases'][case]}-h0.004-n481.csv")
        ax.semilogx(fem.frequency, 20*np.log10(np.abs(archived_fem_impedance(fem, protocol["air"]))), ":", color="#d96625", lw=1.2, label="Production FEM: 481 points")
        counts = audit["by_case"][case]
        ax.set_title(f"{CASE_NAMES[case]}  ·  TMM {counts['tmm_dense']['magnitude_passes']}/{counts['curves']}, FEM {counts['fem_481']['magnitude_passes']}/{counts['curves']}", loc="left")
        ax.set_ylabel("Normalized impedance magnitude (dB)")
        ax.set_xlim(110, 3900)
        ax.set_xticks([110, 300, 1000, 3000], labels=["110", "300", "1,000", "3,000"])
        ax.grid(alpha=.15)
    for ax in axes[-1]:
        ax.set_xlabel("Frequency (Hz)")
    axes[0, 0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Measured pipe impedance versus unfitted predictions\n110–3900 Hz · nominal dimensions · dry-air source constants", fontsize=15)
    fig.text(.5, .005, "Measurements: Ernoult et al., Zenodo 20024938 (CC BY 4.0). Ranges describe repeated files, not confidence intervals. Counts are magnitude agreement only.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .02, 1, .96))
    fig.savefig(out/"reference-comparison.png", dpi=160)
    plt.close(fig)
    totals = audit["totals"]
    cases = [{"Case": CASE_NAMES[c], "Curves": d["curves"], "TMM matches": d["tmm_dense"]["magnitude_passes"],
        "FEM sampled matches": d["fem_481"]["magnitude_passes"],
        "TMM typical p95 error (dB)": d["tmm_dense"]["median_p95_db"],
        "FEM typical p95 error (dB)": d["fem_481"]["median_p95_db"],
        "FEM typical p95 phase error (°)": d["fem_481"]["median_p95_phase_error_deg"]}
        for c, d in audit["by_case"].items()]
    operators = [{"Experimenter": c, "Curves": d["curves"],
        "TMM matches": d["tmm_dense"]["magnitude_passes"],
        "FEM sampled matches": d["fem_481"]["magnitude_passes"]} for c, d in audit["by_operator"].items()]
    convergence = [{"Case": c, "Maximum mesh change (dB)": d["middle_vs_fine"]["max_db"],
        "95th-percentile scaled complex change": d["middle_vs_fine"]["p95_scaled_complex_error"],
        "Meets frozen mesh limits": d["passed"]} for c, d in audit["mesh_convergence"].items()]
    sampling = [{"Frequency points": n, **d} for n, d in audit["frequency_sampling"].items()]
    same_grid_changes = sum(r["fem_481"]["magnitude_agreement"] != r["sampling"]["481"]["magnitude_agreement"]
                            for r in audit["measurements"])
    numerical = [{"Reference file": r["source_member"].split("/")[-1],
        "Samples": r["sample_count"], "Maximum gap (Hz)": r["maximum_frequency_gap_hz"],
        "p95 error (dB)": r["tmm_at_reference_frequencies"]["p95_db"],
        "max error (dB)": r["tmm_at_reference_frequencies"]["max_db"],
        "p95 phase error (°)": r["tmm_at_reference_frequencies"]["p95_phase_error_deg"],
        "Sampled magnitude match": r["tmm_at_reference_frequencies"]["magnitude_agreement"]}
        for r in audit["independent_simulations"]]
    searches = [{"Band (Hz)": f"{c['target']['f_low_hz']:g}–{c['target']['f_high_hz']:g}",
        "Feasible pairs": c["audit"]["feasible_pair_count"],
        "Top-ten recall": c["audit"]["feasible_top_k_recall"],
        "Winner retained": c["audit"]["winner_retained"], "Score regret": c["audit"]["score_regret"]} for c in search["cases"]]
    measured = [{"File": r["source_member"].split("/")[-1], "Case": r["configuration"],
        "Experimenter": r["operator"], "TMM p95 (dB)": r["tmm_dense"]["p95_db"],
        "FEM p95 (dB)": r["fem_481"]["p95_db"], "FEM max (dB)": r["fem_481"]["max_db"],
        "FEM p95 phase (°)": r["fem_481"]["p95_phase_error_deg"],
        "FEM magnitude match": r["fem_481"]["magnitude_agreement"]} for r in audit["measurements"]]
    picture = base64.b64encode((out/"reference-comparison.png").read_bytes()).decode()
    search_heading = "The search passed its benchmark." if search["passed"] else "The search audit failed."
    html = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Horn reference validation — 8 September 2026</title><style>
body{{font:16px/1.55 system-ui,sans-serif;color:#213142;background:#f4f6f8;margin:0}}main{{max-width:1200px;margin:auto;padding:36px}}h1{{font-size:36px;line-height:1.15}}h2{{margin-top:36px}}p{{max-width:1000px}}.verdict{{background:#fff2d9;border-left:5px solid #bc7310;padding:18px 22px}}.cards{{display:flex;flex-wrap:wrap;gap:16px;margin:24px 0}}.card{{background:white;padding:20px;flex:1;min-width:170px;border-radius:8px}}.card b{{font-size:30px;display:block}}table{{border-collapse:collapse;width:100%;font-size:13px;background:white;margin:15px 0}}th,td{{text-align:left;padding:9px;border-bottom:1px solid #dce2e7;overflow-wrap:anywhere}}th{{background:#e7eef3}}.scroll{{overflow:auto}}img{{width:100%;height:auto}}summary{{cursor:pointer;font-weight:650;padding:12px 0}}a{{color:#006979}}code{{overflow-wrap:anywhere}}small{{color:#586b7b}}
</style><main><small>REFERENCE AUDIT · 8 SEPTEMBER 2026</small><h1>{search_heading}<br>The acoustic model is only partly validated.</h1>
<div class="verdict">This fresh audit checks every acquired pipe reference and repeats the finite-grid optimiser benchmark. It does <strong>not</strong> certify complete horn-and-driver response, absolute SPL, off-axis behaviour or physical recommendation ordering. Failed comparisons remain failures.</div>
<div class="cards"><div class="card"><b>{totals['tmm_magnitude_passes']}/299</b>measured curves within dense TMM magnitude limits</div><div class="card"><b>{totals['fem_sampled_magnitude_passes']}/299</b>within production FEM sampled magnitude limits</div><div class="card"><b>{totals['independent_simulation_sampled_magnitude_matches']}/40</b>independent simulations matching at supplied frequencies</div><div class="card"><b>{sum(c['audit']['passed'] for c in search['cases'])}/4</b>optimiser audit bands passed</div></div>
<h2>What was compared</h2><p><a href="https://zenodo.org/records/20024938">Ernoult and colleagues’ reference archive</a> provides 299 measured and 40 independently simulated impedance curves. Eight measured configurations reduce to five distinct nominal rigid-wall air-domain models; ABS and wood compliance are not modeled. All original curve hashes were checked. The verified import inventory contains {imports['curves']} curves ({imports['unique_curve_files']} unique files) from {imports['archives']} checksum-verified archives. Its {imports['diy_curves']} DIY response/electrical-impedance curves remain unmatched to a fully specified supported assembly.</p>
<p>Predictions use the source’s dry-air 25°C constants, the published dimensions, side-wall losses and finite-flange or rigid closed terminations. No damping, level, frequency or geometry parameters were fitted. Magnitude limits remain median absolute error ≤2 dB and 95th-percentile error ≤4 dB across 110–3900 Hz. These are project tolerances, not the paper’s statistical uncertainty test. Phase errors are reported separately, without a phase acceptance gate.</p>
<div class="scroll">{table(cases)}</div><img alt="Eight comparisons showing measured spread, median, dense TMM and production FEM" src="data:image/png;base64,{picture}">
<h2>Numerical reliability and frequency resolution</h2><p>{totals['healthy_fem_runs']}/15 production FEM runs passed residual, solver convergence, nonnegative dissipation, passive input impedance and energy-balance checks. {totals['mesh_converged_cases']}/5 geometries meet the predeclared mesh-change limits: ≤0.5 dB maximum impedance change and ≤5% scaled complex error at the 95th percentile between 4 mm and 3 mm meshes. The mesh comparison uses 121 common frequencies, so it is not a continuous-band guarantee.</p>{table(convergence)}
<p>The following test evaluates the same analytical model and measurements on progressively denser grids. Classification changes expose how a frequency grid can affect a pass/fail claim. FEM’s 481-point count is explicitly a sampled result, not equivalent to the dense 20,001-point TMM assessment.</p>{table(sampling)}
<p>On the same 481-point grid, FEM and TMM differ in {same_grid_changes} measured-curve classifications. A larger sampled match count must not be interpreted as improved model accuracy.</p>
<h2>Differences between experimenters</h2>{table(operators)}<p>These totals mix specimens and configurations; they are descriptive and do not rank experimenter accuracy. Repeated files are not independent assemblies. The historical held-out split is preserved in the raw audit, but those data had already been examined before this rerun.</p>
<p>The <a href="https://doi.org/10.1051/aacus/2026048">authors’ benchmark paper, section 5.3</a> also finds larger discrepancies for wooden pipes and closed ends. It discusses roughness, porosity and capping procedures as possible contributors. Those explanations do not establish the cause of each failure in our run.</p>
<h2>Independent horn datum</h2><p><a href="{escape(resonance['source_url'])}">Sa and Park (2014), section 2.2 / Fig. 13</a> reports the sixth measured impedance resonance of a 535 mm conical horn (18 mm throat, 80 mm mouth) at 1,712 Hz. Our unfitted model predicts {resonance['predictions']['800']['resonances_hz'][5]:,.1f} Hz, a {resonance['predictions']['800']['sixth_resonance_error_hz']:.1f} Hz difference; both 400 and 800 segments give this value. Air properties and lip geometry are assumptions, measurement uncertainty is not established, and this one scalar does not validate magnitude or driver coupling. A separate short-horn candidate has conflicting mouth diameters in the paper’s captions (80 versus 90 mm); it is not treated as a matched assembly.</p>
<h2>Optimiser against exhaustive simulation</h2><p>Freshly evaluated 20 horn geometries × 3 manufacturer-parameter drivers across four bands. Runtime: {search['total_seconds']/60:.2f} minutes. This tests whether screening loses good simulated choices; it cannot establish the best physical driver or a global optimum.</p>{table(searches)}
<h2>What still prevents full physical validation</h2><p>The measured pipe results apply to the opt-in boundary-layer loss configuration used here. The pipeline’s default lossless configuration remains experimental; it does not inherit these results.</p><ol><li>A matched complete assembly with characterized chamber, phase plug, adapter and rear load, plus electrical impedance and calibrated output measurements.</li><li>Resolved treatment of the failed material/closed-end cases and resonance errors, with independently specified uncertainties.</li><li>An appropriate exterior-radiation model for the actual mounting and observation geometry.</li><li>Measured comparisons between competing assemblies to validate recommendation order.</li></ol>
<details><summary>All 40 independent numerical comparisons</summary><p>These include different propagation, loss and wavefront assumptions. Comparisons use only the supplied frequencies. Two 3D references have only 3 and 11 points inside the band; they cannot resolve full response curves or resonances. The initial interpolated comparison was invalid and is preserved in the earlier diagnostic output. A sampled magnitude match is diagnostic and does not imply identical governing models or continuous-band agreement.</p><div class="scroll">{table(numerical)}</div></details>
<details><summary>All 299 measured comparisons, including failures</summary><div class="scroll">{table(measured)}</div></details>
<details><summary>Provenance and reproduction</summary><p>Base revision: <code>{escape(str(protocol['git_revision']))}</code>. Actual source-file hashes, dimensions and tolerances are stored in the frozen protocol. All predicted CSVs are hashed before the runner opens reference values.</p><p>Archive SHA-256: <code>{audit['source_archive_sha256']}</code>.</p><p>Full raw audit: <code>{escape(str(root/'audit.json'))}</code>. Runner: <code>scripts/validate_reference_audit.py</code>. Protocol SHA-256: <code>{audit['protocol_sha256']}</code>.</p></details></main></html>"""
    (out/"reference-report.html").write_text(html)
    compact = {k: v for k, v in audit.items() if k not in ("measurements", "independent_simulations")}
    compact["fem_vs_tmm_same_grid_classification_changes"] = same_grid_changes
    compact["independent_simulation_error_range"] = {
        "p95_db": [min(r["tmm_at_reference_frequencies"]["p95_db"] for r in audit["independent_simulations"]),
                   max(r["tmm_at_reference_frequencies"]["p95_db"] for r in audit["independent_simulations"])],
        "max_db": max(r["tmm_at_reference_frequencies"]["max_db"] for r in audit["independent_simulations"])}
    compact["verified_import_inventory"] = imports
    compact["search"] = {"passed": search["passed"], "total_seconds": search["total_seconds"], "cases": searches}
    compact["published_horn_scalar"] = resonance
    compact["artifact_sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (root/"audit.json", args.search_json, args.resonance_json, out/"reference-comparison.png")}
    (out/"summary.json").write_text(json.dumps(compact, indent=2, allow_nan=False)+"\n")
    print(out/"reference-report.html")


if __name__ == "__main__":
    main()
