"""Auto and single reports regenerated from preserved solver results, with no new solves.

Rankings are recomputed from the preserved FEM responses with the current
coupling and scoring code, then the report is generated. The two live runs
live outside the repository and are skipped when absent.
"""
import csv
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from horn_core.profiles import radius_profile
from horn_drivers.loader import load_driver
from horn_analysis.auto_report import generate_auto_report
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "examples/800-1600-modal"
RESULTS = Path("/Users/timrichardson/Documents/projects/personal/horn-simulation/results")
LIVE_RUNS = ["20261115T181353Z-78a53b73", "20261115T181058Z-d0f72ed3"]


def _example_inputs():
    preserved = json.loads((EXAMPLE / "ranking.json").read_text())[0]
    driver = load_driver(str(ROOT / "data/drivers-curated"), preserved["driver_id"])
    return preserved, driver


def _annotate(rows, geometry):
    for row in rows:
        row.update({k: geometry[k] for k in ("profile", "throat_radius", "mouth_radius", "length")})
    return rows


def _example_report(out, target=None, monkeypatch=None):
    preserved, driver = _example_inputs()
    target = target or TargetSpec(800, 1600)
    label = preserved["horn_label"]
    rows = _annotate(rank_horn_drivers(str(EXAMPLE / "horn-impedance.csv"), label, preserved["throat_radius"],
                                       [driver], target, top_n=1, length=preserved["length"]), preserved)
    feasible = [r for r in rows if r["model_feasible"]]
    rejected = [r for r in rows if not r["model_feasible"]]
    generate_auto_report(feasible, {label: str(EXAMPLE / "horn-impedance.csv")}, {driver.driver_id: driver},
                         None, target, str(out), top_n=5, rejected=rejected, rejection_stage="FEM",
                         step_files={label: str(EXAMPLE / "horn.step")})
    return rows, driver


def _live_report(run, out):
    """Re-rank a preserved auto run's refinement responses and regenerate its report."""
    auto = RESULTS / run / "outputs/auto"
    spec = json.loads((auto.parent / "resolved_specification.json").read_text())["parameters"]
    target = TargetSpec(spec["target_f_low"], spec["target_f_high"], voltage_rms=spec["voltage_rms"],
                        observation_distance_m=spec["observation_distance"], max_ripple_db=spec["max_ripple_db"],
                        max_compression_ratio=spec["max_compression_ratio"])
    preserved = json.loads((auto / "refinement/ranked_results.json").read_text())
    geometry = {}
    with open(auto / "candidates.csv") as handle:
        for row in csv.DictReader(handle):
            geometry[row["candidate_id"]] = {"profile": row["profile"], **{k: float(row[k]) for k in ("throat_radius", "mouth_radius", "length")}}
    for row in preserved["results"] + preserved["rejected"]:
        if row.get("profile") and row.get("length"):
            geometry.setdefault(row["horn_label"], {k: row[k] for k in ("profile", "throat_radius", "mouth_radius", "length")})
    drivers = {}
    for driver_id in dict.fromkeys(r["driver_id"] for r in preserved["results"]):
        try:
            drivers[driver_id] = load_driver(str(ROOT / spec["drivers_db"]), driver_id)
        except (KeyError, ValueError, FileNotFoundError):
            continue  # removed from the cleaned database since the run
    assert drivers
    solver_csvs, steps, rows = {}, {}, []
    for path in sorted((auto / "refinement").glob("*_results.csv")):
        label = path.name[:-len("_results.csv")]
        if label not in geometry:
            continue
        solver_csvs[label] = str(path)
        step = auto / "refinement" / f"{label}.step" if label.startswith("refine_") else auto / "geometry" / f"horn_{label}.step"
        steps[label] = str(step)
        g = geometry[label]
        rows += _annotate(rank_horn_drivers(str(path), label, g["throat_radius"], list(drivers.values()), target,
                                            top_n=len(drivers), length=g["length"]), g)
    feasible = [r for r in rows if r["model_feasible"]]
    rejected = [r for r in rows if not r["model_feasible"]]
    design = json.loads((auto / "design.json").read_text())
    generate_auto_report(feasible, solver_csvs, drivers, None, target, str(out), top_n=5, derived_geometry=design,
                         rejected=rejected, rejection_stage="FEM", step_files=steps)
    return feasible, rejected


def _capture_legends(monkeypatch):
    """Record every legend label of every figure embedded in the HTML report."""
    from horn_analysis import plot_theme
    legends = []
    original = plot_theme.fig_to_b64

    def capture(fig):
        legends.append([t.get_text() for ax in fig.axes if ax.get_legend() for t in ax.get_legend().get_texts()])
        return original(fig)
    monkeypatch.setattr(plot_theme, "fig_to_b64", capture)
    return legends


def _ranking_rows(html_text):
    return {int(rank): ident for rank, ident in re.findall(r'<tr data-rank="(\d+)" data-candidate="([^"]*)">', html_text)}


def _assert_legends_match_rankings(legends, html_text):
    rows = _ranking_rows(html_text)
    assert rows
    labels = [label for figure in legends for label in figure if label != "Target band"]
    assert labels, "no candidate legends found"
    for label in labels:
        match = re.fullmatch(r"#(\d+) (\S+)(?: \+ (\S+))?", label)
        assert match, f"legend {label!r} is not '#rank ID'"
        rank, horn, driver = int(match.group(1)), match.group(2), match.group(3)
        assert rank in rows, f"legend {label!r} names rank {rank}, which is not in the ranking table"
        if driver is None:
            assert rows[rank].split(" + ")[0] == horn, (label, rows[rank])
        else:
            assert rows[rank] == f"{horn} + {driver}", (label, rows[rank])
    # The electrical figure carries one legend entry per ranked candidate.
    candidate_labels = {label for label in labels if " + " in label}
    assert candidate_labels == {f"#{rank} {ident}" for rank, ident in rows.items()}


def _assert_profile_csv(out, row, rank=1):
    table = pd.read_csv(out / f"profile_{rank:02d}_{row['horn_label']}.csv")
    assert list(table.columns) == ["z_mm", "r_mm"]
    assert table.z_mm.iloc[0] == 0 and table.z_mm.iloc[-1] == pytest.approx(row["length"] * 1e3, rel=1e-12)
    expected = radius_profile(table.z_mm.to_numpy() / 1e3, row["throat_radius"], row["mouth_radius"], row["length"], row["profile"]) * 1e3
    np.testing.assert_allclose(table.r_mm.to_numpy(), expected, rtol=1e-12, atol=0)
    assert table.r_mm.iloc[0] == pytest.approx(row["throat_radius"] * 1e3, rel=1e-9)
    assert table.r_mm.iloc[-1] == pytest.approx(row["mouth_radius"] * 1e3, rel=1e-6)


def test_example_report_labels_profile_and_builder_data(tmp_path, monkeypatch):
    legends = _capture_legends(monkeypatch)
    rows, driver = _example_report(tmp_path)
    text = (tmp_path / "auto_report.html").read_text()
    _assert_legends_match_rankings(legends, text)
    _assert_profile_csv(tmp_path, rows[0])
    summary = (tmp_path / "auto_summary.txt").read_text()
    step = str(EXAMPLE / "horn.step")
    assert step in text and f"STEP (relative to this folder): {step}" in summary
    assert "File paths are relative to the folder that holds this report" in text
    assert "exceeded only because of an estimated front chamber is still ranked" in text
    assert f"Throat radius: {rows[0]['throat_radius']*1e3:.2f} mm" in summary
    # Experimental label, reading guide, formula, bound sources, chamber and evidence status.
    assert "Experimental prediction." in text and "EXPERIMENTAL PREDICTION" in summary
    assert "<h2>How to read this report</h2>" in text
    assert "clip((L_mean-60)/60,0,1)" in text and "Geometry bounds:" in text
    assert "front_chamber_estimated" in text
    assert "estimated 1" in text  # chamber volume shown with its source
    assert "feasible, evidence gaps" in text
    assert "|Ze|" in text and "Electrical impedance, excursion and group delay" in text


def test_report_lists_ripple_rejections_with_counts_and_margins(tmp_path, monkeypatch):
    """--max_ripple_db below the best result: every FEM candidate is rejected and the report says why."""
    target = TargetSpec(800, 1600, max_ripple_db=1.0)
    rows, _ = _example_report(tmp_path, target)
    assert rows and not any(r["model_feasible"] for r in rows)
    text = (tmp_path / "auto_report.html").read_text()
    assert "No feasible design in the evaluated set" in text
    assert "<h2>Rejected designs</h2>" in text
    assert 'data-reason="target_band_ripple_exceeded"' in text
    excess = rows[0]["passband_ripple_db"] - 1.0
    row = re.search(r'<tr data-reason="target_band_ripple_exceeded">(.*?)</tr>', text).group(1)
    assert f"<td>1</td><td>{excess:.2f}</td><td>{excess:.2f}</td><td>dB</td>" in row
    saved = json.loads((tmp_path / "auto_rejections.json").read_text())
    assert saved["stage"] == "FEM" and saved["rejected_pairs"] == 1
    assert saved["reasons"][0]["reason"] == "target_band_ripple_exceeded"
    assert saved["reasons"][0]["closest_excess"] == pytest.approx(excess, abs=1e-12)
    assert "target_band_ripple_exceeded: 1" in (tmp_path / "auto_summary.txt").read_text()


def test_empty_path_report_shows_stage_reason_and_prescreen_rejections(tmp_path):
    """report_no_drivers: LEM rows rejected before FEM, and the stop reason."""
    target = TargetSpec(500, 4000, max_compression_ratio=10.0)
    lem_rows = [
        {"driver_id": "a", "horn_label": "auto_con_0000", "rejection_reasons": ["compression_ratio_exceeded"],
         "driver_interface": {"compression_ratio": 12.5}, "model_feasible": False},
        {"driver_id": "b", "horn_label": "auto_con_0001", "rejection_reasons": ["compression_ratio_exceeded", "target_band_ripple_exceeded"],
         "driver_interface": {"compression_ratio": 30.0}, "passband_ripple_db": 9.0, "model_feasible": False},
    ]
    reason = "No searched geometry lies within the selected radiation model domain."
    generate_auto_report([], {}, {}, None, target, str(tmp_path), no_feasible_reason=reason,
                         rejected=lem_rows, rejection_stage="analytical prescreen")
    text = (tmp_path / "auto_report.html").read_text()
    assert reason in text and "at the analytical prescreen stage" in text
    row = re.search(r'<tr data-reason="compression_ratio_exceeded">(.*?)</tr>', text).group(1)
    assert "<td>2</td><td>2.50</td><td>20.00</td><td>Sd/St</td>" in row
    assert 'data-reason="target_band_ripple_exceeded"' in text
    # With nothing scored at all the stop reason is still stated.
    empty = tmp_path / "empty"
    generate_auto_report([], {}, {}, None, target, str(empty), no_feasible_reason="No drivers passed the initial screening for this request.")
    assert "No drivers passed the initial screening" in (empty / "auto_report.html").read_text()
    assert json.loads((empty / "auto_rejections.json").read_text())["rejected_pairs"] == 0


@pytest.mark.parametrize("run", LIVE_RUNS)
def test_live_run_report_regenerates_from_preserved_results(run, tmp_path, monkeypatch):
    if not (RESULTS / run / "outputs/auto/refinement/ranked_results.json").exists():
        pytest.skip(f"preserved run {run} not present")
    legends = _capture_legends(monkeypatch)
    feasible, rejected = _live_report(run, tmp_path)
    text = (tmp_path / "auto_report.html").read_text()
    if feasible:
        _assert_legends_match_rankings(legends, text)
        ranked = sorted(feasible, key=lambda r: r["composite_score"], reverse=True)
        _assert_profile_csv(tmp_path, ranked[0])
        step = re.search(r"STEP file \(acoustic air volume, mm; path relative to this report\)</dt><dd><code>([^<]+)</code>", text).group(1)
        assert Path(step).exists() and ranked[0]["horn_label"] in step
        assert len(list(tmp_path.glob("coupled_*.csv"))) == min(5, len(feasible))
    if rejected:
        assert "<h2>Rejected designs</h2>" in text
        reasons = {reason for r in rejected for reason in r["rejection_reasons"]}
        for reason in reasons:
            assert f'data-reason="{reason}"' in text


@pytest.mark.parametrize("imported", [False, True])
def test_single_mode_matches_auto_mode_levels_and_kpis(tmp_path, imported):
    """Same STEP response and driver: single-mode levels equal auto-mode levels to 1e-9 dB."""
    from horn_analysis.couple_single import couple
    preserved, driver = _example_inputs()
    target = TargetSpec(800, 1600)
    csv_path = str(EXAMPLE / "horn-impedance.csv")
    couple(csv_path, str(ROOT / "data/drivers-curated"), driver.driver_id,
           throat_radius=None if imported else preserved["throat_radius"], voltage=target.voltage_rms,
           output_csv=str(tmp_path / "coupled_spl.csv"), output_png=str(tmp_path / "coupled_spl.png"),
           output_kpis=str(tmp_path / "kpis.json"), imported_geometry=imported, target=target, length=preserved["length"])
    single = pd.read_csv(tmp_path / "coupled_spl.csv")
    kpis = json.loads((tmp_path / "kpis.json").read_text())

    auto_dir = tmp_path / "auto"
    rows, _ = _example_report(auto_dir, target)
    auto = pd.read_csv(auto_dir / f"coupled_01_{driver.driver_id}_{preserved['horn_label']}.csv")
    np.testing.assert_allclose(single.coupled_spl.to_numpy(), auto.spl.to_numpy(), rtol=0, atol=1e-9)
    for column in ("ze_mag_ohm", "excursion_peak_mm", "group_delay_ms"):
        np.testing.assert_allclose(single[column].to_numpy(), auto[column].to_numpy(), rtol=0, atol=1e-12)
    row = rows[0]
    for key in ("passband_ripple_db", "avg_sensitivity_db", "f3_low_hz", "f3_high_hz", "extension_octaves", "composite_score"):
        assert kpis[key] == pytest.approx(row[key], abs=1e-9), key
    assert kpis["eligibility_status"] == row["eligibility_status"]
    assert kpis["front_chamber_estimated"] is True
    assert kpis["output_metric"] == row["output_metric"] == "modal_baffled_on_axis"


def _cli_couple(tmp_path, monkeypatch, *extra, csv_path=None):
    """Run the couple_single command line as main.nf does and return its KPIs and levels."""
    from horn_analysis import couple_single
    preserved, driver = _example_inputs()
    argv = ["couple_single", "--solver-csv", str(csv_path or EXAMPLE / "horn-impedance.csv"),
            "--drivers-db", str(ROOT / "data/drivers-curated"), "--driver-id", driver.driver_id,
            "--throat-radius", str(preserved["throat_radius"]), "--length", str(preserved["length"]),
            "--output-csv", str(tmp_path / "c.csv"), "--output-png", str(tmp_path / "c.png"),
            "--output-kpis", str(tmp_path / "k.json"), *extra]
    monkeypatch.setattr("sys.argv", argv)
    couple_single.main()
    return json.loads((tmp_path / "k.json").read_text()), pd.read_csv(tmp_path / "c.csv")


@pytest.mark.parametrize("weights,min_output", [("ripple=1,extension=1", None), (None, 130.0)])
def test_single_mode_cli_uses_auto_mode_weights_and_output_floor(tmp_path, monkeypatch, weights, min_output):
    """--score-weights and --min-output reach single mode, so its score and eligibility match auto mode."""
    extra = ["--target-f-low", "800", "--target-f-high", "1600"]
    extra += ["--score-weights", weights] if weights else []
    extra += ["--min-output", str(min_output)] if min_output is not None else []
    kpis, _ = _cli_couple(tmp_path, monkeypatch, *extra)
    target = TargetSpec(800, 1600, score_weights=weights, min_output_db=min_output)
    rows, _ = _example_report(tmp_path / "auto", target)
    row = rows[0]
    assert kpis["eligibility_status"] == row["eligibility_status"]
    assert kpis["rejection_reasons"] == row["rejection_reasons"]
    assert kpis["composite_score"] == pytest.approx(row["composite_score"], abs=1e-9)
    if weights:
        default_rows, _ = _example_report(tmp_path / "default", TargetSpec(800, 1600))
        assert abs(row["composite_score"] - default_rows[0]["composite_score"]) > 1e-3
    else:
        assert "minimum_output_not_met" in kpis["rejection_reasons"]


def test_single_mode_cli_applies_distance_without_target_band(tmp_path, monkeypatch):
    """With no target band the sweep is the band, and --distance still moves the observer."""
    (tmp_path / "near").mkdir()
    (tmp_path / "far").mkdir()
    _, near_levels = _cli_couple(tmp_path / "near", monkeypatch, "--distance", "1")
    far, far_levels = _cli_couple(tmp_path / "far", monkeypatch, "--distance", "2")
    f = pd.read_csv(EXAMPLE / "horn-impedance.csv")["frequency"]
    assert far["target_band_hz"] == [f.min(), f.max()]
    assert far["target_band_clipped_to_sweep"] is False
    drop = (near_levels.coupled_spl - far_levels.coupled_spl).to_numpy()
    assert np.all(drop > 3.0)


def test_single_mode_sweep_outside_target_band_assesses_the_sweep(tmp_path):
    """A sweep that misses the target band falls back to the whole sweep instead of raising."""
    from horn_analysis.couple_single import single_assessment
    preserved, driver = _example_inputs()
    df = pd.read_csv(EXAMPLE / "horn-impedance.csv")
    df = df[df.frequency >= 1700].reset_index(drop=True)
    _, _, assessment = single_assessment(df, driver, TargetSpec(800, 1600), preserved["throat_radius"], preserved["length"])
    assert assessment["target_band_outside_sweep"] is True
    assert assessment["target_band_clipped_to_sweep"] is False
    assert assessment["target_band_requested_hz"] == [800, 1600]
    assert assessment["target_band_hz"] == [df.frequency.min(), df.frequency.max()]
    assert np.isfinite(assessment["passband_ripple_db"])


def _single_report(tmp_path, kpis):
    from horn_analysis.single_report import generate_single_report
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    png = tmp_path / "blank.png"
    plt.figure(figsize=(1, 1)); plt.savefig(png); plt.close()
    preserved, _ = _example_inputs()
    return generate_single_report(preserved["throat_radius"], preserved["mouth_radius"], preserved["length"],
                                  preserved.get("profile", "horn"), {}, str(png), str(png), str(png), str(png), str(png),
                                  coupled_png=str(png), coupled_kpis=kpis)


@pytest.mark.parametrize("band,note", [((800, 4000), "(clipped to the sweep)"),
                                       ((10000, 16000), "(outside the sweep; whole sweep assessed)")])
def test_single_report_states_experimental_label_chamber_and_band(tmp_path, band, note):
    """The single report carries the experimental label, evidence status, chamber and band handling."""
    from horn_analysis.couple_single import couple
    preserved, driver = _example_inputs()
    couple(str(EXAMPLE / "horn-impedance.csv"), str(ROOT / "data/drivers-curated"), driver.driver_id,
           throat_radius=preserved["throat_radius"], output_csv=str(tmp_path / "c.csv"),
           output_png=str(tmp_path / "c.png"), output_kpis=str(tmp_path / "k.json"),
           target=TargetSpec(*band), length=preserved["length"])
    kpis = json.loads((tmp_path / "k.json").read_text())
    assert kpis["front_chamber_estimated"] is True
    text = _single_report(tmp_path, kpis)
    assert "<strong>Experimental prediction.</strong>" in text
    assert note in text
    assert "Evidence status: " + kpis["eligibility_status"].replace("_", " ") in text
    volume = kpis["driver_interface"]["front_chamber_volume_cm3"]
    assert f"Front chamber: estimated {volume:.1f} cm3" in text
