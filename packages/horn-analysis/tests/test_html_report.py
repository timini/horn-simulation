"""Tests for the HTML report generator."""

import numpy as np
import pandas as pd
import pytest

from horn_core.parameters import DriverParameters
from horn_analysis.scoring import TargetSpec
from horn_analysis.html_report import generate_html_report


def _make_driver(driver_id, bl_tm=8.0, sd_m2=0.0008, fs_hz=500.0, **kwargs):
    return DriverParameters(
        driver_id=driver_id,
        manufacturer="Test Mfg",
        model_name=driver_id,
        fs_hz=fs_hz,
        re_ohm=kwargs.get("re_ohm", 6.0),
        bl_tm=bl_tm,
        sd_m2=sd_m2,
        mms_kg=kwargs.get("mms_kg", 0.003),
        le_h=kwargs.get("le_h", 0.0005),
        qms=kwargs.get("qms", 5.0),
        qes=kwargs.get("qes", 0.45),
    )


def _make_solver_csv(path, profile="conical"):
    """Create a synthetic solver CSV with SPL, impedance, and phase columns."""
    freq = np.geomspace(500, 8000, 60)
    spl = 90 - 10 * ((np.log10(freq / 2000)) ** 2) * 20
    spl = np.clip(spl, 50, 95)
    z_real = 400 + 50 * np.sin(2 * np.pi * np.log10(freq))
    z_imag = 100 * np.cos(2 * np.pi * np.log10(freq))
    phase_deg = np.degrees(np.arctan2(z_imag, z_real))
    pd.DataFrame({
        "frequency": freq,
        "spl": spl,
        "z_real": z_real,
        "z_imag": z_imag,
        "phase_deg": phase_deg,
    }).to_csv(path, index=False)
    return str(path)


def _make_coupled_csv(path):
    """Create a synthetic coupled-SPL CSV."""
    freq = np.geomspace(500, 8000, 60)
    spl = 85 + 5 * np.sin(np.log10(freq) * 4)
    pd.DataFrame({"frequency": freq, "spl": spl}).to_csv(path, index=False)
    return str(path)


@pytest.fixture
def report_inputs(tmp_path):
    """Build all inputs needed for generate_html_report."""
    solver_csvs = {}
    for profile in ("conical", "exponential", "hyperbolic"):
        csv_path = tmp_path / f"{profile}_results.csv"
        solver_csvs[profile] = _make_solver_csv(csv_path, profile)

    drivers = {
        "d1": _make_driver("d1", bl_tm=8.0),
        "d2": _make_driver("d2", bl_tm=10.0),
    }

    target = TargetSpec(f_low_hz=500, f_high_hz=4000)
    throat_radius = 0.025

    all_ranked = []
    for i, (did, drv) in enumerate(drivers.items()):
        for profile in solver_csvs:
            all_ranked.append({
                "driver_id": did,
                "horn_label": profile,
                "manufacturer": drv.manufacturer,
                "model_name": drv.model_name,
                "composite_score": 0.85 - i * 0.05 - list(solver_csvs).index(profile) * 0.02,
                "extension_octaves": 0.40 - i * 0.1,
                "f3_low_hz": 430.0,
                "f3_high_hz": 4600.0,
                "f3_low_is_bound": False,
                "f3_high_is_bound": i == 0,
                "passband_ripple_db": 2.5 + i * 0.5,
                "avg_sensitivity_db": 92.0 - i * 2.0,
                "kpi": {
                    "f3_low_hz": 520.0,
                    "f3_high_hz": 3800.0,
                    "peak_spl_db": 95.0,
                    "peak_frequency_hz": 2000.0,
                },
            })
    all_ranked.sort(key=lambda r: r["composite_score"], reverse=True)

    csv_pairs = []
    for rank, r in enumerate(all_ranked[:3], 1):
        csv_path = tmp_path / f"coupled_{rank:02d}.csv"
        label = f"#{rank} {r['manufacturer']} {r['model_name']} ({r['horn_label']})"
        csv_pairs.append((_make_coupled_csv(csv_path), label))

    return {
        "all_ranked": all_ranked,
        "solver_csvs": solver_csvs,
        "drivers": drivers,
        "throat_radius": throat_radius,
        "target": target,
        "csv_pairs": csv_pairs,
        "top_n": 5,
    }


class TestHtmlReport:
    def test_returns_string(self, report_inputs):
        """Should return a non-trivial HTML string."""
        result = generate_html_report(**report_inputs)
        assert isinstance(result, str)
        assert len(result) > 1000

    def test_html_contains_key_sections(self, report_inputs):
        """HTML should contain DOCTYPE, title, rankings, drivers, scores."""
        result = generate_html_report(**report_inputs)
        assert "<!DOCTYPE html>" in result
        assert "Horn Auto-Select Report" in result
        assert "Rankings" in result
        assert "Pre-Screened Drivers" in result
        # Check driver data appears
        assert "Test Mfg" in result
        # Check score values appear
        assert "0.85" in result or "0.850" in result

    def test_html_states_the_stress_power(self, report_inputs):
        """With test_power_w the page states the level and each candidate's stress results."""
        assert "Stress check" not in generate_html_report(**report_inputs)
        report_inputs["target"] = TargetSpec(f_low_hz=500, f_high_hz=4000, test_power_w=50.)
        report_inputs["all_ranked"][0].update({
            "test_power_w": 50., "test_voltage_rms": 20., "test_power_band_max_excursion_mm": 0.2,
            "test_power_xmax_margin_mm": 2.8, "test_power_band_max_input_power_w": 72.4,
            "test_power_rating_margin_w": 127.6})
        result = generate_html_report(**report_inputs)
        assert "Stress check: 50 W into each driver's nominal impedance" in result
        assert "in-band excursion 0.20 mm, Xmax margin 2.80 mm; in-band input power up to 72.4 W, rating margin 127.6 W" in result

    def test_html_has_four_embedded_images(self, report_inputs):
        """Report should contain exactly 4 base64-embedded PNG images."""
        result = generate_html_report(**report_inputs)
        count = result.count("data:image/png;base64,")
        assert count == 4, f"Expected 4 embedded images, found {count}"

    def test_html_is_writable(self, report_inputs, tmp_path):
        """Generated HTML should be writable and produce a substantial file."""
        result = generate_html_report(**report_inputs)
        out_file = tmp_path / "report.html"
        out_file.write_text(result)
        assert out_file.stat().st_size > 50_000


def test_variable_throats_appear_in_table_summary_and_candidate_renders(report_inputs, monkeypatch):
    import horn_analysis.html_report as report
    radii = [.012, .018, .024]
    for i, row in enumerate(report_inputs['all_ranked']):
        row.update(throat_radius=radii[i % 3], mouth_radius=.1, length=.2, profile=row['horn_label'])
    rendered = []
    monkeypatch.setattr(report, 'fig_to_b64_3d', lambda **kw: rendered.append(kw['throat_radius']) or 'data:image/png;base64,')
    text = generate_html_report(**report_inputs, mouth_radius=.1, length=.2,
                               derived_geometry={'throat_radius_range':[.012,.024],
                                                 'mouth_radius_range':[.1,.1], 'length_range':[.2,.2]})
    assert 'Throat radius: 0.0120 to 0.0240 m (search range)' in text
    # Builder dimensions in the ranking table are in millimetres.
    assert '<th>Throat R (mm)</th>' in text
    assert '<td>12.0</td>' in text and '<td>18.0</td>' in text
    assert '<dt>Throat radius range</dt><dd>0.0120 to 0.0240 m</dd>' in text
    assert set(rendered) == set(radii)
    assert 'Mouth: 0.1000 m' in text
    assert 'Length: 0.200 m' in text


@pytest.mark.parametrize('mouth_bounds,length_bounds,mouth_text,length_text', [
    ([.1, .1], [.2, .3], '0.1000 m', '0.200 to 0.300 m (search range)'),
    ([.1, .15], [.2, .2], '0.1000 to 0.1500 m (search range)', '0.200 m'),
])
def test_fixed_and_variable_dimensions_display_independently(report_inputs, mouth_bounds, length_bounds, mouth_text, length_text):
    text = generate_html_report(**report_inputs, derived_geometry={
        'mouth_radius_range': mouth_bounds, 'length_range': length_bounds})
    assert f'Mouth: {mouth_text}' in text
    assert f'Length: {length_text}' in text


def test_rankings_use_profile_metadata_with_legacy_label_fallback(report_inputs):
    from horn_analysis.html_report import _render_rankings_rows
    row = dict(report_inputs['all_ranked'][0], horn_label='auto_hyp_0012', profile='hyperbolic')
    text = _render_rankings_rows([row], report_inputs['drivers'])
    assert '>Hyperbolic</span>' in text
    # The candidate ID has its own column; the profile badge names the profile, not the ID.
    assert '<td><code>auto_hyp_0012 + d1</code></td>' in text
    assert '>Auto_hyp_0012</span>' not in text
    row.update(profile='', horn_label='conical')
    assert '>Conical</span>' in _render_rankings_rows([row], report_inputs['drivers'])


def test_design_summary_names_the_scaled_tractrix_profile():
    """The old "lecleach" profile is an axially scaled tractrix, not a Le Cleac'h horn."""
    from horn_analysis.html_report import _render_design_summary
    text = _render_design_summary({})
    assert "axially scaled tractrix" in text
    assert "Cl&eacute;ac" not in text


def test_imported_report_uses_cad_areas_without_parametric_dimensions(tmp_path):
    from horn_analysis.single_report import generate_single_report
    import base64
    png = tmp_path/'image.png'
    png.write_bytes(base64.b64decode('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+aGUcAAAAASUVORK5CYII='))
    csv = tmp_path/'results.csv'
    pd.DataFrame({'frequency':[100,200], 'inlet_area_m2':[.001,.001], 'mouth_area_m2':[.008,.008]}).to_csv(csv,index=False)
    text = generate_single_report(.05,.2,.08,'conical',{},*[str(png)]*5,
                                  final_csv=str(csv),imported_geometry=True)
    assert 'Inlet area (CAD)</dt><dd>0.001 m²' in text
    assert '<dt>Throat radius</dt>' not in text
    assert '<dt>Mouth diameter</dt>' not in text
    assert 'Throat: 100.0 mm' not in text
    assert 'Imported step' in text


def test_report_states_score_formula_weights_and_band_referenced_f3(report_inputs):
    from horn_analysis.scoring import describe_scoring
    target = TargetSpec(500, 4000, min_output_db=90., score_weights="size=0.2")
    scoring = describe_scoring(target, {"length_range": [.1, .2]}, {"max_length": .2})
    text = generate_html_report(**{**report_inputs, "target": target}, scoring=scoring)
    assert "Extension (oct; ≥ = lower bound)" in text and "BW Cov." not in text
    assert "&plusmn;3 dB flat band" in text
    assert "+0.40" in text
    # Rankings show f3 from the target-band reference, not the peak-lobe KPI.
    assert "<td>430 to \u22654600</td>" in text and "<td>520 to" not in text
    assert "clip((L_mean-60)/60,0,1)" in text
    assert "size 0.167" in text
    assert "heuristic range with a user limit" in text
    assert "min_output_db 90.0 (user)" in text


def test_report_embeds_the_electrical_plot_when_one_is_written(report_inputs, tmp_path):
    import base64
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([1, 2], [3, 4])
    png = tmp_path / "auto_electrical.png"
    fig.savefig(png)
    plt.close(fig)
    text = generate_html_report(**report_inputs, electrical_png=str(png))
    encoded = base64.b64encode(png.read_bytes()).decode("ascii")
    assert f'<img src="data:image/png;base64,{encoded}"' in text
    assert "Electrical impedance, excursion and group delay" in text
    # Without an electrical plot the section is left out, not shown empty.
    plain = generate_html_report(**report_inputs)
    assert "Electrical impedance, excursion and group delay" not in plain


def test_single_report_states_flanged_piston_error_only_when_that_model_was_used(tmp_path):
    from horn_analysis.single_report import generate_single_report
    import base64
    png = tmp_path/'image.png'
    png.write_bytes(base64.b64decode('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+aGUcAAAAASUVORK5CYII='))
    texts = {}
    for model in ('flanged_piston', 'modal_baffled'):
        csv = tmp_path/f'{model}.csv'
        pd.DataFrame({'frequency':[100,200], 'inlet_area_m2':[.001,.001], 'mouth_area_m2':[.008,.008],
                      'radiation_model':[model, model]}).to_csv(csv,index=False)
        texts[model] = generate_single_report(.05,.2,.08,'conical',{},*[str(png)]*5,
                                              final_csv=str(csv),imported_geometry=True)
    assert 'Known radiation-model error' in texts['flanged_piston']
    assert '1.05 dB' in texts['flanged_piston']
    assert 'Known radiation-model error' not in texts['modal_baffled']


def test_scoring_section_states_size_compromise_refinement_and_chamber_policy():
    from horn_analysis.html_report import _render_scoring, _chamber_note
    from horn_analysis.scoring import TargetSpec, describe_scoring
    from horn_core.geometry_designer import generate_auto_candidates
    _, design = generate_auto_candidates(800, 1600, [.02], max_length=.1)
    scoring = describe_scoring(TargetSpec(800, 1600), {"size_compromise": design.size_compromise})
    scoring["refinement"] = {"seed_candidate": "auto_os_0001", "seed_source": "least_infeasible", "new_evaluations": 4,
                             "budget": 6, "bound_sides_reached": ["length_high"],
                             "widened_bounds": [{"field": "length", "side": "high", "from": [.1, .2], "to": [.1, .25]}]}
    page = _render_scoring(scoring)
    assert "858 Hz" in page and "0.10 octaves" in page
    assert "least infeasible" in page and "length high to 0.2500" in page
    assert "front_chamber_estimated" in page
    note = _chamber_note({"front_chamber_estimated": True, "without_front_chamber": {"passband_ripple_db": 4.29}})
    assert "est. chamber" in note and "4.3 without" in note
    assert _chamber_note({"front_chamber_estimated": False}) == ""


def test_size_compromise_is_worded_as_a_sizing_rule_and_checked_against_predictions():
    from horn_analysis.html_report import _render_scoring
    from horn_analysis.scoring import TargetSpec, describe_scoring, size_compromise_note
    from horn_core.geometry_designer import generate_auto_candidates
    _, design = generate_auto_candidates(800, 1600, [.02], max_length=.1)
    compromise = design.size_compromise
    text = size_compromise_note(compromise)
    assert "858 Hz" in text and "rule" in text and "not a prediction" in text
    rows = [{"model_feasible": True, "composite_score": .5, "f3_low_hz": 590.},
            {"model_feasible": True, "composite_score": .2, "f3_low_hz": 900.}]
    assert "does not show this shortfall" in size_compromise_note(compromise, rows)
    assert "does not show" not in size_compromise_note(compromise, rows[1:])
    page = _render_scoring(describe_scoring(TargetSpec(800, 1600), {"size_compromise": compromise}), rows)
    assert "590 Hz" in page and "does not show this shortfall" in page
    assert "support about" not in page


def test_evidence_summary_states_each_driver_once_and_header_names_horns(report_inputs):
    """Three ranked horns with one driver used to repeat the same evidence sentence three times."""
    gaps = ["driver_interface_unverified", "front_chamber_volume_estimated"]
    for row in report_inputs["all_ranked"]:
        row["evidence_gaps"] = list(gaps)
    text = generate_html_report(**report_inputs)
    summary = text.split("<ul>", 1)[1].split("</ul>", 1)[0]
    # Two drivers, three horns each: one entry per driver with the same gaps.
    assert summary.count("<li>") == 2
    # The header lists horn candidate IDs, so it must not call them profiles.
    assert "Horns: conical, exponential, hyperbolic" in text
    assert "Profiles: conical" not in text
