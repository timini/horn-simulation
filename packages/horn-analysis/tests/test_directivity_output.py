"""Experimental infinite-baffle directivity in the auto report."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import j1

from horn_analysis.directivity import directivity_table, DIRECTIVITY_LABEL
from horn_analysis.directivity_plot import compute_beamwidth
from horn_analysis.scoring import TargetSpec

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "examples/800-1600-modal"


def uniform_frame():
    frame = pd.DataFrame(dict(frequency=[1000., 4000.], radiation_model="modal_baffled", modal_mode_count=16,
                              mouth_area_m2=np.pi * 0.1**2, air_c_m_s=343., air_rho_kg_m3=1.225))
    for n in range(16):
        frame[f"modal_v_{n}_real"] = 0.02 if n == 0 else 0.
        frame[f"modal_v_{n}_imag"] = 0.
    return frame


def test_uniform_mode_table_is_the_piston_pattern():
    table = directivity_table(uniform_frame(), angles_deg=[30., 60., 90.])
    assert set(table.status) == {DIRECTIVITY_LABEL}
    for f, rows in table.groupby("frequency"):
        x = 2 * np.pi * f / 343. * 0.1 * np.sin(np.radians(rows.theta_deg.to_numpy()))
        expected = np.r_[0., 20 * np.log10(np.abs(2 * j1(x[1:]) / x[1:]))]
        np.testing.assert_allclose(rows.spl_db.to_numpy(), expected, atol=1e-6)


def test_local_radiation_models_have_no_directivity():
    frame = uniform_frame()
    frame["radiation_model"] = "flanged_piston"
    assert directivity_table(frame) is None


def test_preserved_modal_horn_narrows_with_frequency():
    table = directivity_table(pd.read_csv(EXAMPLE / "horn-impedance.csv"))
    assert np.isfinite(table.spl_db).all()
    assert (table[table.theta_deg == 0].spl_db == 0).all()
    f, width = compute_beamwidth(table)
    low, high = width[f < 700].mean(), width[f > 1900].mean()
    assert high < low


def ranked_example(tmp_path):
    from horn_drivers.loader import load_driver
    from horn_analysis.rank_pipeline import rank_horn_drivers
    preserved = json.loads((EXAMPLE / "ranking.json").read_text())[0]
    driver = load_driver(str(ROOT / "data/drivers-curated"), preserved["driver_id"])
    target = TargetSpec(800, 1600)
    label = preserved["horn_label"]
    ranked = rank_horn_drivers(str(EXAMPLE / "horn-impedance.csv"), label, preserved["throat_radius"],
                               [driver], target, top_n=1, length=preserved["length"])
    for key in ("profile", "throat_radius", "mouth_radius", "length"):
        ranked[0][key] = preserved[key]
    return ranked, driver, target, label


def test_auto_report_writes_labelled_directivity_for_modal_horns(tmp_path):
    from horn_analysis.auto_report import generate_auto_report
    ranked, driver, target, label = ranked_example(tmp_path)
    generate_auto_report(ranked, {label: str(EXAMPLE / "horn-impedance.csv")}, {driver.driver_id: driver},
                         ranked[0]["throat_radius"], target, str(tmp_path), top_n=1)
    table = pd.read_csv(tmp_path / f"directivity_{label}.csv")
    assert table.theta_deg.max() == 90. and set(table.status) == {DIRECTIVITY_LABEL}
    assert (tmp_path / "auto_beamwidth.png").stat().st_size > 0
    html = (tmp_path / "auto_report.html").read_text()
    assert "Directivity (experimental, infinite baffle)" in html
    assert "Known radiation-model error" not in html
    assert DIRECTIVITY_LABEL in (tmp_path / "auto_summary.txt").read_text()


def test_flanged_piston_report_states_the_known_impedance_error(tmp_path):
    from horn_analysis.auto_report import generate_auto_report
    ranked, driver, target, label = ranked_example(tmp_path)
    ranked[0]["radiation_model"] = "flanged_piston"
    generate_auto_report(ranked, {label: str(EXAMPLE / "horn-impedance.csv")}, {driver.driver_id: driver},
                         ranked[0]["throat_radius"], target, str(tmp_path), top_n=1)
    for name in ("auto_report.html", "auto_summary.txt"):
        text = (tmp_path / name).read_text()
        assert "1.05 dB" in text and "5.7 degrees" in text
