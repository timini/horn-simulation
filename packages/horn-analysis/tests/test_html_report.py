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
                "bandwidth_coverage": 0.90 - i * 0.1,
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
