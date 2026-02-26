"""Tests for the integrated ranking pipeline."""

import pytest
import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_analysis.scoring import TargetSpec
from horn_analysis.rank_pipeline import rank_horn_drivers


def _make_driver(driver_id, bl_tm=8.0, sd_m2=0.0008, fs_hz=500.0, **kwargs):
    """Helper to create a DriverParameters with sensible defaults."""
    return DriverParameters(
        driver_id=driver_id,
        manufacturer="Test",
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


@pytest.fixture
def solver_csv(tmp_path):
    """Create a synthetic solver CSV with impedance data."""
    csv_path = tmp_path / "solver_results.csv"
    freq = np.geomspace(500, 8000, 50)
    # Simulated bandpass response peaking around 2kHz
    spl = 90 - 10 * ((np.log10(freq / 2000)) ** 2) * 20
    spl = np.clip(spl, 50, 95)
    # Simulated impedance
    z_real = 400 + 50 * np.sin(2 * np.pi * np.log10(freq))
    z_imag = 100 * np.cos(2 * np.pi * np.log10(freq))
    pd.DataFrame({
        "frequency": freq,
        "spl": spl,
        "z_real": z_real,
        "z_imag": z_imag,
    }).to_csv(csv_path, index=False)
    return str(csv_path)


@pytest.fixture
def target():
    return TargetSpec(f_low_hz=500, f_high_hz=4000)


class TestRankHornDrivers:
    def test_returns_ranked_results(self, solver_csv, target):
        """Should return a list of ranked dicts."""
        drivers = [
            _make_driver("d1", bl_tm=8.0),
            _make_driver("d2", bl_tm=10.0),
            _make_driver("d3", bl_tm=6.0),
        ]
        results = rank_horn_drivers(
            solver_csv=solver_csv,
            horn_label="conical",
            throat_radius=0.025,
            drivers=drivers,
            target=target,
            top_n=3,
        )
        assert len(results) == 3
        assert all("composite_score" in r for r in results)
        assert all("driver_id" in r for r in results)
        assert all("kpi" in r for r in results)
        assert all("horn_label" in r for r in results)

    def test_results_sorted_by_score(self, solver_csv, target):
        """Results should be sorted by composite_score descending."""
        drivers = [
            _make_driver("d1", bl_tm=8.0),
            _make_driver("d2", bl_tm=12.0),
            _make_driver("d3", bl_tm=4.0),
        ]
        results = rank_horn_drivers(
            solver_csv=solver_csv,
            horn_label="exponential",
            throat_radius=0.025,
            drivers=drivers,
            target=target,
            top_n=3,
        )
        scores = [r["composite_score"] for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_higher_bl_ranks_higher(self, solver_csv, target):
        """Higher BL (force factor) generally produces higher sensitivity -> higher rank."""
        low_bl = _make_driver("low_bl", bl_tm=3.0)
        high_bl = _make_driver("high_bl", bl_tm=15.0)
        results = rank_horn_drivers(
            solver_csv=solver_csv,
            horn_label="conical",
            throat_radius=0.025,
            drivers=[low_bl, high_bl],
            target=target,
            top_n=2,
        )
        # Higher BL should produce higher sensitivity and rank first
        assert results[0]["driver_id"] == "high_bl"

    def test_top_n_limits_results(self, solver_csv, target):
        """top_n should limit the number of results returned."""
        drivers = [_make_driver(f"d{i}", bl_tm=float(5 + i)) for i in range(10)]
        results = rank_horn_drivers(
            solver_csv=solver_csv,
            horn_label="conical",
            throat_radius=0.025,
            drivers=drivers,
            target=target,
            top_n=3,
        )
        assert len(results) == 3

    def test_result_contains_metadata(self, solver_csv, target):
        """Each result should contain manufacturer and model_name."""
        d = _make_driver("test_driver")
        results = rank_horn_drivers(
            solver_csv=solver_csv,
            horn_label="hyperbolic",
            throat_radius=0.025,
            drivers=[d],
            target=target,
            top_n=1,
        )
        assert results[0]["manufacturer"] == "Test"
        assert results[0]["model_name"] == "test_driver"

    def test_kpi_values_reasonable(self, solver_csv, target):
        """KPIs in results should have reasonable values."""
        d = _make_driver("test_kpi", bl_tm=10.0)
        results = rank_horn_drivers(
            solver_csv=solver_csv,
            horn_label="conical",
            throat_radius=0.025,
            drivers=[d],
            target=target,
            top_n=1,
        )
        kpi = results[0]["kpi"]
        assert kpi["peak_spl_db"] > 0
        assert kpi["peak_frequency_hz"] > 0
