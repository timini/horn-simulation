"""Tests for the driver–horn transfer function module."""

import numpy as np
import pandas as pd
import pytest

from horn_core.parameters import DriverParameters
from horn_analysis.transfer_function import (
    compute_driver_response,
    scale_solver_spl,
    screen_all_drivers,
)


def _make_driver(**overrides):
    """Helper to create a test driver with sensible defaults."""
    defaults = dict(
        driver_id="test",
        manufacturer="Test",
        model_name="T1",
        fs_hz=500.0,
        re_ohm=6.0,
        bl_tm=8.5,
        sd_m2=0.0008,
        mms_kg=0.003,
        le_h=0.0006,
        qms=5.0,
        qes=0.45,
    )
    defaults.update(overrides)
    return DriverParameters(**defaults)


class TestComputeDriverResponse:
    """Tests for the core transfer function computation."""

    def test_returns_correct_shape(self):
        drv = _make_driver()
        freq = np.geomspace(500, 8000, 50)
        z_real = np.full_like(freq, 420.0)  # ρc ≈ 420 Pa·s/m
        z_imag = np.zeros_like(freq)
        throat_area = np.pi * 0.025 ** 2

        p = compute_driver_response(drv, freq, z_real, z_imag, throat_area)
        assert p.shape == freq.shape
        assert np.iscomplexobj(p)

    def test_rigid_piston_in_infinite_tube(self):
        """With Z_horn = ρc (plane wave), verify non-trivial response."""
        drv = _make_driver()
        freq = np.geomspace(500, 8000, 100)
        rho_c = 343.0 * 1.225  # ≈ 420 Pa·s/m
        z_real = np.full_like(freq, rho_c)
        z_imag = np.zeros_like(freq)
        throat_area = np.pi * 0.025 ** 2

        p = compute_driver_response(drv, freq, z_real, z_imag, throat_area)
        spl = 20 * np.log10(np.abs(p) / 20e-6 + 1e-30)

        # SPL should be in a reasonable range for a compression driver
        assert np.all(spl > 50), f"SPL too low: min={spl.min():.1f}"
        assert np.all(spl < 160), f"SPL too high: max={spl.max():.1f}"

    def test_doubling_bl_increases_spl_by_6db(self):
        """Doubling BL should increase SPL by ~6 dB (linear force, quadratic power)."""
        freq = np.geomspace(500, 4000, 50)
        rho_c = 343.0 * 1.225
        z_real = np.full_like(freq, rho_c)
        z_imag = np.zeros_like(freq)
        throat_area = np.pi * 0.025 ** 2

        drv1 = _make_driver(bl_tm=4.0)
        drv2 = _make_driver(bl_tm=8.0)

        p1 = compute_driver_response(drv1, freq, z_real, z_imag, throat_area)
        p2 = compute_driver_response(drv2, freq, z_real, z_imag, throat_area)

        spl1 = 20 * np.log10(np.abs(p1) + 1e-30)
        spl2 = 20 * np.log10(np.abs(p2) + 1e-30)
        delta = np.mean(spl2 - spl1)

        # At mid frequencies where motional impedance doesn't dominate,
        # doubling BL gives ~6 dB. Allow wider tolerance since the
        # relationship is approximate (depends on Z_mot vs Z_e balance).
        assert 3.0 < delta < 12.0, f"BL doubling gave {delta:.1f} dB change (expected ~6)"

    def test_zero_le_purely_resistive(self):
        """With Le=0, the electrical impedance is purely resistive."""
        drv = _make_driver(le_h=0.0)
        freq = np.geomspace(500, 8000, 50)
        rho_c = 343.0 * 1.225
        z_real = np.full_like(freq, rho_c)
        z_imag = np.zeros_like(freq)
        throat_area = np.pi * 0.025 ** 2

        p = compute_driver_response(drv, freq, z_real, z_imag, throat_area)
        # Should still produce valid output
        assert p.shape == freq.shape
        assert np.all(np.isfinite(p))
        assert np.all(np.abs(p) > 0)


class TestScaleSolverSpl:
    def test_unit_pressure_gives_no_change(self):
        """p_throat = 1 Pa → 0 dB correction → SPL unchanged."""
        solver_spl = np.array([90.0, 92.0, 88.0])
        p_throat = np.array([1.0 + 0j, 1.0 + 0j, 1.0 + 0j])
        result = scale_solver_spl(solver_spl, p_throat)
        np.testing.assert_allclose(result, solver_spl, atol=1e-10)

    def test_10x_pressure_adds_20db(self):
        solver_spl = np.array([90.0])
        p_throat = np.array([10.0 + 0j])
        result = scale_solver_spl(solver_spl, p_throat)
        assert result[0] == pytest.approx(110.0, abs=1e-10)

    def test_handles_complex_pressure(self):
        solver_spl = np.array([90.0])
        # |3 + 4j| = 5
        p_throat = np.array([3.0 + 4.0j])
        result = scale_solver_spl(solver_spl, p_throat)
        expected = 90.0 + 20.0 * np.log10(5.0)
        assert result[0] == pytest.approx(expected, abs=1e-10)


class TestScreenAllDrivers:
    def test_batch_screening(self, tmp_path):
        """Test screening multiple drivers against a solver CSV."""
        # Create synthetic solver output
        freq = np.geomspace(500, 8000, 20)
        rho_c = 343.0 * 1.225
        csv_path = tmp_path / "solver.csv"
        pd.DataFrame({
            "frequency": freq,
            "spl": np.full_like(freq, 90.0),
            "z_real": np.full_like(freq, rho_c),
            "z_imag": np.zeros_like(freq),
        }).to_csv(csv_path, index=False)

        drivers = [
            _make_driver(driver_id="d1", bl_tm=8.0),
            _make_driver(driver_id="d2", bl_tm=12.0),
        ]

        results = screen_all_drivers(str(csv_path), drivers, throat_radius=0.025)
        assert len(results) == 2
        assert list(results.columns) == [
            "driver_id", "manufacturer", "model_name",
            "avg_spl_db", "peak_spl_db", "min_spl_db", "ripple_db",
        ]
        # Higher BL driver should have higher average SPL
        d1_avg = results.loc[results["driver_id"] == "d1", "avg_spl_db"].iloc[0]
        d2_avg = results.loc[results["driver_id"] == "d2", "avg_spl_db"].iloc[0]
        assert d2_avg > d1_avg
