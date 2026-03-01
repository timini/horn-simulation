"""Tests for driver pre-screening module."""

import math

import pytest
import numpy as np

from horn_core.parameters import DriverParameters
from horn_analysis.prescreen import PrescreenConfig, PrescreenResult, prescreen_drivers


def _make_driver(
    driver_id="test",
    fs_hz=500.0,
    qes=0.45,
    qms=5.0,
    sd_m2=0.0008,
    driver_type="compression",
    **kwargs,
):
    """Helper to create a DriverParameters with sensible defaults."""
    return DriverParameters(
        driver_id=driver_id,
        manufacturer="Test",
        model_name=driver_id,
        fs_hz=fs_hz,
        re_ohm=kwargs.get("re_ohm", 6.0),
        bl_tm=kwargs.get("bl_tm", 8.0),
        sd_m2=sd_m2,
        mms_kg=kwargs.get("mms_kg", 0.003),
        le_h=kwargs.get("le_h", 0.0005),
        qms=qms,
        qes=qes,
        driver_type=driver_type,
    )


@pytest.fixture
def config():
    return PrescreenConfig(
        target_f_low_hz=500,
        target_f_high_hz=4000,
        mouth_radius_m=0.2,
        length_m=0.5,
    )


class TestPrescreenDrivers:
    def test_basic_filtering(self, config):
        """Drivers with suitable fs, EBP, and type should pass."""
        good = _make_driver("good", fs_hz=400, qes=0.4, driver_type="compression")
        # EBP = 400/0.4 = 1000 > 50 -- passes
        result = prescreen_drivers([good], config)
        assert result.count == 1
        assert result.drivers[0].driver_id == "good"

    def test_fs_too_high_filtered(self, config):
        """Driver with fs above target_f_low * 1.5 should be filtered."""
        # target_f_low = 500, threshold = 750
        bad = _make_driver("bad_fs", fs_hz=800, qes=0.4)
        result = prescreen_drivers([bad], config)
        assert result.count == 0

    def test_ebp_too_low_filtered(self, config):
        """Driver with EBP below threshold should be filtered."""
        # EBP = 400 / 10.0 = 40 < 50
        bad = _make_driver("bad_ebp", fs_hz=400, qes=10.0)
        result = prescreen_drivers([bad], config)
        assert result.count == 0

    def test_large_cone_filtered_for_high_freq(self, config):
        """Large cone (e.g. 15") should be filtered when f_piston*factor < target_f_high."""
        # 15" driver: Sd ~ 855 cm² = 0.0855 m²
        # f_piston = 343/(2π·√(0.0855/π)) = ~331 Hz, ×10 = 3310 < 4000 → filtered
        big_cone = _make_driver("big_cone", fs_hz=400, qes=0.4, sd_m2=0.0855, driver_type="cone")
        result = prescreen_drivers([big_cone], config)
        assert result.count == 0

    def test_small_cone_passes_high_freq(self, config):
        """Small cone (e.g. 6") should pass when f_piston*factor >= target_f_high."""
        # 6" driver: Sd ~ 130 cm² = 0.0130 m²
        # f_piston = 343/(2π·√(0.0130/π)) = ~849 Hz, ×10 = 8490 > 4000 → passes
        small_cone = _make_driver("small_cone", fs_hz=400, qes=0.4, sd_m2=0.0130, driver_type="cone")
        result = prescreen_drivers([small_cone], config)
        assert result.count == 1

    def test_compression_passes_low_freq(self):
        """Compression drivers should pass low-freq targets (high f_piston)."""
        config = PrescreenConfig(
            target_f_low_hz=100,
            target_f_high_hz=1000,
            mouth_radius_m=0.3,
            length_m=0.8,
        )
        # Sd=20cm²=0.0020 m², f_piston ~2186 Hz, ×10 = 21860 > 1000
        comp = _make_driver("comp", fs_hz=80, qes=0.4, sd_m2=0.0020, driver_type="compression")
        result = prescreen_drivers([comp], config)
        assert result.count == 1

    def test_both_types_pass_mid_band(self):
        """Both types should pass mid-band when Sd is appropriate."""
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=1500,
            mouth_radius_m=0.2,
            length_m=0.5,
        )
        comp = _make_driver("comp", fs_hz=200, qes=0.4, sd_m2=0.0020, driver_type="compression")
        cone = _make_driver("cone", fs_hz=200, qes=0.4, sd_m2=0.0130, driver_type="cone")
        result = prescreen_drivers([comp, cone], config)
        assert result.count == 2

    def test_representative_throat_radius(self, config):
        """Throat radius should be median of sqrt(Sd/pi)."""
        d1 = _make_driver("d1", fs_hz=400, qes=0.4, sd_m2=0.0008)
        d2 = _make_driver("d2", fs_hz=400, qes=0.4, sd_m2=0.0012)
        d3 = _make_driver("d3", fs_hz=400, qes=0.4, sd_m2=0.0010)

        result = prescreen_drivers([d1, d2, d3], config)
        assert result.count == 3

        expected_radii = sorted([math.sqrt(s / math.pi) for s in [0.0008, 0.0012, 0.0010]])
        expected_median = expected_radii[1]  # median of 3
        assert result.throat_radius_m == pytest.approx(expected_median, rel=1e-3)

    def test_sd_ratio_filtering(self, config):
        """Drivers with extreme Sd ratios should be filtered."""
        # Create 3 normal drivers and 1 with very different Sd
        normal = [_make_driver(f"n{i}", fs_hz=400, qes=0.4, sd_m2=0.001) for i in range(3)]
        outlier = _make_driver("outlier", fs_hz=400, qes=0.4, sd_m2=0.05)  # 50x larger

        result = prescreen_drivers(normal + [outlier], config)
        ids = [d.driver_id for d in result.drivers]
        assert "outlier" not in ids

    def test_empty_driver_list(self, config):
        """Empty input should return empty result."""
        result = prescreen_drivers([], config)
        assert result.count == 0
        assert result.throat_radius_m == 0.0
        assert result.drivers == []

    def test_all_filtered_out(self, config):
        """When all drivers fail criteria, return empty result."""
        bad = _make_driver("bad", fs_hz=10000, qes=0.4)
        result = prescreen_drivers([bad], config)
        assert result.count == 0

    def test_to_dict(self, config):
        """PrescreenResult.to_dict should contain driver IDs."""
        d = _make_driver("test1", fs_hz=400, qes=0.4)
        result = prescreen_drivers([d], config)
        d_dict = result.to_dict()
        assert "drivers" in d_dict
        assert "throat_radius_m" in d_dict
        assert "count" in d_dict
        assert d_dict["count"] == 1
        assert "test1" in d_dict["drivers"]
