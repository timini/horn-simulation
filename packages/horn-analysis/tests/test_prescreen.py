"""Tests for driver pre-screening module."""

import math

import pytest
import numpy as np

from horn_core.parameters import DriverParameters
from horn_analysis.prescreen import (
    PrescreenConfig,
    PrescreenResult,
    prescreen_drivers,
    _parse_diameter_inches,
)

_SPEED_OF_SOUND = 343.0


def _make_driver(
    driver_id="test",
    fs_hz=500.0,
    qes=0.45,
    qms=5.0,
    sd_m2=0.0008,
    driver_type="compression",
    nominal_diameter=None,
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
        nominal_diameter=nominal_diameter,
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
        """Drivers with suitable fs and EBP should pass."""
        good = _make_driver("good", fs_hz=400, qes=0.4, driver_type="compression")
        result = prescreen_drivers([good], config)
        assert result.count == 1
        assert result.drivers[0].driver_id == "good"

    def test_fs_too_high_filtered(self, config):
        """Driver with fs above target_f_low * 1.5 should be filtered."""
        bad = _make_driver("bad_fs", fs_hz=800, qes=0.4)
        result = prescreen_drivers([bad], config)
        assert result.count == 0

    def test_ebp_too_low_filtered(self, config):
        """Driver with EBP below threshold should be filtered."""
        bad = _make_driver("bad_ebp", fs_hz=400, qes=10.0)
        result = prescreen_drivers([bad], config)
        assert result.count == 0

    def test_driver_larger_than_mouth_filtered(self, config):
        """Driver with effective radius > mouth_radius should be filtered.

        config has mouth_radius=0.2m.
        18" cone: sqrt(0.1/π) = 0.178m < 0.2m → passes.
        Giant driver: sqrt(0.15/π) = 0.219m > 0.2m → filtered.
        """
        giant = _make_driver("giant", fs_hz=400, qes=0.4, sd_m2=0.15, driver_type="cone")
        result = prescreen_drivers([giant], config)
        assert result.count == 0

    def test_large_cone_now_passes_with_decoupled_throat(self, config):
        """15" cone now passes since throat is decoupled from driver.

        Old behavior: max_throat = 0.2/2.0 = 0.1m, 15" filtered.
        New behavior: max_driver = 0.2m (mouth), 15" cone radius=0.165m < 0.2m → passes.
        """
        big_cone = _make_driver("big_cone", fs_hz=400, qes=0.4, sd_m2=0.0855, driver_type="cone")
        result = prescreen_drivers([big_cone], config)
        assert result.count == 1

    def test_small_cone_passes(self, config):
        """6" cone should pass: fits in mouth easily."""
        small_cone = _make_driver("small_cone", fs_hz=400, qes=0.4, sd_m2=0.0130, driver_type="cone")
        result = prescreen_drivers([small_cone], config)
        assert result.count == 1

    def test_8inch_cone_passes_low_freq_horn(self):
        """8" cone should pass for 200-2kHz horn with mouth_radius=0.25m."""
        config = PrescreenConfig(
            target_f_low_hz=200,
            target_f_high_hz=2000,
            mouth_radius_m=0.25,
            length_m=0.4,
        )
        cone_8 = _make_driver("cone_8", fs_hz=100, qes=0.4, sd_m2=0.0214, driver_type="cone")
        result = prescreen_drivers([cone_8], config)
        assert result.count == 1

    def test_8inch_cone_passes_small_mouth(self):
        """8" cone now passes even with smaller mouth (decoupled throat).

        mouth_radius=0.15m → max_driver = 0.15m.
        8" cone: sqrt(0.0214/π) = 0.0825m < 0.15m → passes.
        """
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=7000,
            mouth_radius_m=0.15,
            length_m=0.3,
        )
        cone_8 = _make_driver("cone_8", fs_hz=200, qes=0.4, sd_m2=0.0214, driver_type="cone")
        result = prescreen_drivers([cone_8], config)
        assert result.count == 1

    def test_compression_passes_low_freq(self):
        """Compression drivers should pass low-freq targets."""
        config = PrescreenConfig(
            target_f_low_hz=100,
            target_f_high_hz=1000,
            mouth_radius_m=0.3,
            length_m=0.8,
        )
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
        """PrescreenResult.to_dict should contain driver IDs and throat_radii_m."""
        d = _make_driver("test1", fs_hz=400, qes=0.4)
        result = prescreen_drivers([d], config)
        d_dict = result.to_dict()
        assert "drivers" in d_dict
        assert "throat_radius_m" in d_dict
        assert "throat_radii_m" in d_dict
        assert "count" in d_dict
        assert d_dict["count"] == 1
        assert "test1" in d_dict["drivers"]

    def test_10inch_cone_passes_with_adequate_mouth(self):
        """10" cone passes when mouth is large enough (decoupled throat).

        mouth_radius=0.2m → max_driver = 0.2m.
        10" cone: sqrt(0.0346/π) = 0.105m < 0.2m → passes.
        """
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=7000,
            mouth_radius_m=0.2,
            length_m=0.3,
        )
        cone_10 = _make_driver("cone_10", fs_hz=50, qes=0.35, sd_m2=0.0346, driver_type="cone")
        result = prescreen_drivers([cone_10], config)
        assert result.count == 1

    def test_compression_driver_small_exit_passes(self, config):
        """Compression driver with small exit area passes."""
        comp = _make_driver("comp_small", fs_hz=400, qes=0.4, sd_m2=0.0020,
                            driver_type="compression")
        comp.exit_area_m2 = 0.001
        result = prescreen_drivers([comp], config)
        assert result.count == 1

    def test_6inch_cone_passes_mid_freq(self):
        """6" cone should pass for a mid-frequency target with adequate mouth."""
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=3000,
            mouth_radius_m=0.2,
            length_m=0.4,
        )
        cone_6 = _make_driver("cone_6", fs_hz=80, qes=0.4, sd_m2=0.0130, driver_type="cone")
        result = prescreen_drivers([cone_6], config)
        assert result.count == 1

    def test_8inch_cone_passes_low_freq_large_mouth(self):
        """8" cone passes for low-frequency target with large mouth."""
        config = PrescreenConfig(
            target_f_low_hz=200,
            target_f_high_hz=2000,
            mouth_radius_m=0.30,
            length_m=0.5,
        )
        cone_8 = _make_driver("cone_8", fs_hz=80, qes=0.4, sd_m2=0.0214, driver_type="cone")
        result = prescreen_drivers([cone_8], config)
        assert result.count == 1

    def test_max_driver_derived_from_target_freq(self):
        """When mouth_radius is None, max_driver is derived from ideal mouth at f_low.

        f_low=300 Hz → ideal_mouth = 343/(2π·300) = 0.182m
        6" cone: radius=0.0643m < 0.182m → passes.
        Giant: radius=0.219m > 0.182m → filtered.
        """
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=3000,
            mouth_radius_m=None,
        )
        cone_6 = _make_driver("cone_6", fs_hz=80, qes=0.4, sd_m2=0.0130, driver_type="cone")
        giant = _make_driver("giant", fs_hz=80, qes=0.4, sd_m2=0.15, driver_type="cone")
        result = prescreen_drivers([cone_6, giant], config)
        ids = [d.driver_id for d in result.drivers]
        assert "cone_6" in ids
        assert "giant" not in ids

    def test_effective_throat_area_with_exit_area(self, config):
        """Driver with exit_area_m2 should use it for driver radius calculation."""
        d = _make_driver("phase_plug", fs_hz=400, qes=0.4, sd_m2=0.0214)
        d.exit_area_m2 = 0.001
        result = prescreen_drivers([d], config)
        assert result.count == 1


class TestAcousticThroatRadii:
    """Test that throat radii are derived from acoustic constraints."""

    def test_throat_radii_from_acoustics(self):
        """Throat radii should be based on ka_max at f_high, not driver Sd."""
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=10000,
            mouth_radius_m=0.2,
        )
        # ka_max = 2π, f_high = 10kHz
        # a_acoustic_max = 343 * 2π / (2π * 10000) = 343/10000 = 0.0343m
        a_max = _SPEED_OF_SOUND * config.ka_max / (2 * math.pi * config.target_f_high_hz)

        d = _make_driver("comp", fs_hz=400, qes=0.4, sd_m2=0.0008)
        result = prescreen_drivers([d], config)

        # Acoustic radii at 0.3, 0.65, 1.0 of a_max should be present
        assert len(result.throat_radii_m) >= 3
        assert result.throat_radii_m[-1] == pytest.approx(round(a_max, 6), rel=1e-3)

    def test_throat_radii_include_acoustic_fractions(self, ):
        """throat_radii_m should include 0.3, 0.65, 1.0 fractions of a_acoustic_max."""
        config = PrescreenConfig(
            target_f_low_hz=500,
            target_f_high_hz=4000,
            mouth_radius_m=0.2,
        )
        a_max = _SPEED_OF_SOUND * config.ka_max / (2 * math.pi * config.target_f_high_hz)
        expected_acoustic = [round(a_max * f, 6) for f in [0.3, 0.65, 1.0]]

        d = _make_driver("comp", fs_hz=400, qes=0.4, sd_m2=0.0008)
        result = prescreen_drivers([d], config)

        for r in expected_acoustic:
            assert r in result.throat_radii_m

    def test_throat_radii_include_direct_coupling(self):
        """Small drivers that fit within a_acoustic_max should appear in throat_radii."""
        config = PrescreenConfig(
            target_f_low_hz=500,
            target_f_high_hz=4000,
            mouth_radius_m=0.2,
        )
        a_max = _SPEED_OF_SOUND * config.ka_max / (2 * math.pi * config.target_f_high_hz)

        # Small compression driver: radius well below a_max
        d = _make_driver("small_comp", fs_hz=400, qes=0.4, sd_m2=0.0008)
        drv_radius = round(math.sqrt(0.0008 / math.pi), 6)
        assert drv_radius < a_max  # sanity check

        result = prescreen_drivers([d], config)
        assert drv_radius in result.throat_radii_m

    def test_large_driver_radius_excluded_from_throat_radii(self):
        """Driver radii exceeding a_acoustic_max should not appear in throat_radii."""
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=10000,
            mouth_radius_m=0.2,
        )
        a_max = _SPEED_OF_SOUND * config.ka_max / (2 * math.pi * config.target_f_high_hz)
        # 0.0343m max. 6" cone radius = 0.0643m > a_max
        cone_6 = _make_driver("cone_6", fs_hz=200, qes=0.4, sd_m2=0.0130, driver_type="cone")
        drv_radius = round(math.sqrt(0.0130 / math.pi), 6)
        assert drv_radius > a_max

        result = prescreen_drivers([cone_6], config)
        assert drv_radius not in result.throat_radii_m

    def test_representative_is_middle_of_radii(self):
        """Representative radius should be the middle element of throat_radii_m."""
        config = PrescreenConfig(
            target_f_low_hz=500,
            target_f_high_hz=4000,
            mouth_radius_m=0.2,
        )
        d = _make_driver("comp", fs_hz=400, qes=0.4, sd_m2=0.0008)
        result = prescreen_drivers([d], config)
        mid_idx = len(result.throat_radii_m) // 2
        assert result.throat_radius_m == result.throat_radii_m[mid_idx]

    def test_throat_radii_capped_at_five(self):
        """At most 5 throat radii to limit combinatorial explosion."""
        config = PrescreenConfig(
            target_f_low_hz=500,
            target_f_high_hz=4000,
            mouth_radius_m=0.2,
        )
        # Create many drivers with varied Sd to generate many direct_radii
        drivers = [
            _make_driver(f"d{i}", fs_hz=400, qes=0.4, sd_m2=0.0001 * (i + 1))
            for i in range(20)
        ]
        result = prescreen_drivers(drivers, config)
        assert len(result.throat_radii_m) <= 5


class TestDiameterFilter:
    """Test optional nominal diameter filtering."""

    def test_min_diameter_filter(self):
        """Drivers below min_nominal_diameter_in should be filtered."""
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=3000,
            mouth_radius_m=0.2,
            min_nominal_diameter_in=4.0,
        )
        small = _make_driver("small", fs_hz=200, qes=0.4, sd_m2=0.0020,
                             nominal_diameter="2in")
        big = _make_driver("big", fs_hz=200, qes=0.4, sd_m2=0.0130,
                           nominal_diameter="6in")
        result = prescreen_drivers([small, big], config)
        ids = [d.driver_id for d in result.drivers]
        assert "small" not in ids
        assert "big" in ids

    def test_max_diameter_filter(self):
        """Drivers above max_nominal_diameter_in should be filtered."""
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=3000,
            mouth_radius_m=0.3,
            max_nominal_diameter_in=6.0,
        )
        small = _make_driver("small", fs_hz=200, qes=0.4, sd_m2=0.0020,
                             nominal_diameter="2in")
        big = _make_driver("big", fs_hz=200, qes=0.4, sd_m2=0.0346,
                           nominal_diameter="10in")
        result = prescreen_drivers([small, big], config)
        ids = [d.driver_id for d in result.drivers]
        assert "small" in ids
        assert "big" not in ids

    def test_diameter_filter_skips_missing_diameter(self):
        """Drivers without nominal_diameter should not be filtered by diameter."""
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=3000,
            mouth_radius_m=0.2,
            min_nominal_diameter_in=4.0,
        )
        no_dia = _make_driver("no_dia", fs_hz=200, qes=0.4, sd_m2=0.0020)
        result = prescreen_drivers([no_dia], config)
        assert result.count == 1

    def test_no_diameter_filter_by_default(self):
        """Without diameter config, all drivers pass diameter check."""
        config = PrescreenConfig(
            target_f_low_hz=300,
            target_f_high_hz=3000,
            mouth_radius_m=0.2,
        )
        d = _make_driver("d", fs_hz=200, qes=0.4, sd_m2=0.0020, nominal_diameter="1in")
        result = prescreen_drivers([d], config)
        assert result.count == 1


class TestParseDiameterInches:
    def test_simple_integer(self):
        assert _parse_diameter_inches("4in") == 4.0

    def test_float_value(self):
        assert _parse_diameter_inches("6.5in") == 6.5

    def test_bare_number(self):
        assert _parse_diameter_inches("8") == 8.0

    def test_none_input(self):
        assert _parse_diameter_inches(None) is None

    def test_empty_string(self):
        assert _parse_diameter_inches("") is None

    def test_no_match(self):
        assert _parse_diameter_inches("abc") is None
