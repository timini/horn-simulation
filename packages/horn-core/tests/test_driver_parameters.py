"""Tests for DriverParameters dataclass and derived parameter logic."""

import json
import math

import numpy as np
import pytest

from horn_core.parameters import DriverParameters


class TestParameterDerivation:
    """Verify that derived T-S parameters are computed correctly."""

    def _make_driver(self, **overrides):
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

    def test_cms_derived_from_fs_and_mms(self):
        d = self._make_driver()
        omega_s = 2.0 * np.pi * d.fs_hz
        expected_cms = 1.0 / (d.mms_kg * omega_s ** 2)
        assert d.cms_m_per_n == pytest.approx(expected_cms, rel=1e-9)

    def test_qts_derived_from_qms_and_qes(self):
        d = self._make_driver(qms=5.0, qes=0.45, qts=None)
        expected_qts = (5.0 * 0.45) / (5.0 + 0.45)
        assert d.qts == pytest.approx(expected_qts, rel=1e-9)

    def test_qes_derived_from_qms_and_qts(self):
        qms, qts = 5.0, 0.4
        d = self._make_driver(qms=qms, qes=None, qts=qts)
        expected_qes = (qms * qts) / (qms - qts)
        assert d.qes == pytest.approx(expected_qes, rel=1e-9)

    def test_qms_derived_from_qes_and_qts(self):
        qes, qts = 0.45, 0.4
        d = self._make_driver(qms=None, qes=qes, qts=qts)
        expected_qms = (qes * qts) / (qes - qts)
        assert d.qms == pytest.approx(expected_qms, rel=1e-9)

    def test_rms_derived_from_qms(self):
        d = self._make_driver(qms=5.0)
        omega_s = 2.0 * np.pi * d.fs_hz
        expected_rms = omega_s * d.mms_kg / 5.0
        assert d.rms_kg_per_s == pytest.approx(expected_rms, rel=1e-9)

    def test_q_roundtrip(self):
        """Given Qms and Qes, derive Qts, then verify 1/Qts = 1/Qms + 1/Qes."""
        d = self._make_driver(qms=6.0, qes=0.5)
        reciprocal_sum = 1.0 / d.qms + 1.0 / d.qes
        assert 1.0 / d.qts == pytest.approx(reciprocal_sum, rel=1e-9)


class TestQDerivationFromElectricalParams:
    """Verify Qes derivation from BL, Mms, Re when Q factors not given."""

    def test_qes_derived_when_no_q_given(self):
        d = DriverParameters(
            driver_id="test",
            manufacturer="Test",
            model_name="T1",
            fs_hz=500.0,
            re_ohm=6.0,
            bl_tm=8.5,
            sd_m2=0.0008,
            mms_kg=0.003,
            le_h=0.0006,
        )
        omega_s = 2.0 * np.pi * 500.0
        expected_qes = (omega_s * 0.003 * 6.0) / (8.5 ** 2)
        assert d.qes == pytest.approx(expected_qes, rel=1e-9)


@pytest.mark.parametrize("field", ["qms", "qes", "qts", "cms_m_per_n", "exit_area_m2", "xmax_m", "nominal_impedance_ohm", "power_w", "peak_power_w", "usable_f_low_hz", "usable_f_high_hz", "mmd_kg", "rms_kg_per_s", "rear_load_mass_kg"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1.])
def test_invalid_optional_numbers_are_rejected(field, value):
    with pytest.raises(ValueError, match=field):
        DriverParameters("test", "Test", "Motor", 100., 6., 5., .002, .004, .0001, **{field: value})


@pytest.mark.parametrize("limits", [{"usable_f_low_hz": 1000., "usable_f_high_hz": 100.}, {"usable_f_low_hz": 100., "usable_f_high_hz": 100.}, {"power_w": 0.}, {"qms": 1., "qts": 1.}, {"qes": 1., "qts": 2.}])
def test_invalid_limit_order_and_zero_limits_are_rejected(limits):
    with pytest.raises(ValueError):
        DriverParameters("test", "Test", "Motor", 100., 6., 5., .002, .004, .0001, **limits)


def test_zero_mechanical_loss_and_rear_load_remain_valid_limiting_cases():
    driver = DriverParameters("test", "Test", "Motor", 100., 6., 5., .002, .004, 0., rms_kg_per_s=0., mmd_kg=.003, rear_load_mass_kg=0.)
    assert driver.coupled_moving_mass_kg == .003
