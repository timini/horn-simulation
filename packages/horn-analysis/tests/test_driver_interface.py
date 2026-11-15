"""Driver interface (front/rear chamber) and electrical outputs, checked against
independent physics rather than the implementation."""
import dataclasses
import math

import numpy as np
import pandas as pd
import pytest

from horn_core.acoustics import C0, RHO0, baffled_piston_on_axis, pressure_level
from horn_core.parameters import DriverParameters
from horn_analysis.evaluation import evaluate_response
from horn_analysis.electrical_output import electrical_table, system_group_delay, stress_assessment
from horn_analysis.scoring import TargetSpec
from horn_analysis.transfer_function import compute_driver_operating_point


def cone(**overrides):
    values = dict(driver_id="cone", manufacturer="Test", model_name="6.5in", fs_hz=110., re_ohm=5.3,
                  bl_tm=9., sd_m2=.013, mms_kg=.0085, le_h=.0001, qms=2.7, qes=.38,
                  driver_type="cone", xmax_m=.003, nominal_impedance_ohm=8., power_w=200.)
    values.update(overrides)
    return DriverParameters(**values)


def ideal_transformer_reference(driver, f, z, area, v_g):
    """The operating point exactly as computed before chambers were added."""
    omega = 2 * np.pi * f
    z_e = driver.re_ohm + 1j * omega * driver.le_h
    rms = driver.rms_kg_per_s if driver.rms_kg_per_s is not None else 0.0
    z_mech_driver = rms + 1j * omega * driver.coupled_moving_mass_kg + 1.0 / (1j * omega * driver.cms_m_per_n)
    z_mech_load = z * (driver.sd_m2 ** 2) / area
    z_mech_total = z_mech_driver + z_mech_load
    denominator = z_e * z_mech_total + driver.bl_tm**2
    current = v_g * z_mech_total / denominator
    velocity = driver.bl_tm * v_g / denominator
    return {"throat_pressure": z * velocity * (driver.sd_m2 / area), "velocity_rms": velocity,
            "current_rms": current, "electrical_impedance": v_g / current,
            "displacement_peak_m": np.sqrt(2) * np.abs(velocity) / omega,
            "input_power_w": np.real(v_g * np.conjugate(current)),
            "copper_power_w": np.abs(current)**2 * driver.re_ohm,
            "mechanical_loss_w": np.abs(velocity)**2 * rms,
            "horn_power_w": np.abs(velocity)**2 * np.real(z_mech_load)}


F = np.geomspace(20., 20000., 400)
Z = 300. + 150j * np.sin(F / 700.)  # arbitrary complex specific throat impedance
AREA = np.pi * .0325**2


@pytest.mark.parametrize("driver", [cone(front_chamber_volume_m3=0.), cone(driver_type=None), cone(driver_type="compression")])
def test_zero_front_chamber_reproduces_ideal_transformer_exactly(driver):
    point = compute_driver_operating_point(driver, F, Z.real, Z.imag, AREA, 2.83)
    reference = ideal_transformer_reference(driver, F, Z, AREA, 2.83)
    for key, expected in reference.items():
        np.testing.assert_array_equal(point[key], expected, err_msg=key)


def test_cone_without_chamber_data_uses_flagged_cone_volume_estimate():
    driver = cone()
    diameter = 2 * math.sqrt(.013 / math.pi)
    expected_volume = .013 * .25 * diameter / 3  # right cone, depth a quarter of the diameter
    point = compute_driver_operating_point(driver, F, Z.real, Z.imag, AREA, 2.83)
    assert point["front_chamber_source"] == "estimated"
    assert point["front_chamber_volume_m3"] == pytest.approx(expected_volume, rel=1e-12)
    assert 1.3e-4 < expected_volume < 1.5e-4  # about 140 cm3 for a 6.5 inch cone
    result = evaluate_response(np.array([800., 1600.]), np.array([100., 100.]), TargetSpec(800, 1600), driver, AREA)
    assert "front_chamber_volume_estimated" in result["evidence_gaps"]
    assert result["driver_interface"]["front_chamber_source"] == "estimated"
    assert result["driver_interface"]["front_chamber_volume_cm3"] == pytest.approx(expected_volume * 1e6)


def test_interface_flags_for_expansion_and_unknown_chamber():
    target = TargetSpec(500, 1000)
    unknown = evaluate_response(np.array([500., 1000.]), np.full(2, 100.), target, cone(driver_type=None), AREA)
    assert "front_chamber_volume_unknown" in unknown["evidence_gaps"]
    expanding = evaluate_response(np.array([500., 1000.]), np.full(2, 100.), target, cone(), .02)
    assert "diaphragm_smaller_than_throat" in expanding["evidence_gaps"]
    # Gaps flag the assumption; they do not reject the design.
    assert unknown["model_feasible"] and expanding["model_feasible"]


def test_front_chamber_corner_on_resistive_load():
    """Shunt compliance C = V/(rho c^2) feeding Z_A = rho c / S_t.

    Volume velocity divides as U_t = U_d / (1 + j w C Z_A), so |U_t/U_d|^2 = 1/2
    where w = 1/(C Z_A) = (rho c^2 / V)(S_t / (rho c)) = c S_t / V, giving
    f = c S_t / (2 pi V).
    """
    volume = 4e-4  # large chamber: corner well inside the sweep
    driver = cone(front_chamber_volume_m3=volume)
    corner = C0 * AREA / (2 * math.pi * volume)
    f = np.geomspace(corner / 20, corner * 20, 20001)
    z = np.full(f.shape, RHO0 * C0)
    point = compute_driver_operating_point(driver, f, z, 0 * z, AREA, 2.83)
    ratio_db = 20 * np.log10(np.abs(point["throat_volume_velocity"] / (driver.sd_m2 * point["velocity_rms"])))
    np.testing.assert_allclose(ratio_db, -10 * np.log10(1 + (f / corner)**2), atol=1e-9)
    crossing = f[np.argmin(np.abs(ratio_db + 10 * np.log10(2)))]
    assert crossing == pytest.approx(corner, rel=1e-3)
    assert ratio_db[-1] < -20  # first-order roll-off: -26 dB a decade and a bit above


@pytest.mark.parametrize("chamber", [dict(front_chamber_volume_m3=2e-4),
                                     dict(front_chamber_volume_m3=2e-4, front_chamber_area_m2=.006),
                                     dict(front_chamber_volume_m3=1e-4, rear_chamber_volume_m3=3e-3)])
def test_energy_balance_with_chambers(chamber):
    driver = cone(**chamber)
    point = compute_driver_operating_point(driver, F, Z.real, Z.imag, AREA, 2.83)
    acoustic = np.real(point["throat_pressure"] * np.conj(point["throat_volume_velocity"]))
    np.testing.assert_allclose(point["horn_power_w"], acoustic, rtol=1e-12)
    np.testing.assert_allclose(point["input_power_w"],
                               point["copper_power_w"] + point["mechanical_loss_w"] + acoustic, rtol=1e-9)


def two_port_impedance(driver, f, z_specific, area, volume, chamber_area=None, rear_volume=None):
    """Input impedance from chained ABCD matrices: Ze, gyrator, Zm, Sd transformer, chamber, horn."""
    out = []
    for frequency, z in zip(f, z_specific):
        w = 2 * np.pi * frequency
        k = w / C0
        series = lambda impedance: np.array([[1, impedance], [0, 1]], complex)
        zm = driver.rms_kg_per_s + 1j * w * driver.mms_kg + 1 / (1j * w * driver.cms_m_per_n)
        if rear_volume:
            zm += driver.sd_m2**2 * RHO0 * C0**2 / (1j * w * rear_volume)
        gyrator = np.array([[0, driver.bl_tm], [1 / driver.bl_tm, 0]], complex)
        transformer = np.array([[driver.sd_m2, 0], [0, 1 / driver.sd_m2]], complex)  # F = Sd p, v = U / Sd
        if chamber_area is None:
            chamber = np.array([[1, 0], [1j * w * volume / (RHO0 * C0**2), 1]], complex)
        else:
            length, zc = volume / chamber_area, RHO0 * C0 / chamber_area
            chamber = np.array([[np.cos(k * length), 1j * zc * np.sin(k * length)],
                                [1j * np.sin(k * length) / zc, np.cos(k * length)]])
        chain = series(driver.re_ohm + 1j * w * driver.le_h) @ gyrator @ series(zm) @ transformer @ chamber
        load = z / area
        out.append((chain[0, 0] * load + chain[0, 1]) / (chain[1, 0] * load + chain[1, 1]))
    return np.array(out)


@pytest.mark.parametrize("chamber_area,rear", [(None, None), (.004, None), (None, 2e-3)])
def test_electrical_impedance_matches_independent_two_port(chamber_area, rear):
    driver = cone(front_chamber_volume_m3=1.5e-4, front_chamber_area_m2=chamber_area, rear_chamber_volume_m3=rear)
    point = compute_driver_operating_point(driver, F, Z.real, Z.imag, AREA, 2.83)
    expected = two_port_impedance(driver, F, Z, AREA, 1.5e-4, chamber_area, rear)
    np.testing.assert_allclose(point["electrical_impedance"], expected, rtol=1e-10)


def test_impedance_peak_at_horn_loaded_resonance():
    """A mass-like horn load and a sealed rear chamber move the |Ze| peak.

    With Le = 0, |Ze| = |Re + Bl^2 / Zm| peaks where Zm is real, at
    w0 = sqrt(k / M) with k = 1/Cms + Sd^2 rho c^2 / Vrc and
    M = Mms + Sd^2 rho l / S_t for a load z = R + j w rho l.
    """
    driver = cone(le_h=0., front_chamber_volume_m3=0., rear_chamber_volume_m3=5e-3)
    air_length = .05
    f = np.geomspace(40, 400, 40001)
    z = 200. + 1j * 2 * np.pi * f * RHO0 * air_length
    point = compute_driver_operating_point(driver, f, z.real, z.imag, AREA, 2.83)
    stiffness = 1 / driver.cms_m_per_n + driver.sd_m2**2 * RHO0 * C0**2 / 5e-3
    mass = driver.mms_kg + driver.sd_m2**2 * RHO0 * air_length / AREA
    loaded = math.sqrt(stiffness / mass) / (2 * math.pi)
    peak = f[np.argmax(np.abs(point["electrical_impedance"]))]
    assert peak == pytest.approx(loaded, rel=1e-4)
    assert abs(loaded - driver.fs_hz) > 10  # the shift is real, not fs by default


@pytest.mark.parametrize("rear", [None, 4e-3])
def test_low_frequency_excursion_limit(rear):
    """Far below resonance the stiffness dominates: x_peak = sqrt(2) Bl I / k, I = V / Re."""
    driver = cone(rear_chamber_volume_m3=rear)
    f = np.array([.5])
    point = compute_driver_operating_point(driver, f, np.array([RHO0 * C0]), np.zeros(1), AREA, 2.83)
    stiffness = 1 / driver.cms_m_per_n + (0 if rear is None else driver.sd_m2**2 * RHO0 * C0**2 / rear)
    expected = math.sqrt(2) * driver.bl_tm * abs(point["current_rms"][0]) / stiffness
    assert point["displacement_peak_m"][0] == pytest.approx(expected, rel=2e-4)
    assert point["displacement_peak_m"][0] == pytest.approx(math.sqrt(2) * driver.bl_tm * 2.83 / driver.re_ohm / stiffness, rel=2e-3)


def test_stress_power_is_separate_from_drive_voltage():
    driver = cone(front_chamber_volume_m3=0.)
    f = np.geomspace(50, 4000, 300)
    z = np.full(f.shape, RHO0 * C0)
    levels = np.full(f.shape, 100.)
    point = compute_driver_operating_point(driver, f, z, 0 * z, AREA, 2.83)
    base = evaluate_response(f, levels, TargetSpec(800, 1600), driver, AREA, operating_point=point)
    stressed = evaluate_response(f, levels, TargetSpec(800, 1600, test_power_w=100.), driver, AREA, operating_point=point)
    # Drive-level results and the score do not change with the stress option.
    for key in ("composite_score", "passband_ripple_db", "avg_sensitivity_db", "max_displacement_peak_m"):
        assert stressed[key] == base[key]
    volts = math.sqrt(100 * 8.)
    excursion = point["displacement_peak_m"] * volts / 2.83
    assert stressed["test_voltage_rms"] == pytest.approx(volts)
    assert stressed["test_power_sweep_max_excursion_mm"] == pytest.approx(excursion.max() * 1e3)
    # Excursion falls with frequency above resonance, so the in-band maximum is at the 800 Hz edge.
    at_band_edge = np.interp(np.log(800.), np.log(f), excursion)
    assert stressed["test_power_xmax_margin_mm"] == pytest.approx(3. - at_band_edge * 1e3, rel=1e-12)
    # 100 W drives this cone past Xmax only below the band; the crossover must sit above that.
    assert stressed["model_feasible"]
    assert 50 < stressed["test_power_xmax_exceeded_below_hz"] < 800
    assert excursion[f <= stressed["test_power_xmax_exceeded_below_hz"]].max() > .003
    loud = evaluate_response(f, levels, TargetSpec(800, 1600, test_power_w=5e4), driver, AREA, operating_point=point)
    assert "test_power_excursion_exceeded" in loud["rejection_reasons"] and loud["composite_score"] == 0


def test_group_delay_removes_time_of_flight():
    f = np.geomspace(100, 10000, 2000)
    w = 2 * np.pi * f
    corner = 2 * np.pi * 1000
    tof = .0035
    pressure = np.exp(-1j * w * tof) / (1 + 1j * w / corner)
    delay = system_group_delay(f, pressure, tof)
    expected = (1 / corner) / (1 + (w / corner)**2)  # first-order low-pass group delay
    np.testing.assert_allclose(delay, expected, rtol=2e-3, atol=2e-7)


def test_electrical_table_columns_and_values():
    f = np.geomspace(300, 3000, 60)
    w = 2 * np.pi * f
    length, distance, mouth_area = .2, 1., .02
    u = 1e-4 * np.exp(-1j * w * length / C0)  # mouth flow per pascal at the throat
    frame = pd.DataFrame(dict(frequency=f, spl=80., z_real=RHO0 * C0, z_imag=0., schema_version=2,
                              phasor_convention="exp(+iwt)_rms", bc_mode="dirichlet", inlet_area_m2=AREA,
                              mouth_area_m2=mouth_area, mouth_u_real=u.real, mouth_u_imag=u.imag))
    driver = cone()
    target = TargetSpec(500, 2000, test_power_w=50.)
    table = electrical_table(frame, driver, target, None, length)
    point = compute_driver_operating_point(driver, f, frame.z_real, frame.z_imag, AREA, 2.83)
    for column in ("frequency", "spl", "ze_real_ohm", "ze_imag_ohm", "ze_mag_ohm", "excursion_peak_mm",
                   "input_power_w", "group_delay_ms", "test_excursion_peak_mm", "test_xmax_margin_mm"):
        assert column in table and np.isfinite(table[column]).all(), column
    np.testing.assert_allclose(table.ze_mag_ohm, np.abs(point["electrical_impedance"]), rtol=1e-12)
    np.testing.assert_allclose(table.excursion_peak_mm, point["displacement_peak_m"] * 1e3, rtol=1e-12)
    observer = baffled_piston_on_axis(f, u * point["throat_pressure"], mouth_area, distance)
    np.testing.assert_allclose(table.spl, pressure_level(observer), atol=1e-9)
    assert table.time_of_flight_removed_ms.iloc[0] == pytest.approx((length + distance) / C0 * 1e3)
    # With the straight-through delay removed, what is left is the driver and chamber delay:
    # small compared with the 3.5 ms time of flight.
    assert np.abs(table.group_delay_ms).max() < 1.
    np.testing.assert_allclose(table.test_xmax_margin_mm, 3. - table.test_excursion_peak_mm, rtol=1e-12)


def test_auto_report_writes_electrical_csv_for_ranked_candidates(tmp_path):
    """Regenerate the report from the preserved 800-1600 Hz modal example, no new solves."""
    import json
    from pathlib import Path
    from horn_drivers.loader import load_driver
    from horn_analysis.auto_report import generate_auto_report

    root = Path(__file__).resolve().parents[3]
    example = root / "examples/800-1600-modal"
    ranked = json.loads((example / "ranking.json").read_text())[:1]
    driver = load_driver(str(root / "data/drivers-curated"), ranked[0]["driver_id"])
    label = ranked[0]["horn_label"]
    target = TargetSpec(800, 1600, test_power_w=50.)
    generate_auto_report(ranked, {label: str(example / "horn-impedance.csv")}, {driver.driver_id: driver},
                         ranked[0]["throat_radius"], target, str(tmp_path), top_n=1)
    table = pd.read_csv(tmp_path / f"coupled_01_{driver.driver_id}_{label}.csv")
    assert list(table.columns[:2]) == ["frequency", "spl"]
    for column in ("ze_real_ohm", "ze_imag_ohm", "ze_mag_ohm", "excursion_peak_mm", "input_power_w",
                   "group_delay_ms", "test_xmax_margin_mm"):
        assert np.isfinite(table[column]).all(), column
    # Above resonance a cone's |Ze| sits above Re; time of flight is the horn length plus 1 m.
    assert (table.ze_mag_ohm > driver.re_ohm).all()
    assert table.time_of_flight_removed_ms.iloc[0] == pytest.approx((ranked[0]["length"] + 1.) / C0 * 1e3)
    assert (tmp_path / "auto_electrical.png").stat().st_size > 0
