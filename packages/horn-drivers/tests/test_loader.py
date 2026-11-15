"""Tests for the driver database loader."""

import json

import pytest

from horn_drivers.loader import load_drivers, load_driver


@pytest.fixture
def v2_db(tmp_path):
    """Create a v2-format driver database."""
    db = {
        "schema_version": 2,
        "drivers": [
            {
                "driver_id": "test1",
                "manufacturer": "Acme",
                "model_name": "T1",
                "driver_type": "compression",
                "parameters": {
                    "fs_hz": 500.0,
                    "re_ohm": 6.0,
                    "bl_tm": 8.0,
                    "sd_m2": 0.0008,
                    "mms_kg": 0.003,
                    "le_h": 0.0005,
                    "qms": 5.0,
                    "qes": 0.4,
                },
            },
            {
                "driver_id": "test2",
                "manufacturer": "Acme",
                "model_name": "T2",
                "parameters": {
                    "fs_hz": 600.0,
                    "re_ohm": 5.5,
                    "bl_tm": 7.0,
                    "sd_m2": 0.001,
                    "mms_kg": 0.002,
                    "le_h": 0.0004,
                    "qms": 4.0,
                    "qes": 0.5,
                },
            },
        ],
    }
    path = tmp_path / "drivers.json"
    path.write_text(json.dumps(db))
    return str(path)


@pytest.fixture
def v1_db(tmp_path):
    """Create a v1-format (legacy) driver database."""
    db = {
        "drv1": {
            "driver_id": "drv1",
            "manufacturer": "Old",
            "model_name": "V1",
            "fs_hz": 500.0,
            "re_ohms": 6.0,
            "bl_tm": 8.0,
            "sd_sq_meters": 0.0008,
            "mms_kg": 0.003,
            "le_mh": 0.5,
            "qms": 5.0,
            "qes": 0.4,
            "xmax_mm": 0.6,
        }
    }
    path = tmp_path / "drivers_v1.json"
    path.write_text(json.dumps(db))
    return str(path)


class TestLoadDrivers:
    def test_load_v2_all_drivers(self, v2_db):
        drivers = load_drivers(v2_db)
        assert len(drivers) == 2
        assert drivers[0].driver_id == "test1"
        assert drivers[1].driver_id == "test2"

    def test_all_fields_populated(self, v2_db):
        drivers = load_drivers(v2_db)
        d = drivers[0]
        assert d.fs_hz == 500.0
        assert d.re_ohm == 6.0
        assert d.bl_tm == 8.0
        assert d.sd_m2 == 0.0008
        assert d.mms_kg == 0.003
        assert d.le_h == 0.0005
        # Derived
        assert d.cms_m_per_n is not None
        assert d.qts is not None
        assert d.rms_kg_per_s is not None

    def test_load_v1_format(self, v1_db):
        drivers = load_drivers(v1_db)
        assert len(drivers) == 1
        d = drivers[0]
        assert d.driver_id == "drv1"
        assert d.re_ohm == 6.0
        assert d.le_h == pytest.approx(0.0005, rel=1e-6)
        assert d.sd_m2 == 0.0008
        assert d.xmax_m == pytest.approx(0.0006, rel=1e-6)


class TestLoadDriver:
    def test_load_by_id(self, v2_db):
        d = load_driver(v2_db, "test2")
        assert d.driver_id == "test2"
        assert d.fs_hz == 600.0

    def test_missing_id_raises(self, v2_db):
        with pytest.raises(KeyError, match="nonexistent"):
            load_driver(v2_db, "nonexistent")


class TestRealDatabase:
    """Smoke test against the actual project database."""

    def test_load_project_drivers(self):
        from pathlib import Path

        db_path = Path(__file__).resolve().parents[3] / "data" / "drivers"
        if not db_path.exists():
            pytest.skip("Project drivers database not found")
        drivers = load_drivers(str(db_path))
        assert len(drivers) >= 2
        for d in drivers:
            assert d.fs_hz > 0
            assert d.re_ohm > 0
            assert d.sd_m2 > 0
            assert d.mms_kg > 0
            assert d.cms_m_per_n is not None


def test_invalid_optional_limit_is_rejected_without_discarding_valid_drivers(v2_db, tmp_path):
    from pathlib import Path
    raw = json.loads(Path(v2_db).read_text())
    raw['drivers'][0]['parameters']['power_w'] = float('nan')
    path = tmp_path/'mixed.json'
    path.write_text(json.dumps(raw))
    invalid_id = raw['drivers'][0]['driver_id']
    with pytest.warns(RuntimeWarning, match='Rejected driver '+invalid_id+'.*power_w'):
        drivers = load_drivers(str(path))
    assert [driver.driver_id for driver in drivers] == [raw['drivers'][1]['driver_id']]
    with pytest.raises(ValueError, match='power_w'):
        load_driver(str(path), invalid_id)


@pytest.mark.parametrize("primary,alternate", [("exit_area_m2", "exit_area_cm2"), ("xmax_m", "xmax_mm"), ("sd_m2", "sd_sq_meters"), ("re_ohm", "re_ohms")])
@pytest.mark.parametrize("has_alternate", [False, True])
def test_explicit_zero_is_not_replaced_by_missing_or_alternate_units(v2_db, primary, alternate, has_alternate):
    from pathlib import Path
    from horn_drivers.loader import _driver_from_dict
    raw = json.loads(Path(v2_db).read_text())['drivers'][0]
    raw['parameters'][primary] = 0.
    if has_alternate:
        raw['parameters'][alternate] = 1.
    with pytest.raises(ValueError, match=primary):
        _driver_from_dict(raw)


def test_explicit_zero_inductance_takes_precedence_over_alternate_units(v2_db):
    from pathlib import Path
    from horn_drivers.loader import _driver_from_dict
    raw = json.loads(Path(v2_db).read_text())['drivers'][0]
    raw['parameters'].update(le_h=0., le_mh=5.)
    assert _driver_from_dict(raw).le_h == 0.


# ---------------------------------------------------------------------------
# Datasheet units, unit mistakes and Q consistency
# ---------------------------------------------------------------------------

def _record(**params):
    base = {"fs_hz": 45.0, "re_ohm": 5.6, "bl_tm": 15.2, "le_mh": 0.9, "qms": 5.1, "qes": 0.38}
    base.update(params)
    return {"driver_id": "acme-ac10", "manufacturer": "Acme", "model_name": "AC10", "parameters": base}


def test_datasheet_units_with_vas_instead_of_mms():
    """A typical datasheet gives Sd in cm^2 and Vas in litres, not Mms."""
    from horn_drivers.loader import _driver_from_dict
    d = _driver_from_dict(_record(sd_cm2=346.0, vas_l=52.0))
    assert d.sd_m2 == pytest.approx(0.0346)
    # Hand calculation: Cms = Vas/(rho c^2 Sd^2) = 0.052/(1.204*343^2*0.0346^2)
    # = 3.067e-4 m/N; Mms = 1/((2 pi 45)^2 Cms) = 0.04079 kg.
    assert d.mms_kg == pytest.approx(0.04079, rel=2e-3)
    assert d.fs_hz == pytest.approx(1 / (2 * 3.141592653589793 * (d.mms_kg * d.cms_m_per_n) ** 0.5))


def test_brief_table_aliases_are_accepted():
    from horn_drivers.loader import _driver_from_dict
    d = _driver_from_dict(_record(sd_sq_cm=1225.0, vas_liters=180.0, power_aes_watts=600.0, fs_hz=38.0))
    assert d.sd_m2 == pytest.approx(0.1225)
    assert d.power_w == 600.0
    assert 0.05 < d.mms_kg < 0.5


@pytest.mark.parametrize("field,value,hint", [
    ("sd_m2", 1225.0, "sd_cm2"),       # cm^2 typed under the m^2 key
    ("mms_kg", 180.0, "mms_g"),        # grams under the kg key
    ("le_h", 1.9, "le_mh"),            # mH under the H key
])
def test_order_of_magnitude_unit_mistakes_are_rejected_with_a_hint(field, value, hint):
    from horn_drivers.loader import _driver_from_dict
    record = _record(sd_m2=0.0346, mms_kg=0.04)
    record["parameters"].pop("le_mh")
    record["parameters"].setdefault("le_h", 0.0009)
    record["parameters"][field] = value
    with pytest.raises(ValueError, match=hint):
        _driver_from_dict(record)


def test_inconsistent_three_q_record_is_rejected():
    """lavoce-waf102-50 shipped with Qts 0.38 while Qms 4.01 and Qes 0.12 imply 0.117."""
    from horn_drivers.loader import _driver_from_dict
    record = _record(sd_m2=0.0346, mms_kg=0.04, qms=4.01, qes=0.12, qts=0.38)
    with pytest.raises(ValueError, match="Q factors are inconsistent"):
        _driver_from_dict(record)
    # Two-decimal rounding on a consistent driver is accepted.
    record["parameters"].update(qms=6.6, qes=0.56, qts=0.52)
    assert _driver_from_dict(record).qts == 0.52


def test_supplied_rms_fixes_qms_instead_of_ill_conditioned_derivation():
    """With only Qts given, Qms = Qes*Qts/(Qes-Qts) explodes when the derived
    Qes is close to Qts (sica-12cx3pl8-8 came out at Qms ~ 107), or the driver
    is rejected outright. A supplied Rms must be used instead."""
    from horn_drivers.loader import _driver_from_dict
    params = {"fs_hz": 30.3, "re_ohm": 6.1, "bl_tm": 12.38, "sd_m2": 0.0531, "mms_kg": 0.049,
              "le_h": 0.001, "qts": 0.37}
    with pytest.warns(RuntimeWarning, match="implausibly high"):
        ill = _driver_from_dict({"driver_id": "x", "parameters": dict(params)})
    assert ill.qms > 30
    d = _driver_from_dict({"driver_id": "x", "parameters": {**params, "rms_kg_per_s": 1.5}})
    omega = 2 * 3.141592653589793 * 30.3
    assert d.rms_kg_per_s == 1.5
    assert d.qms == pytest.approx(omega * 0.049 / 1.5)
    assert 1 / d.qts == pytest.approx(1 / d.qms + 1 / d.qes)
