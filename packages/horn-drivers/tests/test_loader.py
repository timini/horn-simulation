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
