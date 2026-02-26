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
