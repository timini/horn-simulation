"""Validation tests for the project driver database.

These tests run the validator against data/drivers.json to ensure
every driver has essential T-S parameters and that values fall within
plausible physical ranges for compression drivers.
"""

from pathlib import Path

import pytest

from horn_drivers.validator import validate_driver, ESSENTIAL_FIELDS, RANGES
from horn_drivers.loader import load_drivers

DB_PATH = Path(__file__).resolve().parents[3] / "data" / "drivers.json"


@pytest.fixture
def db_raw():
    """Load the raw JSON database (as dicts, not DriverParameters)."""
    import json
    if not DB_PATH.exists():
        pytest.skip("Project drivers.json not found")
    raw = json.loads(DB_PATH.read_text())
    return raw["drivers"] if "drivers" in raw else list(raw.values())


@pytest.fixture
def db_drivers():
    """Load drivers as DriverParameters objects."""
    if not DB_PATH.exists():
        pytest.skip("Project drivers.json not found")
    return load_drivers(str(DB_PATH))


class TestDriverDatabaseCompleteness:
    """Every driver must have all essential T-S fields populated and non-zero."""

    def test_all_drivers_have_essential_fields(self, db_raw):
        for driver in db_raw:
            issues = validate_driver(driver)
            errors = [msg for level, msg in issues if level == "ERROR"]
            assert errors == [], f"Validation errors: {errors}"

    def test_minimum_driver_count(self, db_raw):
        assert len(db_raw) >= 3, "Database should have at least 3 drivers"


class TestDriverDatabasePlausibility:
    """All parameter values must be within physically plausible ranges."""

    def test_fs_in_range(self, db_drivers):
        lo, hi = RANGES["fs_hz"]
        for d in db_drivers:
            assert lo <= d.fs_hz <= hi, (
                f"{d.driver_id}: fs_hz={d.fs_hz} outside [{lo}, {hi}]"
            )

    def test_re_in_range(self, db_drivers):
        lo, hi = RANGES["re_ohm"]
        for d in db_drivers:
            assert lo <= d.re_ohm <= hi, (
                f"{d.driver_id}: re_ohm={d.re_ohm} outside [{lo}, {hi}]"
            )

    def test_bl_in_range(self, db_drivers):
        lo, hi = RANGES["bl_tm"]
        for d in db_drivers:
            assert lo <= d.bl_tm <= hi, (
                f"{d.driver_id}: bl_tm={d.bl_tm} outside [{lo}, {hi}]"
            )

    def test_sd_in_range(self, db_drivers):
        lo, hi = RANGES["sd_m2"]
        for d in db_drivers:
            assert lo <= d.sd_m2 <= hi, (
                f"{d.driver_id}: sd_m2={d.sd_m2} outside [{lo}, {hi}]"
            )

    def test_mms_in_range(self, db_drivers):
        lo, hi = RANGES["mms_kg"]
        for d in db_drivers:
            assert lo <= d.mms_kg <= hi, (
                f"{d.driver_id}: mms_kg={d.mms_kg} outside [{lo}, {hi}]"
            )

    def test_le_in_range(self, db_drivers):
        lo, hi = RANGES["le_h"]
        for d in db_drivers:
            assert lo <= d.le_h <= hi, (
                f"{d.driver_id}: le_h={d.le_h} outside [{lo}, {hi}]"
            )


class TestDriverDatabaseDerivedParams:
    """Derived parameters must be computed correctly for every driver."""

    def test_all_have_cms(self, db_drivers):
        for d in db_drivers:
            assert d.cms_m_per_n is not None and d.cms_m_per_n > 0, (
                f"{d.driver_id}: missing Cms"
            )

    def test_all_have_qts(self, db_drivers):
        for d in db_drivers:
            assert d.qts is not None and d.qts > 0, (
                f"{d.driver_id}: missing Qts"
            )

    def test_q_factor_consistency(self, db_drivers):
        """1/Qts should equal 1/Qms + 1/Qes."""
        for d in db_drivers:
            if d.qms and d.qes and d.qts:
                reciprocal = 1.0 / d.qms + 1.0 / d.qes
                assert abs(1.0 / d.qts - reciprocal) < 1e-6, (
                    f"{d.driver_id}: Q factor inconsistency"
                )

    def test_unique_driver_ids(self, db_drivers):
        ids = [d.driver_id for d in db_drivers]
        assert len(ids) == len(set(ids)), "Duplicate driver IDs found"
