"""Validation tests for the project driver database.

These tests run the validator against data/drivers.json to ensure
every driver has essential T-S parameters, values fall within plausible
physical ranges, and cross-parameter physics checks pass.
"""

import math
from pathlib import Path

import pytest

from horn_drivers.validator import (
    ESSENTIAL_FIELDS,
    RANGES,
    RANGES_BY_TYPE,
    RANGES_BY_DIAMETER,
    ValidationResult,
    compute_derived_params,
    validate_driver,
)
from horn_drivers.loader import load_drivers

DB_PATH = Path(__file__).resolve().parents[3] / "data" / "drivers"


@pytest.fixture
def db_raw():
    """Load the raw driver database (as dicts, not DriverParameters)."""
    if not DB_PATH.exists():
        pytest.skip("Project drivers database not found")
    from horn_drivers.loader import load_drivers_raw
    drivers = load_drivers_raw(str(DB_PATH))
    if not drivers:
        pytest.skip("Driver database is empty")
    return drivers


@pytest.fixture
def db_drivers():
    """Load drivers as DriverParameters objects."""
    if not DB_PATH.exists():
        pytest.skip("Project drivers database not found")
    drivers = load_drivers(str(DB_PATH))
    if not drivers:
        pytest.skip("Driver database is empty")
    return drivers


# -----------------------------------------------------------------------
# Completeness
# -----------------------------------------------------------------------

class TestDriverDatabaseCompleteness:
    """Every driver must have all essential T-S fields populated and non-zero."""

    def test_all_drivers_have_essential_fields(self, db_raw):
        for driver in db_raw:
            result = validate_driver(driver)
            assert result.errors == [], (
                f"[{result.driver_id}] Validation errors: {result.errors}"
            )

    def test_minimum_driver_count(self, db_raw):
        assert len(db_raw) >= 30, (
            f"Database should have at least 30 drivers, got {len(db_raw)}"
        )


# -----------------------------------------------------------------------
# Plausibility — per-type range checks
# -----------------------------------------------------------------------

class TestDriverDatabasePlausibility:
    """All parameter values must be within physically plausible ranges."""

    def test_fs_in_range(self, db_drivers):
        for d in db_drivers:
            ranges = RANGES_BY_TYPE.get(d.driver_type, RANGES)
            lo, hi = ranges.get("fs_hz", RANGES["fs_hz"])
            assert lo <= d.fs_hz <= hi, (
                f"{d.driver_id}: fs_hz={d.fs_hz} outside [{lo}, {hi}]"
            )

    def test_re_in_range(self, db_drivers):
        for d in db_drivers:
            ranges = RANGES_BY_TYPE.get(d.driver_type, RANGES)
            lo, hi = ranges.get("re_ohm", RANGES["re_ohm"])
            assert lo <= d.re_ohm <= hi, (
                f"{d.driver_id}: re_ohm={d.re_ohm} outside [{lo}, {hi}]"
            )

    def test_bl_in_range(self, db_drivers):
        for d in db_drivers:
            ranges = RANGES_BY_TYPE.get(d.driver_type, RANGES)
            lo, hi = ranges.get("bl_tm", RANGES["bl_tm"])
            assert lo <= d.bl_tm <= hi, (
                f"{d.driver_id}: bl_tm={d.bl_tm} outside [{lo}, {hi}]"
            )

    def test_sd_in_range(self, db_drivers):
        for d in db_drivers:
            ranges = RANGES_BY_TYPE.get(d.driver_type, RANGES)
            lo, hi = ranges.get("sd_m2", RANGES["sd_m2"])
            assert lo <= d.sd_m2 <= hi, (
                f"{d.driver_id}: sd_m2={d.sd_m2} outside [{lo}, {hi}]"
            )

    def test_mms_in_range(self, db_drivers):
        for d in db_drivers:
            ranges = RANGES_BY_TYPE.get(d.driver_type, RANGES)
            lo, hi = ranges.get("mms_kg", RANGES["mms_kg"])
            assert lo <= d.mms_kg <= hi, (
                f"{d.driver_id}: mms_kg={d.mms_kg} outside [{lo}, {hi}]"
            )

    def test_le_in_range(self, db_drivers):
        for d in db_drivers:
            ranges = RANGES_BY_TYPE.get(d.driver_type, RANGES)
            lo, hi = ranges.get("le_h", RANGES["le_h"])
            assert lo <= d.le_h <= hi, (
                f"{d.driver_id}: le_h={d.le_h} outside [{lo}, {hi}]"
            )


# -----------------------------------------------------------------------
# Derived parameters
# -----------------------------------------------------------------------

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


# -----------------------------------------------------------------------
# Physics consistency (Tier 2 cross-checks)
# -----------------------------------------------------------------------

class TestPhysicsConsistency:
    """Cross-parameter physics checks on the real database."""

    def test_qes_cross_check(self, db_raw):
        """Qes computed from BL, Mms, Re, fs should be within 15% of stated.

        Compression drivers are excluded — their published Qes often accounts
        for horn throat loading which the free-air formula doesn't capture.
        """
        for driver in db_raw:
            if driver.get("driver_type") == "compression":
                continue
            params = driver.get("parameters", driver)
            qes = params.get("qes")
            if qes is None or qes <= 0:
                continue

            fs = params.get("fs_hz", 0)
            mms = params.get("mms_kg", 0)
            re = params.get("re_ohm", 0)
            bl = params.get("bl_tm", 0)
            if not all(v > 0 for v in [fs, mms, re, bl]):
                continue

            omega = 2.0 * math.pi * fs
            qes_comp = (omega * mms * re) / (bl ** 2)
            ratio = qes_comp / qes
            assert 0.85 <= ratio <= 1.15, (
                f"[{driver['driver_id']}] Qes cross-check: "
                f"stated={qes:.3f}, computed={qes_comp:.3f}, ratio={ratio:.2f}"
            )

    def test_dual_efficiency_agreement(self, db_raw):
        """Two independent eta0 formulas must agree within 1 dB.

        Compression drivers are excluded — these formulas assume
        direct-radiator free-air conditions.
        """
        for driver in db_raw:
            if driver.get("driver_type") == "compression":
                continue
            params = driver.get("parameters", driver)
            derived = compute_derived_params(params)
            if "eta0_A" not in derived or "eta0_B" not in derived:
                continue
            if derived["eta0_A"] <= 0 or derived["eta0_B"] <= 0:
                continue

            spl_a = 10.0 * math.log10(derived["eta0_A"])
            spl_b = 10.0 * math.log10(derived["eta0_B"])
            delta = abs(spl_a - spl_b)
            assert delta <= 1.0, (
                f"[{driver['driver_id']}] Efficiency mismatch: "
                f"method A={spl_a + 112.2:.1f} dB, "
                f"method B={spl_b + 112.2:.1f} dB, "
                f"delta={delta:.2f} dB"
            )


# -----------------------------------------------------------------------
# Derived parameter computation
# -----------------------------------------------------------------------

class TestDerivedParameters:
    """Test that derived parameter computation produces sensible values."""

    def test_vas_computation(self, db_raw):
        """Vas must be positive for all drivers with valid params."""
        for driver in db_raw:
            params = driver.get("parameters", driver)
            derived = compute_derived_params(params)
            if "Vas_m3" in derived:
                assert derived["Vas_m3"] > 0, (
                    f"[{driver['driver_id']}] Vas should be positive"
                )

    def test_eta0_computation(self, db_raw):
        """Efficiency eta0 must be between 0 and 1 (0-100%)."""
        for driver in db_raw:
            params = driver.get("parameters", driver)
            derived = compute_derived_params(params)
            if "eta0_A" in derived:
                assert 0 < derived["eta0_A"] < 1.0, (
                    f"[{driver['driver_id']}] eta0={derived['eta0_A']:.4f} "
                    f"should be 0 < eta0 < 1"
                )

    def test_sensitivity_range(self, db_raw):
        """SPL_1W should be within expected range.

        Cone drivers: 85-105 dB (direct radiator).
        Compression drivers: lower free-air SPL is expected (75-105 dB)
        since they are designed for horn loading.
        """
        for driver in db_raw:
            params = driver.get("parameters", driver)
            derived = compute_derived_params(params)
            if "SPL_1W" not in derived:
                continue
            if driver.get("driver_type") == "compression":
                lo, hi = 75.0, 105.0
            else:
                lo, hi = 85.0, 105.0
            assert lo <= derived["SPL_1W"] <= hi, (
                f"[{driver['driver_id']}] SPL_1W={derived['SPL_1W']:.1f} "
                f"dB outside [{lo}, {hi}]"
            )

    def test_ebp_for_cone_drivers(self, db_raw):
        """Cone drivers intended for horn loading should have EBP > 50."""
        for driver in db_raw:
            if driver.get("driver_type") != "cone":
                continue
            params = driver.get("parameters", driver)
            derived = compute_derived_params(params)
            if "EBP" in derived:
                assert derived["EBP"] > 50.0, (
                    f"[{driver['driver_id']}] EBP={derived['EBP']:.1f}, "
                    f"expected > 50 for horn-suitable drivers"
                )


# -----------------------------------------------------------------------
# Database coverage
# -----------------------------------------------------------------------

class TestDatabaseCoverage:
    """Database must have adequate coverage across sizes and manufacturers."""

    def test_at_least_3_per_size_category(self, db_raw):
        """At least 3 drivers for each nominal diameter."""
        from collections import Counter
        sizes = Counter()
        for driver in db_raw:
            diameter = driver.get("nominal_diameter")
            if diameter:
                sizes[diameter] += 1
        for size in ["6in", "8in", "10in", "12in", "15in", "18in"]:
            assert sizes.get(size, 0) >= 3, (
                f"Need at least 3 drivers for {size}, got {sizes.get(size, 0)}"
            )

    def test_at_least_5_manufacturers(self, db_raw):
        """Database should represent at least 5 distinct manufacturers."""
        manufacturers = {d.get("manufacturer", "").lower() for d in db_raw}
        manufacturers.discard("")
        assert len(manufacturers) >= 5, (
            f"Need at least 5 manufacturers, got {len(manufacturers)}: "
            f"{manufacturers}"
        )

    def test_has_both_driver_types(self, db_raw):
        """Database must contain both compression and cone drivers."""
        types = {d.get("driver_type") for d in db_raw}
        assert "compression" in types, "No compression drivers in database"
        assert "cone" in types, "No cone drivers in database"

    def test_compression_driver_count(self, db_raw):
        """Should have at least 3 compression drivers."""
        count = sum(1 for d in db_raw if d.get("driver_type") == "compression")
        assert count >= 3, f"Need at least 3 compression drivers, got {count}"
