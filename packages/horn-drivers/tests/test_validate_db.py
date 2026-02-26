"""Validation tests for the project driver database.

These tests run the validator against data/drivers/ to ensure the
database has good coverage, the vast majority of drivers have essential
T-S parameters, values fall within plausible physical ranges, and
cross-parameter physics checks pass.

Since the database contains scraped data, some drivers may have data
quality issues from the source.  Tests use statistical thresholds
(e.g. "at least 95% pass") rather than requiring perfection.
"""

import math
from collections import Counter
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
from horn_drivers.loader import load_drivers, load_drivers_raw

DB_PATH = Path(__file__).resolve().parents[3] / "data" / "drivers"


@pytest.fixture
def db_raw():
    """Load the raw driver database (as dicts, not DriverParameters)."""
    if not DB_PATH.exists():
        pytest.skip("Project drivers database not found")
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
    """Driver database must have essential T-S fields and good coverage."""

    def test_minimum_driver_count(self, db_raw):
        assert len(db_raw) >= 100, (
            f"Database should have at least 100 drivers, got {len(db_raw)}"
        )

    def test_unique_driver_ids(self, db_drivers):
        ids = [d.driver_id for d in db_drivers]
        assert len(ids) == len(set(ids)), "Duplicate driver IDs found"

    def test_no_validation_errors_above_threshold(self, db_raw):
        """At least 95% of drivers should have zero validation errors."""
        error_free = sum(
            1 for d in db_raw if not validate_driver(d).errors
        )
        rate = error_free / len(db_raw)
        assert rate >= 0.95, (
            f"Only {error_free}/{len(db_raw)} ({rate:.1%}) error-free, "
            f"need >= 95%"
        )


# -----------------------------------------------------------------------
# Plausibility — broad range checks
# -----------------------------------------------------------------------

class TestDriverDatabasePlausibility:
    """Parameter values must be within physically plausible ranges.

    Uses the broadest RANGES (covering all driver types) and requires
    at least 99% pass rate since these are very generous bounds.
    """

    @pytest.mark.parametrize("field", [
        "fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg", "le_h",
    ])
    def test_parameter_in_broad_range(self, db_drivers, field):
        lo, hi = RANGES[field]
        failures = []
        for d in db_drivers:
            val = getattr(d, field)
            if val is not None and val > 0 and not (lo <= val <= hi):
                failures.append(f"{d.driver_id}: {field}={val}")
        rate = 1.0 - len(failures) / len(db_drivers)
        assert rate >= 0.99, (
            f"{field}: {len(failures)}/{len(db_drivers)} outside [{lo}, {hi}] "
            f"({rate:.1%} pass). First failures: {failures[:5]}"
        )


# -----------------------------------------------------------------------
# Derived parameters
# -----------------------------------------------------------------------

class TestDriverDatabaseDerivedParams:
    """Derived parameters must be computed correctly."""

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
        """For drivers with Qms+Qes+Qts, the identity should hold.

        Allows 10% tolerance since scraped data has rounding.
        """
        checked = 0
        failures = []
        for d in db_drivers:
            if d.qms and d.qes and d.qts:
                checked += 1
                reciprocal = 1.0 / d.qms + 1.0 / d.qes
                diff = abs(1.0 / d.qts - reciprocal)
                if diff > 0.1:
                    failures.append(f"{d.driver_id}: diff={diff:.4f}")

        if checked > 0:
            rate = 1.0 - len(failures) / checked
            assert rate >= 0.95, (
                f"Q factor: {len(failures)}/{checked} have diff>0.1 "
                f"({rate:.1%} pass). First: {failures[:5]}"
            )


# -----------------------------------------------------------------------
# Physics consistency (Tier 2 cross-checks)
# -----------------------------------------------------------------------

class TestPhysicsConsistency:
    """Cross-parameter physics checks on the real database."""

    def test_qes_cross_check(self, db_raw):
        """Qes computed from BL, Mms, Re, fs should be within 20% of stated.

        Compression drivers excluded (horn loading affects Qes).
        """
        checked = 0
        failures = []
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

            checked += 1
            omega = 2.0 * math.pi * fs
            qes_comp = (omega * mms * re) / (bl ** 2)
            ratio = qes_comp / qes
            if not (0.80 <= ratio <= 1.20):
                failures.append(
                    f"{driver['driver_id']}: ratio={ratio:.2f}"
                )

        if checked > 0:
            rate = 1.0 - len(failures) / checked
            assert rate >= 0.90, (
                f"Qes cross-check: {len(failures)}/{checked} outside 20% "
                f"({rate:.1%} pass). First: {failures[:5]}"
            )

    def test_dual_efficiency_agreement(self, db_raw):
        """Two independent eta0 formulas should agree within 2 dB.

        Compression drivers excluded.
        """
        checked = 0
        failures = []
        for driver in db_raw:
            if driver.get("driver_type") == "compression":
                continue
            params = driver.get("parameters", driver)
            derived = compute_derived_params(params)
            if "eta0_A" not in derived or "eta0_B" not in derived:
                continue
            if derived["eta0_A"] <= 0 or derived["eta0_B"] <= 0:
                continue

            checked += 1
            spl_a = 10.0 * math.log10(derived["eta0_A"])
            spl_b = 10.0 * math.log10(derived["eta0_B"])
            delta = abs(spl_a - spl_b)
            if delta > 2.0:
                failures.append(
                    f"{driver['driver_id']}: delta={delta:.2f} dB"
                )

        if checked > 0:
            rate = 1.0 - len(failures) / checked
            assert rate >= 0.95, (
                f"Efficiency: {len(failures)}/{checked} delta>2dB "
                f"({rate:.1%} pass). First: {failures[:5]}"
            )


# -----------------------------------------------------------------------
# Derived parameter computation
# -----------------------------------------------------------------------

class TestDerivedParameters:
    """Derived parameter computation produces sensible values."""

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
        failures = []
        for driver in db_raw:
            params = driver.get("parameters", driver)
            derived = compute_derived_params(params)
            if "eta0_A" in derived:
                if not (0 < derived["eta0_A"] < 1.0):
                    failures.append(
                        f"{driver['driver_id']}: eta0={derived['eta0_A']:.4f}"
                    )
        rate = 1.0 - len(failures) / len(db_raw)
        assert rate >= 0.99, (
            f"eta0: {len(failures)} invalid. First: {failures[:5]}"
        )

    def test_sensitivity_range(self, db_raw):
        """SPL_1W should be within 75-110 dB for the vast majority."""
        checked = 0
        failures = []
        for driver in db_raw:
            params = driver.get("parameters", driver)
            derived = compute_derived_params(params)
            if "SPL_1W" not in derived:
                continue
            checked += 1
            if not (75.0 <= derived["SPL_1W"] <= 110.0):
                failures.append(
                    f"{driver['driver_id']}: SPL_1W={derived['SPL_1W']:.1f}"
                )
        if checked > 0:
            rate = 1.0 - len(failures) / checked
            assert rate >= 0.99, (
                f"Sensitivity: {len(failures)}/{checked} outside [75, 110] dB "
                f"({rate:.1%} pass). First: {failures[:5]}"
            )


# -----------------------------------------------------------------------
# Database coverage
# -----------------------------------------------------------------------

class TestDatabaseCoverage:
    """Database must have adequate coverage across sizes and manufacturers."""

    def test_at_least_10_per_size_category(self, db_raw):
        """At least 10 drivers for each nominal diameter."""
        sizes = Counter()
        for driver in db_raw:
            diameter = driver.get("nominal_diameter")
            if diameter:
                sizes[diameter] += 1
        for size in ["6in", "8in", "10in", "12in", "15in", "18in"]:
            assert sizes.get(size, 0) >= 10, (
                f"Need at least 10 drivers for {size}, got {sizes.get(size, 0)}"
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
