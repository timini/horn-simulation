"""Validation tests for the project driver database.

These tests run the validator against data/drivers/ to ensure the
database has good coverage, the vast majority of drivers have essential
T-S parameters, values fall within plausible physical ranges, and
cross-parameter physics checks pass.

Hard requirements: the validator reports zero errors on the shipped
databases, no record is a copy of another page or a synthetic seed, every
record carries provenance, and the loader rejects nothing. Plausibility
warnings on scraped values use statistical thresholds.
"""

import json
import math
import warnings
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
    find_duplicate_groups,
    validate_driver,
)
from horn_drivers.loader import load_drivers, load_drivers_raw

REPO = Path(__file__).resolve().parents[3]
DB_PATH = REPO / "data" / "drivers"
CURATED_PATH = REPO / "data" / "drivers-curated"
SHIPPED_DBS = [DB_PATH, CURATED_PATH]
SEED_FIXTURE = Path(__file__).parent / "fixtures" / "synthetic_seed_drivers.json"


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


@pytest.mark.parametrize("db", SHIPPED_DBS, ids=lambda p: p.name)
class TestShippedDatabaseHygiene:
    """What CI also enforces with horn-validate-drivers, as direct checks."""

    def test_validator_reports_no_errors(self, db):
        from horn_drivers.validator import validate_database
        errors = [(r.driver_id, e) for r in validate_database(str(db)) for e in r.errors]
        assert errors == [], errors[:10]

    def test_loader_rejects_no_record(self, db):
        raw = load_drivers_raw(str(db))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = load_drivers(str(db))
        rejected = [str(w.message) for w in caught if "Rejected driver" in str(w.message)]
        assert rejected == [] and len(loaded) == len(raw), rejected[:10]

    def test_no_copied_pages(self, db):
        """A wedged origin once served one page for every URL: 376 records
        shared three parameter tuples. Distinct tuples must track the record
        count, with only small variant families sharing values."""
        raw = load_drivers_raw(str(db))
        assert find_duplicate_groups(raw) == []
        keys = {tuple(d["parameters"].get(k) for k in ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg"))
                for d in raw}
        assert len(keys) >= 0.97 * len(raw), f"{len(keys)} distinct tuples for {len(raw)} records"

    def test_no_synthetic_seed_records(self, db):
        seeds = json.loads(SEED_FIXTURE.read_text())["drivers"]
        seed_params = {json.dumps(s["parameters"], sort_keys=True) for s in seeds}
        copies = [d["driver_id"] for d in load_drivers_raw(str(db))
                  if json.dumps(d["parameters"], sort_keys=True) in seed_params]
        assert copies == []

    def test_every_record_has_provenance(self, db):
        """Version 1: a source URL and date. Version 2: a valid provenance entry per stored field."""
        from horn_drivers.records import check_structure
        missing = []
        for d in load_drivers_raw(str(db)):
            provenance = d.get("provenance") or {}
            if d.get("record_schema_version") == 2:
                errors, _ = check_structure(d)
                if errors or set(d["parameters"]) - set(d["field_provenance"]):
                    missing.append(d["driver_id"])
            elif not (str(d.get("parameter_source", "")).startswith("http")
                      and provenance.get("retrieved_on")
                      and d.get("record_schema_version") == 1):
                missing.append(d["driver_id"])
        assert missing == [], missing[:10]

    def test_no_mojibake_and_one_manufacturer_name_per_folder(self, db):
        names = {}
        for path in db.glob("*/*.json"):
            d = json.loads(path.read_text())
            assert "\u00ce\u00a9" not in d["model_name"] and "\u00c3" not in d["model_name"], d["model_name"]
            names.setdefault(path.parent.name, set()).add(d["manufacturer"])
        assert {folder: n for folder, n in names.items() if len(n) > 1} == {}

    def test_driver_categories_are_sourced(self, db):
        """Sd alone once labelled small cones 'compression'."""
        unsourced = [d["driver_id"] for d in load_drivers_raw(str(db))
                     if d.get("driver_type") == "compression" and not d.get("type_source")]
        assert unsourced == []


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

    def test_has_cone_drivers(self, db_raw):
        types = {d.get("driver_type") for d in db_raw}
        assert "cone" in types, "No cone drivers in database"


def test_curated_compression_drivers_are_sourced_with_exit_area():
    """The scraped catalogue has no compression drivers (the source site lists
    woofers only), so the sourced ones live in the curated set. Each must load
    with its datasheet exit area as the throat, not Sd."""
    drivers = [d for d in load_drivers(str(CURATED_PATH)) if d.driver_type == "compression"]
    assert len(drivers) >= 3
    for d in drivers:
        assert d.parameter_source.startswith("http")
        assert d.exit_area_m2 is not None and d.exit_area_m2 < d.sd_m2
        assert d.effective_throat_area == d.exit_area_m2
        assert d.usable_f_low_hz is not None and d.power_w is not None
