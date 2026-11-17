"""Validator checks that protect the database from bad imports."""
import pytest

from horn_drivers.validator import find_duplicate_groups, validate_driver


def record(driver_id, **params):
    base = {"fs_hz": 55.0, "re_ohm": 3.3, "bl_tm": 10.2, "sd_m2": 0.0549, "mms_kg": 0.051,
            "le_h": 0.00044, "qms": 6.6, "qes": 0.56, "qts": 0.52}
    base.update(params)
    return {"driver_id": driver_id, "manufacturer": "Test", "model_name": driver_id,
            "driver_type": "cone", "nominal_diameter": "12in", "parameters": base}


def test_page_copied_across_many_models_is_flagged():
    copies = [record(f"faital-{n}") for n in range(8)]
    variants = [record("a-15", re_ohm=5.0), record("a-15-be", re_ohm=5.0)]
    assert find_duplicate_groups(copies + variants) == [sorted(r["driver_id"] for r in copies)]


def test_small_and_large_driver_sharing_parameters_is_flagged():
    """bc-6md38-16 carried 18-inch woofer values."""
    six = {**record("bc-6md38-16"), "nominal_diameter": "6in"}
    eighteen = {**record("bc-18ps76"), "nominal_diameter": "18in"}
    assert find_duplicate_groups([six, eighteen]) == [["bc-18ps76", "bc-6md38-16"]]


@pytest.mark.parametrize("qts,level", [(0.52, None), (0.54, "warning"), (0.38, "error")])
def test_q_identity_uses_relative_tolerance(qts, level):
    """Implied Qts is 0.4997. The old absolute reciprocal threshold (0.1)
    allowed about 5% here but about 30% for a Qts of 2."""
    result = validate_driver(record("x", qts=qts))
    q_errors = [e for e in result.errors if "Q factor" in e]
    q_warnings = [w for w in result.warnings if "Q factor" in w]
    assert bool(q_errors) == (level == "error")
    assert bool(q_warnings) == (level == "warning")


def test_high_q_driver_inconsistency_is_an_error():
    # Qms=10.31, Qes=5.41 imply Qts 3.548. 3.0 is 15% off but only 0.05 apart
    # in reciprocal terms, which the old check accepted.
    result = validate_driver(record("x", qms=10.31, qes=5.41, qts=3.0))
    assert any("Q factor identity" in e for e in result.errors)


def test_mojibake_and_loader_rejection_are_errors():
    bad = record("x")
    bad["model_name"] = "10NMB420_16Î©"
    assert any("mojibake" in e for e in validate_driver(bad).errors)
    unloadable = record("y", sd_m2=1225.0)
    assert any("Loader rejects" in e for e in validate_driver(unloadable).errors)


def test_implausible_derived_qms_is_warned():
    params = dict(fs_hz=30.3, re_ohm=6.1, bl_tm=12.38, sd_m2=0.0531, mms_kg=0.049, qts=0.37)
    r = record("sica", **params)
    for q in ("qms", "qes"):
        r["parameters"].pop(q)
    assert any("Qms=" in w and "implausible" in w for w in validate_driver(r).warnings)
