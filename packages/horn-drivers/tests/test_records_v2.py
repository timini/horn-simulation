"""Record schema version 2, identity and the loader's reading of both versions."""

import json

import pytest

from horn_drivers.identity import canonical_manufacturer, match_key, normalise_model, record_identity
from horn_drivers.loader import _driver_from_dict, load_drivers
from horn_drivers.records import check_structure, driver_id_for, make_record, parameter_source, write_incomplete, write_record
from horn_drivers.validator import validate_driver

URL = "https://www.bcspeakers.com/en/products/lf-driver/8-0/8/8NDL51"
PARAMS = {"fs_hz": 61.0, "re_ohm": 5.1, "bl_tm": 12.3, "sd_m2": 0.022, "mms_kg": 0.022,
          "le_h": 0.0004, "qms": 6.0, "qes": 0.36, "nominal_impedance_ohm": 8.0,
          "power_w": 200.0, "peak_power_w": 400.0, "sensitivity_db_1w1m": 94.0, "xmax_m": 0.0045}


def provenance(fields=PARAMS, method="page"):
    return {name: {"source_id": "bc", "url": URL, "retrieved_on": "2026-11-17", "method": method} for name in fields}


def v2_record(**overrides):
    record = make_record(manufacturer="B&C Speakers", model="8NDL51", parameters=PARAMS,
                         field_provenance={**provenance(), **provenance(("usable_f_low_hz", "usable_f_high_hz"))},
                         driver_type="cone", aliases=["8NDL51-8"], nominal_diameter="8in",
                         usable_f_low_hz=70.0, usable_f_high_hz=6000.0)
    record.update(overrides)
    return record


def test_manufacturer_aliases_resolve_to_canonical_folder_names():
    assert canonical_manufacturer("Eighteen Sound") == "18Sound"
    assert canonical_manufacturer("eighteensound") == "18Sound"
    assert canonical_manufacturer("B&C Speakers") == "BC"
    assert canonical_manufacturer("Faital Pro") == "Faital"
    assert canonical_manufacturer("Nobody Audio") is None


def test_model_normalisation_strips_brand_and_punctuation():
    assert normalise_model("Eminence Beta-12A", "Eminence") == "BETA12A"
    assert normalise_model("BETA 12A", "Eminence") == "BETA12A"
    assert normalise_model("18 Sound 8NMB420", "18Sound") == "8NMB420"
    assert match_key("eighteen sound", "8NMB420", 8) == ("18Sound", "8NMB420", 8.0)
    with pytest.raises(ValueError):
        match_key("Nobody", "X", 8)


def test_make_record_builds_a_valid_v2_record():
    record = v2_record()
    assert record["record_schema_version"] == 2
    assert record["driver_id"] == "bc-8ndl51"
    assert record["manufacturer"] == record["identity"]["manufacturer"] == "BC"
    assert record["identity"]["variant"]["impedance_ohm"] == 8.0
    assert record["parameter_source"] == URL
    assert check_structure(record) == ([], [])
    assert driver_id_for("BC", "8NDL51", 8) == "bc-8ndl51-8ohm"


def test_make_record_rejects_unknown_manufacturer_and_bad_method():
    with pytest.raises(ValueError, match="Unknown manufacturer"):
        make_record(manufacturer="Nobody", model="X", parameters=PARAMS, field_provenance=provenance())
    bad = provenance()
    bad["fs_hz"]["method"] = "guess"
    with pytest.raises(ValueError, match="method"):
        make_record(manufacturer="BC", model="X", parameters=PARAMS, field_provenance=bad)


def test_structure_check_requires_provenance_for_essential_fields():
    record = v2_record()
    del record["field_provenance"]["bl_tm"]
    errors, _ = check_structure(record)
    assert errors == ["bl_tm has no field_provenance entry"]
    record = v2_record()
    del record["field_provenance"]["xmax_m"]
    errors, warnings = check_structure(record)
    assert errors == [] and warnings == ["xmax_m has no field_provenance entry"]
    record = v2_record()
    record["identity"]["manufacturer"] = "B&C"
    assert any("canonical" in e for e in check_structure(record)[0])
    assert check_structure({"record_schema_version": 3})[0]


def test_loader_reads_v2_fields_sensitivity_power_and_band():
    driver = _driver_from_dict(v2_record())
    assert driver.manufacturer == "BC" and driver.model_name == "8NDL51"
    assert driver.sensitivity_db_1w1m == 94.0
    assert driver.power_w == 200.0 and driver.peak_power_w == 400.0
    assert (driver.usable_f_low_hz, driver.usable_f_high_hz) == (70.0, 6000.0)
    assert driver.parameter_source == URL


def test_loader_takes_identity_when_top_level_names_are_missing():
    record = v2_record()
    del record["manufacturer"], record["model_name"]
    record["parameters"] = {k: v for k, v in PARAMS.items() if k != "nominal_impedance_ohm"}
    driver = _driver_from_dict(record)
    assert (driver.manufacturer, driver.model_name, driver.nominal_impedance_ohm) == ("BC", "8NDL51", 8.0)


def test_peak_power_never_stands_in_for_nominal_power():
    params = {k: v for k, v in PARAMS.items() if k != "power_w"}
    driver = _driver_from_dict({"driver_id": "x", "parameters": params})
    assert driver.power_w is None and driver.peak_power_w == 400.0
    driver = _driver_from_dict({"driver_id": "x", "parameters": {**params, "program_power_w": 500.0, "peak_power_w": None}})
    assert driver.peak_power_w == 500.0 and driver.power_w is None


def test_band_inside_parameters_is_read():
    driver = _driver_from_dict({"driver_id": "x", "parameters": {**PARAMS, "usable_f_low_hz": 80.0, "usable_f_high_hz": 5000.0}})
    assert (driver.usable_f_low_hz, driver.usable_f_high_hz) == (80.0, 5000.0)


def test_implausible_sensitivity_is_rejected():
    with pytest.raises(ValueError, match="sensitivity_db_1w1m"):
        _driver_from_dict({"driver_id": "x", "parameters": {**PARAMS, "sensitivity_db_1w1m": 940.0}})


def test_unsupported_schema_version_is_rejected():
    with pytest.raises(ValueError, match="record_schema_version"):
        _driver_from_dict({"driver_id": "x", "record_schema_version": 9, "parameters": PARAMS})


def test_v1_records_still_load_with_sensitivity(tmp_path):
    v1 = {"driver_id": "old", "manufacturer": "Peerless", "model_name": "Old", "record_schema_version": 1,
          "parameters": PARAMS, "parameter_source": "https://example.com/sheet.pdf"}
    (tmp_path / "Peerless").mkdir()
    (tmp_path / "Peerless" / "old.json").write_text(json.dumps(v1))
    write_record(tmp_path, v2_record())
    drivers = {d.driver_id: d for d in load_drivers(str(tmp_path))}
    assert set(drivers) == {"old", "bc-8ndl51"}
    assert drivers["old"].sensitivity_db_1w1m == 94.0
    assert drivers["old"].parameter_source == "https://example.com/sheet.pdf"
    assert record_identity(v1)["manufacturer"] == "Peerless"


def test_validator_reports_schema_errors_and_sensitivity_gap():
    record = v2_record()
    record["field_provenance"]["fs_hz"]["retrieved_on"] = "yesterday"
    result = validate_driver(record)
    assert any("retrieved_on" in e for e in result.errors)
    loud = v2_record(parameters={**PARAMS, "sensitivity_db_1w1m": 120.0})
    assert any("Sensitivity cross-check" in w for w in validate_driver(loud).warnings)
    assert not any("Sensitivity" in w for w in validate_driver(v2_record()).warnings)


def test_incomplete_writer_keeps_reason(tmp_path):
    path = write_incomplete(tmp_path, v2_record(), "no Bl on the page")
    assert json.loads(path.read_text())["incomplete_reason"] == "no Bl on the page"
    assert parameter_source({"field_provenance": {"qts": {"url": "https://a"}}}) == "https://a"
