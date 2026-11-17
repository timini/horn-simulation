"""driver-base import, tested offline against a small facts-only extract of the dataset."""

import csv
import json
from pathlib import Path

from horn_drivers.fetch import FetchResult
from horn_drivers.loader import _driver_from_dict
from horn_drivers.sources import get_adapter_class
from horn_drivers.sources.driver_base import (DATASET_URL, DriverBaseSource, compare_records, convert_entry,
                                              fill_missing, prefix_candidate, same_driver_by_values)

FIXTURE = Path(__file__).parent / "fixtures" / "driver_base_sample.json"


def entries():
    return {e["canonical_id"]: e for e in json.loads(FIXTURE.read_text())["drivers"]}


class FixtureFetcher:
    source_id = "driver-base"

    def __init__(self):
        self.calls = []

    def get(self, url, max_age_days=None):
        self.calls.append(url)
        if url != DATASET_URL:
            return FetchResult(url, 404, b"", "text/plain", "2026-11-17T10:00:00+00:00", False)
        return FetchResult(url, 200, FIXTURE.read_bytes(), "text/plain", "2026-11-17T10:00:00+00:00", False)

    def status(self):
        return {"source_id": self.source_id}


def legacy_10cx650():
    """A scraped version 1 record for the same driver, without power, sensitivity or band."""
    return {"driver_id": "18sound-10cx650", "manufacturer": "18Sound", "model_name": "10CX650",
            "driver_type": "cone",
            "parameters": {"fs_hz": 54.0, "re_ohm": 5.1, "bl_tm": 14.0, "sd_m2": 0.0346, "mms_kg": 0.0395,
                           "le_h": 0.0009, "qts": 0.35, "qes": 0.39, "qms": 5.9, "nominal_impedance_ohm": 8.0},
            "parameter_source": "https://loudspeakerdatabase.com/18Sound/10CX650", "record_schema_version": 1,
            "provenance": {"known_limitations": ["Band, power handling, Rms and sensitivity were not captured."]}}


def test_registered_under_the_registry_id():
    assert get_adapter_class("driver-base") is DriverBaseSource


def test_convert_keeps_read_values_in_si_and_drops_derived_ones():
    result = convert_entry(entries()["eighteensound__10cx650__8ohm"])
    record = result.record
    assert result.reason is None and record["record_schema_version"] == 2
    params = record["parameters"]
    assert params["sd_m2"] == 0.0346 and params["mms_kg"] == 0.0395 and params["le_h"] == 0.0009
    assert params["power_w"] == 200.0 and params["sensitivity_db_1w1m"] == 94.0
    assert "peak_power_w" not in params  # program power is derived (2 x AES) upstream
    assert (record["usable_f_low_hz"], record["usable_f_high_hz"]) == (60.0, 4800.0)
    assert record["driver_type"] == "cone" and record["driver_id"] == "18sound-10cx650"
    fs = record["field_provenance"]["fs_hz"]
    assert fs["url"].startswith("https://www.eighteensound.it/") and fs["method"] == "page"
    assert fs["via"] == DATASET_URL and fs["retrieved_on"] == "2026-10-18"
    assert any("low-frequency section" in note for note in record["provenance"]["notes"])
    _driver_from_dict(record)


def test_derived_sensitivity_is_not_imported():
    record = convert_entry(entries()["bcspeakers__10cla64__8ohm"]).record
    assert "sensitivity_db_1w1m" not in record["parameters"]
    assert record["manufacturer"] == "BC"


def test_compression_driver_without_ts_is_incomplete_and_others_are_skipped():
    data = entries()
    compression = convert_entry(data["faital__fd371__8ohm"])
    assert compression.incomplete and compression.record["driver_type"] == "compression"
    assert compression.record["parameters"]["sensitivity_db_1w1m"] == 107.0
    assert "missing Thiele/Small fields" in compression.reason
    passive = convert_entry(data["dayton__dma105_pr__dma105_pr_4_dma_series_passive_radiator"])
    assert passive.record is None and "passive radiator" in passive.reason
    jensen = convert_entry(data["jensen__10_blackbird_100__8ohm"])
    assert jensen.record is None and "canonical list" in jensen.reason


def test_fill_needs_matching_core_values_and_never_changes_existing_ones():
    new = convert_entry(entries()["eighteensound__10cx650__8ohm"]).record
    updated, filled = fill_missing(legacy_10cx650(), new)
    assert set(filled) == {"power_w", "sensitivity_db_1w1m", "usable_f_low_hz", "usable_f_high_hz"}
    assert updated["parameters"]["fs_hz"] == 54.0 and updated["field_provenance"]["power_w"]["source_id"] == "driver-base"
    assert updated["provenance"]["known_limitations"] == ["Rms not captured."]
    other = legacy_10cx650()
    other["parameters"]["bl_tm"] = 18.0
    assert fill_missing(other, new) == (None, [])
    assert [row["field"] for row in compare_records(other, new)] == ["bl_tm"]


def test_run_dedupes_fills_reports_conflicts_and_skips_ambiguous(tmp_path):
    db, curated, incomplete, reports = (tmp_path / name for name in ("db", "curated", "incomplete", "reports"))
    (db / "18Sound").mkdir(parents=True)
    (db / "18Sound" / "18sound-10cx650.json").write_text(json.dumps(legacy_10cx650(), indent=4))
    # A hand-entered record for the 8 ohm 10CLA64 whose Fs disagrees: reported, not touched.
    curated_record = {"driver_id": "bc-10cla64", "manufacturer": "BC", "model_name": "10CLA64",
                      "parameters": {"fs_hz": 60.0, "re_ohm": 5.5, "bl_tm": 14.5, "sd_m2": 0.032, "mms_kg": 0.036,
                                     "le_h": 0.0009, "nominal_impedance_ohm": 8.0}}
    (curated / "BC").mkdir(parents=True)
    (curated / "BC" / "bc-10cla64.json").write_text(json.dumps(curated_record))

    fetcher = FixtureFetcher()
    summary = DriverBaseSource(fetcher).run(db, incomplete, existing_dirs=[db, curated], report_dir=reports)
    result = summary.to_dict()

    assert set(fetcher.calls) == {DATASET_URL}
    assert result["matched_existing"] == 2 and result["filled_existing"] == 1
    assert result["written"] == 1 and result["incomplete"] == 1 and result["ambiguous_in_dataset"] == 2
    filled = json.loads((db / "18Sound" / "18sound-10cx650.json").read_text())
    assert filled["parameters"]["power_w"] == 200.0 and filled["record_schema_version"] == 1
    assert json.loads((curated / "BC" / "bc-10cla64.json").read_text()) == curated_record
    new = json.loads((db / "BC" / "bc-10cla64-16ohm.json").read_text())
    assert new["identity"]["variant"]["impedance_ohm"] == 16.0
    assert (incomplete / "Faital" / "faital-fd371.json").exists()
    assert not (db / "Monacor").exists()

    rows = list(csv.DictReader((reports / "conflicts.csv").open()))
    assert [(r["driver_id"], r["field"]) for r in rows] == [("bc-10cla64", "fs_hz")]
    report = json.loads((reports / "summary.json").read_text())
    assert report["license"] == "CC BY 4.0" and "driver-base" in report["attribution"]


def test_inconsistent_mms_sends_the_record_to_incomplete():
    # Beyma 18QLEX1600Fe 4 ohm as driver-base parsed it: Mms 36 g where Qes, Bl, Fs and Re imply about 336 g.
    entry = {**entries()["eighteensound__10cx650__8ohm"], "fs_hz": 31.0, "re_ohm": 3.6, "bl_tm": 30.7,
             "mms_g": 36.0, "cms_mm_per_n": 0.08, "qes": 0.25, "qms": 11.6, "qts": 0.25, "sd_cm2": 1255.0}
    entry["spec_source"] = {**entry["spec_source"], "cms_mm_per_n": "html_table"}
    result = convert_entry(entry)
    assert result.incomplete and result.reason.startswith("inconsistent Thiele/Small values")
    assert "Qes 0.25" in result.reason


def test_cms_that_disagrees_with_fs_is_dropped_when_qes_agrees():
    # Celestion TSQ1545: Cms 0.88 mm/N where Fs and Mms imply about 0.089 mm/N.
    entry = {**entries()["eighteensound__10cx650__8ohm"], "fs_hz": 39.0, "re_ohm": 5.05, "bl_tm": 25.0,
             "mms_g": 187.8, "cms_mm_per_n": 0.88, "qes": 0.372, "qms": 5.612, "qts": 0.349, "sd_cm2": 855.3}
    entry["spec_source"] = {**entry["spec_source"], "cms_mm_per_n": "html_table"}
    result = convert_entry(entry)
    assert result.reason is None and not result.incomplete
    assert "cms_m_per_n" not in result.record["parameters"]
    assert "cms_m_per_n" not in result.record["field_provenance"]
    assert any(note.startswith("Cms not imported") for note in result.record["provenance"]["notes"])


def test_warn_flags_are_kept_and_non_aes_power_is_marked():
    entry = {**entries()["eighteensound__10cx650__8ohm"], "warn_flags": ["power_std_non_aes", "range_reject:x=1"]}
    record = convert_entry(entry).record
    assert "not AES" in record["field_provenance"]["power_w"]["power_standard"]
    assert any("power_std_non_aes, range_reject:x=1" in note for note in record["provenance"]["notes"])


def test_page_nominal_size_is_preferred_over_the_sd_guess():
    # SB Acoustics MW19TX: a 7 inch driver whose Sd alone suggests 6 inch.
    entry = {**entries()["eighteensound__10cx650__8ohm"], "nominal_size_mm": 180.0, "sd_cm2": 150.0}
    entry["spec_source"] = {**entry["spec_source"], "nominal_size_mm": "html_table"}
    assert convert_entry(entry).record["nominal_diameter"] == "7in"
    entry["spec_source"] = {**entry["spec_source"], "nominal_size_mm": "derived"}
    assert convert_entry(entry).record["nominal_diameter"] == "6in"


def test_only_impedance_suffixes_match_by_values():
    new = convert_entry(entries()["eighteensound__10cx650__8ohm"]).record
    coil = {**legacy_10cx650(), "model_name": "10CX650-8"}
    revision = {**legacy_10cx650(), "model_name": "10CX650EVO"}
    numbered = {**legacy_10cx650(), "model_name": "10CX650-2"}  # "-2" on an 8 ohm driver is a revision
    assert same_driver_by_values(coil, new) and prefix_candidate(coil, new) is None
    for other, suffix in ((revision, "EVO"), (numbered, "2")):
        assert not same_driver_by_values(other, new)
        assert prefix_candidate(other, new) == suffix


def test_fill_never_removes_an_existing_peak_rating():
    new = convert_entry(entries()["eighteensound__10cx650__8ohm"]).record
    old = legacy_10cx650()
    old["parameters"]["peak_power_w"] = 150.0  # below the 200 W that will be filled
    updated, filled = fill_missing(old, new)
    assert "power_w" in filled and updated["parameters"]["peak_power_w"] == 150.0
    new["parameters"]["peak_power_w"] = 100.0
    new["field_provenance"]["peak_power_w"] = new["field_provenance"]["power_w"]
    updated, filled = fill_missing(legacy_10cx650(), new)
    assert "peak_power_w" not in filled and "peak_power_w" not in updated["parameters"]


def test_sd_that_disagrees_with_the_page_vas_sends_the_record_to_incomplete():
    # RCF MB12N251 as driver-base parsed it: Sd 5000 cm2 on a 12 inch driver with Vas 83 L.
    entry = {**entries()["eighteensound__10cx650__8ohm"], "fs_hz": 50.0, "re_ohm": 5.1, "bl_tm": 15.0,
             "mms_g": 49.0, "qes": 0.36, "qms": 5.4, "qts": 0.34, "sd_cm2": 5000.0, "vas_liters": 83.0}
    entry["spec_source"] = {**entry["spec_source"], "vas_liters": "html_div_pairs"}
    result = convert_entry(entry)
    assert result.incomplete and "page Vas 83 L" in result.reason
    entry["sd_cm2"] = 500.0
    assert not convert_entry(entry).incomplete


def test_page_nominal_size_is_ignored_when_sd_contradicts_it():
    # Dayton CFC180 is listed with a 2 inch size (the voice coil) but Sd of a 6.5 inch cone.
    entry = {**entries()["eighteensound__10cx650__8ohm"], "nominal_size_mm": 50.8}
    entry["spec_source"] = {**entry["spec_source"], "nominal_size_mm": "html_table"}
    assert convert_entry(entry).record["nominal_diameter"] == "10in"
