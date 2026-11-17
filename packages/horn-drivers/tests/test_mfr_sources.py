"""Manufacturer adapters (B&C, 18 Sound, Celestion, Beyma) and the merge step, offline.

The fixture pages under fixtures/sources are cut down from real responses to
the values the parsers read (no product prose or images).
"""

import csv
import json
from pathlib import Path

import pytest

from horn_drivers.fetch import FetchResult
from horn_drivers.loader import _driver_from_dict
from horn_drivers.merge import merge_staged
from horn_drivers.records import check_structure, write_record
from horn_drivers.sources import ProductRef
from horn_drivers.sources.bc import BCSource
from horn_drivers.sources.beyma import BeymaSource, model_from_title
from horn_drivers.sources.celestion import CelestionSource
from horn_drivers.sources.common import (sensitivity_1w, split_value_unit, stated_drive, states_one_watt, to_si)
from horn_drivers.sources.eighteensound import EighteenSoundSource

FIXTURES = Path(__file__).parent / "fixtures" / "sources"


class CannedFetcher:
    def __init__(self, source_id, pages=None):
        self.source_id, self.pages = source_id, pages or {}

    def get(self, url, max_age_days=None):
        status, text = self.pages.get(url, (404, ""))
        return FetchResult(url, status, text.encode(), "text/html", "2026-11-17T10:00:00+00:00", False)

    def status(self):
        return {}


def page(name, url):
    return FetchResult(url, 200, (FIXTURES / name).read_bytes(), "text/html", "2026-11-17T10:00:00+00:00", False)


def adapter(cls, pages=None):
    return cls(CannedFetcher(cls.source_id, pages))


def assert_loadable(record):
    errors, _ = check_structure(record)
    assert errors == []
    _driver_from_dict(record)


# -- helpers -----------------------------------------------------------------

@pytest.mark.parametrize("text, expected", [
    ("346.0 cm2 (53.63 in2)", (346.0, "cm2")),
    ("530.93cm2 / 82.29in2", (530.93, "cm2")),
    ("44.51g/1.57oz", (44.51, "g")),
    ("0.18mm/N", (0.18, "mm/N")),
    ("0.37 mH", (0.37, "mH")),
    ("14.5 Txm", (14.5, "Txm")),
    ("163 µm/N", (163.0, "µm/N")),
    ("5.0", (5.0, "")),
])
def test_split_value_unit(text, expected):
    assert split_value_unit(text) == expected


def test_to_si_and_unknown_unit():
    assert to_si(346.0, "cm2") == pytest.approx(0.0346)
    assert to_si(29.5, "g") == pytest.approx(0.0295)
    assert to_si(0.37, "mH") == pytest.approx(0.00037)
    assert to_si(163, "µm/N") == pytest.approx(163e-6)
    with pytest.raises(KeyError):
        to_si(1.0, "furlong")


def test_sensitivity_drive_levels():
    assert stated_drive("Applied RMS Voltage is set to 4V for 16 ohm Nominal Impedance") == (4.0, 16.0)
    assert stated_drive("Applied RMS Voltage is set to 2.83 V for 8 ohms Nominal Impedance.") == (2.83, 8.0)
    assert stated_drive("Measured on axis at 1W, 1m") is None
    assert sensitivity_1w(94.0, 4.0, 16.0) == 94.0
    assert sensitivity_1w(94.0, 2.83, 8.0) == 94.0
    assert sensitivity_1w(97.0, 2.83, 16.0) == 100.0     # 2.83 V into 16 ohm is 0.5 W
    assert states_one_watt("100 dB 1 W @ 1 m @ ZN")
    assert states_one_watt("Measured on axis at 1W, 1m in 2pi anechoic environment.")
    assert not states_one_watt("100 dB 2.83 V @ 1 m")


# -- B&C ---------------------------------------------------------------------

def test_bc_listing_keeps_pilot_sizes_and_every_variant():
    refs = list(BCSource.refs_from_listing((FIXTURES / "bc-listing.html").read_text()))
    assert [(r.model, r.hint("inches"), r.hint("ohms")) for r in refs] == [
        ("6MD38", "6.5", "8"), ("6MD38", "6.5", "16"), ("6PEV13", "6.5", "8"), ("8PE21", "8", "8"),
        ("8PS21", "8", "8")]
    assert refs[1].url == "https://www.bcspeakers.com/en/products/lf-driver/6.5/16/6MD38"
    assert len(list(BCSource.refs_from_listing((FIXTURES / "bc-listing.html").read_text(), limit=2))) == 2


def test_bc_product_page_gives_si_values_power_band_and_1w_sensitivity():
    url = "https://www.bcspeakers.com/en/products/lf-driver/8/16/8NDL51"
    result = adapter(BCSource).parse_driver(page("bc-8ndl51-16.html", url), ProductRef(url))
    record = result.record
    assert result.reason is None
    assert record["driver_id"] == "bc-8ndl51-16ohm"
    assert record["manufacturer"] == "BC" and record["model_name"] == "8NDL51"
    params = record["parameters"]
    assert params["fs_hz"] == 70.0 and params["re_ohm"] == 10.6 and params["bl_tm"] == 15.2
    assert params["sd_m2"] == pytest.approx(0.022) and params["mms_kg"] == pytest.approx(0.026)
    assert params["le_h"] == pytest.approx(0.0006) and params["xmax_m"] == pytest.approx(0.007)
    assert params["nominal_impedance_ohm"] == 16.0
    assert (params["power_w"], params["peak_power_w"]) == (200.0, 400.0)
    assert params["sensitivity_db_1w1m"] == 94.0        # stated at 4 V into 16 ohm, which is 1 W
    assert record["field_provenance"]["sensitivity_db_1w1m"]["method"] == "derived"
    assert record["field_provenance"]["fs_hz"] == {"source_id": "bc", "url": url, "retrieved_on": "2026-11-17",
                                                   "method": "page"}
    assert (record["usable_f_low_hz"], record["usable_f_high_hz"]) == (70.0, 3000.0)
    assert record["nominal_diameter"] == "8in"
    assert_loadable(record)


def test_bc_sensitivity_with_mismatched_drive_note_is_not_stored():
    html = (FIXTURES / "bc-8ndl51-16.html").read_text().replace("4V for 16 ohm", "2.83 V for 8 ohms")
    url = "https://www.bcspeakers.com/en/products/lf-driver/8/16/8NDL51"
    fetched = FetchResult(url, 200, html.encode(), "text/html", "2026-11-17T10:00:00+00:00", False)
    record = adapter(BCSource).parse_driver(fetched, ProductRef(url)).record
    assert "sensitivity_db_1w1m" not in record["parameters"]
    assert any("1 W figure is unknown" in note for note in record["provenance"]["notes"])


def test_bc_page_without_product_data_is_skipped():
    url = "https://www.bcspeakers.com/en/products/lf-driver/8/8/X"
    fetched = FetchResult(url, 200, b"<html>maintenance</html>", "text/html", "2026-11-17T10:00:00+00:00", False)
    result = adapter(BCSource).parse_driver(fetched, ProductRef(url))
    assert result.record is None and "no product data" in result.reason


# -- 18 Sound ----------------------------------------------------------------

def test_eighteensound_listing_dedupes_exact_paths_and_filters_sizes():
    refs = list(EighteenSoundSource.refs_from_listing((FIXTURES / "eighteensound-listing.html").read_text()))
    # 10NTLW2000 and 10ntlw2000 are different pages on the site (10NTLS2000 and
    # 10NTLW2000), so paths that differ only in case are both kept.
    assert [r.url.rsplit("/lf-driver/", 1)[1] for r in refs] == [
        "15-0/8/15LW1401", "15-0/4/15lw1401%20-%204ohm", "10-0/8/10nmb420", "10-0/8/10NTLW2000",
        "10-0/8/10ntlw2000", "10-0/16/10NMB420", "6-5/8/6NMB420", "5-0/8/5W430"]
    assert refs[6].hint("inches") == "6.5"


def test_eighteensound_product_page():
    url = "https://www.eighteensound.it/en/products/lf-driver/10-0/8/10nmb420"
    ref = ProductRef(url, hints=(("inches", "10"),))
    record = adapter(EighteenSoundSource).parse_driver(page("eighteensound-10nmb420-8.html", url), ref).record
    params = record["parameters"]
    assert record["driver_id"] == "18sound-10nmb420-8ohm" and record["manufacturer"] == "18Sound"
    assert params["fs_hz"] == 56.0 and params["bl_tm"] == 14.5 and params["sd_m2"] == pytest.approx(0.0346)
    assert params["mms_kg"] == pytest.approx(0.0295) and params["le_h"] == pytest.approx(0.00037)
    assert params["xmax_m"] == pytest.approx(0.004) and params["qts"] == 0.25
    assert (params["power_w"], params["peak_power_w"], params["sensitivity_db_1w1m"]) == (350.0, 500.0, 99.0)
    assert (record["usable_f_low_hz"], record["usable_f_high_hz"]) == (65.0, 5000.0)
    assert_loadable(record)


# -- Celestion ---------------------------------------------------------------

def test_celestion_product_page():
    url = "https://celestion.com/product/tf1220/"
    record = adapter(CelestionSource).parse_driver(page("celestion-tf1220.html", url), ProductRef(url)).record
    params = record["parameters"]
    assert record["driver_id"] == "celestion-tf1220-8ohm"
    assert params["sd_m2"] == pytest.approx(0.053093) and params["mms_kg"] == pytest.approx(0.04451)
    assert params["bl_tm"] == 13.91 and params["le_h"] == pytest.approx(0.00071)
    assert params["cms_m_per_n"] == pytest.approx(0.00018) and params["rms_kg_per_s"] == 4.69
    assert (params["power_w"], params["peak_power_w"]) == (150.0, 300.0)     # AES and continuous, not EIA
    assert params["sensitivity_db_1w1m"] == 97.0
    assert record["field_provenance"]["sensitivity_db_1w1m"]["method"] == "page"
    assert record["nominal_diameter"] == "12in"
    assert_loadable(record)


def test_wrong_unit_and_inconsistent_cms_are_not_stored():
    # Two real page errors: Rms shown in mm/N (celestion-cf1025bmb) and a Cms
    # ten times too high (celestion-tsq1545).
    url = "https://celestion.com/product/tf1220/"
    html = (FIXTURES / "celestion-tf1220.html").read_text()
    html = html.replace("<div>4.69kg/s</div>", "<div>4.69mm/N</div>").replace("<div>0.18mm/N</div>",
                                                                            "<div>1.8mm/N</div>")
    fetched = FetchResult(url, 200, html.encode(), "text/html", "2026-11-17T10:00:00+00:00", False)
    record = adapter(CelestionSource).parse_driver(fetched, ProductRef(url)).record
    params, notes = record["parameters"], record["provenance"]["notes"]
    assert "rms_kg_per_s" not in params and "rms_kg_per_s" not in record["field_provenance"]
    assert "cms_m_per_n" not in params and "cms_m_per_n" not in record["field_provenance"]
    assert any("rms_kg_per_s shown as 4.69 mm/N" in note for note in notes)
    assert any("page Cms 0.0018 m/N gives Fs 17.8 Hz" in note for note in notes)
    assert params["fs_hz"] == 55.8
    assert_loadable(record)


def test_inconsistent_rms_is_not_stored():
    url = "https://celestion.com/product/tf1220/"
    html = (FIXTURES / "celestion-tf1220.html").read_text().replace("<div>4.69kg/s</div>", "<div>6.1kg/s</div>")
    fetched = FetchResult(url, 200, html.encode(), "text/html", "2026-11-17T10:00:00+00:00", False)
    record = adapter(CelestionSource).parse_driver(fetched, ProductRef(url)).record
    assert "rms_kg_per_s" not in record["parameters"]
    assert record["parameters"]["cms_m_per_n"] == pytest.approx(0.00018)
    assert any("page Rms 6.1 kg/s" in note for note in record["provenance"]["notes"])


def test_run_counts_parser_errors(tmp_path):
    class Broken(BCSource):
        def parse_driver(self, page, ref):
            raise KeyError("bug")
    listing_url = "https://www.bcspeakers.com/en/products/lf-driver"
    pages = {listing_url: (200, (FIXTURES / "bc-listing.html").read_text())}
    pages.update({ref.url: (200, "x") for ref in BCSource.refs_from_listing(pages[listing_url][1])})
    summary = adapter(Broken, pages).run(tmp_path / "db", tmp_path / "incomplete")
    assert summary.parse_errors == summary.pages == 5
    assert summary.too_many_parse_errors() and summary.to_dict()["parse_errors"] == 5


def test_celestion_outside_size_range_is_skipped():
    url = "https://celestion.com/product/tf1220/"
    html = (FIXTURES / "celestion-tf1220.html").read_text().replace("305mm / 12in", "460mm / 18in")
    fetched = FetchResult(url, 200, html.encode(), "text/html", "2026-11-17T10:00:00+00:00", False)
    result = adapter(CelestionSource).parse_driver(fetched, ProductRef(url))
    assert result.record is None and "outside" in result.reason


# -- Beyma -------------------------------------------------------------------

def test_beyma_title_and_listing():
    assert model_from_title("LOUDSPEAKER 12MI100 8 Ohm") == ("12MI100", 8.0)
    assert model_from_title("LOUDSPEAKER 10BR60/V2 8 Ohm") == ("10BR60/V2", 8.0)
    assert model_from_title("LOUDSPEAKER 5P200Nd/N 8") == ("5P200Nd/N", 8.0)
    refs = list(BeymaSource.refs_from_listing((FIXTURES / "beyma-listing.html").read_text()))
    assert len(refs) == len({r.sku for r in refs}) > 3


def test_beyma_product_page():
    url = "https://www.beyma.com/en/products/c/low-mid-frequency/112MI108/altavoz-12mi100-8-oh/"
    record = adapter(BeymaSource).parse_driver(page("beyma-12mi100.html", url), ProductRef(url, sku="112MI108")).record
    params = record["parameters"]
    assert record["driver_id"] == "beyma-12mi100-8ohm" and record["identity"]["aliases"] == ["112MI108"]
    assert params["bl_tm"] == 21.6 and params["mms_kg"] == 0.047 and params["sd_m2"] == 0.053
    assert params["cms_m_per_n"] == pytest.approx(163e-6) and params["le_h"] == pytest.approx(0.001)
    assert (params["power_w"], params["peak_power_w"], params["sensitivity_db_1w1m"]) == (450.0, 900.0, 100.0)
    assert_loadable(record)


# -- run() and merge ---------------------------------------------------------

def test_run_writes_loadable_records(tmp_path):
    listing_url = "https://www.bcspeakers.com/en/products/lf-driver"
    product_url = "https://www.bcspeakers.com/en/products/lf-driver/6.5/8/6MD38"
    pages = {listing_url: (200, (FIXTURES / "bc-listing.html").read_text()),
             product_url: (200, (FIXTURES / "bc-8ndl51-16.html").read_text())}
    summary = adapter(BCSource, pages).run(tmp_path / "db", tmp_path / "incomplete")
    assert summary.pages == 5 and summary.written == 1 and len(summary.skipped) == 4    # four 404s
    assert (tmp_path / "db" / "BC" / "bc-8ndl51-16ohm.json").exists()


def legacy(driver_id, model, **params):
    base = {"fs_hz": 66.0, "re_ohm": 5.3, "bl_tm": 12.4, "sd_m2": 0.022, "mms_kg": 0.028, "le_h": 0.0005,
            "qms": 3.6, "qes": 0.41, "qts": 0.37, "nominal_impedance_ohm": 8.0, "xmax_m": 0.007}
    base.update(params)
    return {"driver_id": driver_id, "manufacturer": "BC", "model_name": model, "driver_type": "cone",
            "parameters": base, "parameter_source": f"https://loudspeakerdatabase.com/BC/{model}",
            "record_schema_version": 1,
            "provenance": {"retrieved_on": "2026-02-26",
                           "known_limitations": ["mms_kg may be the dry moving mass Mmd where the page showed no Mms."]}}


def staged_bc(html_name="bc-8ndl51-16.html", **replace):
    url = "https://www.bcspeakers.com/en/products/lf-driver/8/16/8NDL51"
    html = (FIXTURES / html_name).read_text()
    for old, new in replace.items():
        html = html.replace(old, new)
    fetched = FetchResult(url, 200, html.encode(), "text/html", "2026-11-17T10:00:00+00:00", False)
    return adapter(BCSource).parse_driver(fetched, ProductRef(url)).record


def test_merge_keeps_id_records_conflicts_and_writes_report(tmp_path):
    db = tmp_path / "db"
    # Old record for the 16 ohm variant: Fs 3 % off (kept quiet), Bl 20 % off, Le 8 % off (under 10 %),
    # rms only in the old record (carried over).
    write_record(db, legacy("bc-8ndl51-16", "8NDL51", nominal_impedance_ohm=16.0, fs_hz=72.0, bl_tm=12.2,
                            le_h=0.00065, re_ohm=10.6, rms_kg_per_s=1.5))
    report = tmp_path / "conflicts.csv"
    summary = merge_staged([staged_bc()], db, report)
    assert summary.updated == ["bc-8ndl51-16"] and summary.added == []
    merged = json.loads((db / "BC" / "bc-8ndl51-16.json").read_text())
    assert merged["record_schema_version"] == 2
    assert merged["parameters"]["bl_tm"] == 15.2                  # manufacturer wins
    conflicts = merged["field_provenance"]["bl_tm"]["conflicts"]
    assert conflicts == [{"source_id": "loudspeakerdatabase", "url": "https://loudspeakerdatabase.com/BC/8NDL51",
                          "value": 12.2, "retrieved_on": "2026-02-26"}]
    assert "conflicts" not in merged["field_provenance"]["fs_hz"]
    assert "conflicts" not in merged["field_provenance"]["le_h"]
    assert merged["parameters"]["rms_kg_per_s"] == 1.5
    assert merged["field_provenance"]["rms_kg_per_s"]["source_id"] == "loudspeakerdatabase"
    assert "Mmd" in merged["field_provenance"]["mms_kg"]["note"]
    assert_loadable(merged)
    rows = list(csv.DictReader(report.open()))
    assert {(r["field"], r["previous_value"]) for r in rows} >= {("bl_tm", "12.2"), ("mms_kg", "0.028")}
    assert all(r["driver_id"] == "bc-8ndl51-16" for r in rows)

    # A rerun of the same adapter keeps the earlier conflicts and adds none.
    again = merge_staged([staged_bc()], db, tmp_path / "again.csv")
    assert again.conflicts == [] and again.updated == ["bc-8ndl51-16"]
    rerun = json.loads((db / "BC" / "bc-8ndl51-16.json").read_text())
    assert rerun["field_provenance"]["bl_tm"]["conflicts"] == conflicts
    # ... and keeps the field carried over from the older source, with its provenance and note.
    assert rerun["parameters"]["rms_kg_per_s"] == 1.5
    assert rerun["field_provenance"]["rms_kg_per_s"] == merged["field_provenance"]["rms_kg_per_s"]
    assert any(note.startswith("Replaces the loudspeakerdatabase record")
               for note in rerun["provenance"]["notes"])
    assert_loadable(rerun)


def test_merge_ids_for_new_models(tmp_path):
    db = tmp_path / "db"
    single = staged_bc()
    eight = staged_bc(**{'"nominalImpedance": 16': '"nominalImpedance": 8'})
    # One variant of a model not in the db: plain id.
    summary = merge_staged([single], db)
    assert summary.added == ["bc-8ndl51"]
    # Two variants of a new model: each keeps its impedance suffix.
    db2 = tmp_path / "db2"
    other = json.loads(json.dumps(eight))
    other["driver_id"] = "bc-8ndl51-8ohm"
    other["identity"]["variant"]["impedance_ohm"] = 8.0
    other["parameters"]["nominal_impedance_ohm"] = 8.0
    summary = merge_staged([single, other], db2)
    assert sorted(summary.added) == ["bc-8ndl51-16ohm", "bc-8ndl51-8ohm"]


def test_merge_leaves_unmatched_old_records(tmp_path):
    db = tmp_path / "db"
    write_record(db, legacy("bc-12ps100", "12PS100"))
    summary = merge_staged([staged_bc()], db)
    assert summary.untouched_previous == 1
    assert json.loads((db / "BC" / "bc-12ps100.json").read_text())["record_schema_version"] == 1


def test_stored_conflict_is_a_structure_warning(tmp_path):
    write_record(tmp_path, legacy("bc-8ndl51-16", "8NDL51", nominal_impedance_ohm=16.0, bl_tm=12.2, re_ohm=10.6))
    merge_staged([staged_bc()], tmp_path)
    merged = json.loads((tmp_path / "BC" / "bc-8ndl51-16.json").read_text())
    errors, warnings = check_structure(merged)
    assert errors == []
    assert any(w.startswith("bl_tm: bc value 15.2 kept; loudspeakerdatabase gave 12.2") for w in warnings)


def test_scraper_does_not_overwrite_manufacturer_records(tmp_path):
    from horn_drivers.scraper import _save_driver
    record = staged_bc()
    write_record(tmp_path, record)
    scraped = legacy(record["driver_id"], "8NDL51")
    assert _save_driver(tmp_path, scraped) is False
    kept = json.loads((tmp_path / "BC" / f"{record['driver_id']}.json").read_text())
    assert kept["record_schema_version"] == 2
    assert _save_driver(tmp_path, legacy("bc-other", "OTHER")) is True
