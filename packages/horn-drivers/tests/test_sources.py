"""Source adapter interface, registries and the retailer matcher, with offline fakes."""

import csv
import json
from pathlib import Path

import pytest

from horn_drivers.fetch import FetchResult, SourceBlocked
from horn_drivers.prices import load_offers
from horn_drivers.records import make_record
from horn_drivers.sources import (ADAPTERS, DriverIndex, DriverResult, DriverSource, Listing, PriceSource, ProductRef,
                                  create_adapter, load_source_registry, register_adapter)
from horn_drivers.sources.registry import SourceEntry, save_source_registry

ROOT = Path(__file__).resolve().parents[3]
PARAMS = {"fs_hz": 61.0, "re_ohm": 5.1, "bl_tm": 12.3, "sd_m2": 0.022, "mms_kg": 0.022, "le_h": 0.0004,
          "qms": 6.0, "qes": 0.36, "nominal_impedance_ohm": 8.0}


class FakeFetcher:
    """Stands in for PoliteFetcher: serves canned pages, can block."""

    def __init__(self, source_id, pages, block_after=None):
        self.source_id, self.pages, self.block_after, self.calls = source_id, pages, block_after, []

    def get(self, url, max_age_days=None):
        if self.block_after is not None and len(self.calls) >= self.block_after:
            raise SourceBlocked("two 429 in a row")
        self.calls.append(url)
        status, text = self.pages.get(url, (404, ""))
        return FetchResult(url, status, text.encode(), "text/html", "2026-11-17T10:00:00+00:00", False)

    def status(self):
        return {"source_id": self.source_id}


class FakeMaker(DriverSource):
    source_id = "fake-maker"
    kind = "manufacturer"

    def discover(self, limit=None):
        refs = [ProductRef("https://maker.test/a", "BC", "8NDL51"), ProductRef("https://maker.test/bad", "BC", "Bad"),
                ProductRef("https://maker.test/missing", "BC", "Gone")]
        return iter(refs[:limit])

    def parse_driver(self, page, ref):
        values = json.loads(page.text)
        provenance = {k: {"source_id": self.source_id, "url": page.url, "retrieved_on": page.retrieved_on, "method": "page"}
                      for k in values}
        return DriverResult(ref, make_record(manufacturer=ref.manufacturer, model=ref.model, parameters=values,
                                             field_provenance=provenance, driver_type="cone"))


class FakeShop(PriceSource):
    source_id = "fake-shop"
    region, currency, vat_included = "UK", "GBP", True

    def discover(self, limit=None):
        return iter([ProductRef(f"https://shop.test/{n}") for n in ("exact", "fuzzy", "soft404", "other")])

    def parse_listing(self, page, ref):
        if "Unknown Item" in page.text:
            return None
        data = json.loads(page.text)
        return Listing(url=page.url, currency="GBP", **data)


def test_driver_source_writes_valid_records_and_routes_rejects(tmp_path):
    pages = {"https://maker.test/a": (200, json.dumps(PARAMS)),
             "https://maker.test/bad": (200, json.dumps({**PARAMS, "qts": 0.9}))}
    adapter = FakeMaker(FakeFetcher("fake-maker", pages))
    summary = adapter.run(tmp_path / "db", tmp_path / "incomplete")
    assert (summary.written, summary.incomplete, summary.pages) == (1, 1, 3)
    assert summary.skipped == [("https://maker.test/missing", "HTTP 404")]
    written = json.loads((tmp_path / "db" / "BC" / "bc-8ndl51.json").read_text())
    assert written["record_schema_version"] == 2
    rejected = json.loads((tmp_path / "incomplete" / "BC" / "bc-bad.json").read_text())
    assert "loader rejects record" in rejected["incomplete_reason"]


def test_blocked_source_stops_the_run_and_reports_why(tmp_path):
    adapter = FakeMaker(FakeFetcher("fake-maker", {"https://maker.test/a": (200, json.dumps(PARAMS))}, block_after=1))
    summary = adapter.run(tmp_path / "db", tmp_path / "incomplete")
    assert summary.written == 1 and "two 429" in summary.blocked_reason


def test_adapter_needs_fetcher_for_its_own_source():
    with pytest.raises(ValueError):
        FakeMaker(FakeFetcher("someone-else", {}))


def records():
    eight = make_record(manufacturer="Eminence", model="Beta-12A", parameters=PARAMS,
                        field_provenance={k: {"source_id": "x", "url": "https://x.test", "retrieved_on": "2026-11-17",
                                              "method": "page"} for k in PARAMS}, aliases=["Beta 12 A"])
    single = make_record(manufacturer="BC", model="8NDL51", parameters={k: v for k, v in PARAMS.items() if k != "nominal_impedance_ohm"},
                         field_provenance={k: {"source_id": "x", "url": "https://x.test", "retrieved_on": "2026-11-17",
                                               "method": "page"} for k in PARAMS if k != "nominal_impedance_ohm"})
    return [eight, single]


def test_index_matches_exactly_or_sends_to_review():
    index = DriverIndex(records())
    listing = lambda **kw: Listing(**{"url": "https://s.test", "title": "Eminence Beta-12A 8 Ohm", "manufacturer": "Eminence",
                                      "model": "Beta-12A", "price": 99.0, "currency": "GBP", "impedance_ohm": 8, **kw})
    assert index.match(listing())[:2] == ("eminence-beta-12a", "exact")
    assert index.match(listing(model="Beta 12 A", title="Eminence Beta 12 A"))[:2] == ("eminence-beta-12a", "exact")
    assert index.match(listing(impedance_ohm=16))[:2] == (None, None)
    assert index.match(listing(title="Speaker Beta-12A"))[2] == "title does not name the manufacturer"
    assert index.match(listing(title="Eminence Kappa"))[2] == "title does not name the model"
    assert index.match(listing(manufacturer="BC", model="8NDL51", title="B&C 8NDL51", impedance_ohm=None))[:2] == ("bc-8ndl51", "fuzzy")


def test_price_source_writes_exact_offers_and_review_csv(tmp_path):
    pages = {
        "https://shop.test/exact": (200, json.dumps({"title": "Eminence Beta-12A 8 Ohm", "manufacturer": "Eminence",
                                                       "model": "Beta-12A", "impedance_ohm": 8, "price": 99.0, "sku": "S1", "in_stock": True})),
        "https://shop.test/fuzzy": (200, json.dumps({"title": "B&C 8NDL51", "manufacturer": "B&C", "model": "8NDL51", "price": 150.0})),
        "https://shop.test/soft404": (200, "Unknown Item"),
        "https://shop.test/other": (200, json.dumps({"title": "Eminence Kappa Pro 15", "manufacturer": "Eminence",
                                                       "model": "Kappa Pro 15", "impedance_ohm": 8, "price": 200.0})),
    }
    shop = FakeShop(FakeFetcher("fake-shop", pages))
    review = tmp_path / "review.csv"
    summary = shop.run(DriverIndex(records()), tmp_path / "prices", review)
    assert (summary.written, summary.review) == (1, 2)
    offers = {o.driver_id: o for o in load_offers(tmp_path / "prices")}
    assert offers["eminence-beta-12a"].match_confidence == "exact"
    assert offers["eminence-beta-12a"].vat_included is True and offers["eminence-beta-12a"].region == "UK"
    assert "bc-8ndl51" not in offers           # below exact: review CSV only, never the price file
    rows = list(csv.DictReader(review.open()))
    assert {r["note"] for r in rows} == {"impedance not confirmed on both sides", "model not in the database"}
    assert (tmp_path / "prices" / "fake-shop" / "2026-11-17.jsonl").exists()


def test_registry_file_is_valid_and_rates_meet_the_floors():
    entries = load_source_registry(ROOT / "data" / "sources.json")
    assert entries["loudspeakerdatabase"].min_delay_s >= 30
    assert entries["celestion"].min_delay_s >= 10
    assert all(e.min_delay_s >= 5 for e in entries.values())
    assert all(e.region and e.currency for e in entries.values() if e.kind == "retailer")
    assert entries["loudspeakerfreaks"].approval_status == "excluded"


def test_create_adapter_refuses_unapproved_sources(tmp_path):
    register_adapter(FakeShop)
    try:
        path = tmp_path / "sources.json"
        entry = SourceEntry(id="fake-shop", name="Fake", kind="retailer", base_url="https://shop.test",
                            approval_status="pending", min_delay_s=5, cache_days=7, region="UK", currency="GBP")
        save_source_registry({"fake-shop": entry}, path)
        with pytest.raises(PermissionError):
            create_adapter("fake-shop", path, tmp_path / "cache")
        entry.approval_status = "approved"
        save_source_registry({"fake-shop": entry}, path)
        adapter = create_adapter("fake-shop", path, tmp_path / "cache")
        assert adapter.fetcher.min_delay_s == 5 and adapter.fetcher.cache_root == tmp_path / "cache"
        entry.min_delay_s = 1
        save_source_registry({"fake-shop": entry}, path)
        with pytest.raises(ValueError, match="floor"):
            load_source_registry(path)
    finally:
        ADAPTERS.pop("fake-shop", None)
