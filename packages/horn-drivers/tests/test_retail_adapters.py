"""Retailer adapters: parsers against recorded fixture pages, matching, robots rules. No network."""

import csv
import json
from pathlib import Path

import pytest

from horn_drivers.fetch import FetchResult, robots_rules, rules_allow
from horn_drivers.prices import load_offers
from horn_drivers.records import make_record
from horn_drivers.sources import DriverIndex, Listing, ProductRef
from horn_drivers.sources.bax_shop import BaxShopUK
from horn_drivers.sources.dayton_msrp import DaytonMSRP
from horn_drivers.sources.lautsprechershop import Lautsprechershop
from horn_drivers.sources.madisound import Madisound
from horn_drivers.sources.retail import (brand_and_rest, impedance_from_text, parse_money, sitemap_locs,
                                         split_model_impedance)
from horn_drivers.sources.solen import Solen
from horn_drivers.sources.thomann import ThomannDE, ThomannUK

FIXTURES = Path(__file__).parent / "fixtures" / "retail"
PARAMS = {"fs_hz": 61.0, "re_ohm": 5.1, "bl_tm": 12.3, "sd_m2": 0.022, "mms_kg": 0.022, "le_h": 0.0004,
          "qms": 6.0, "qes": 0.36}


class FakeFetcher:
    def __init__(self, source_id, pages=None):
        self.source_id, self.pages, self.calls = source_id, pages or {}, []

    def get(self, url, max_age_days=None, accept=None):
        self.calls.append((url, accept))
        status, body = self.pages.get(url, (404, b""))
        return page(url, body, status)

    def status(self):
        return {"source_id": self.source_id}


def page(url, body, status=200):
    if isinstance(body, str):
        body = body.encode()
    return FetchResult(url, status, body, "text/html", "2026-11-17T10:00:00+00:00", True)


def fixture(name):
    return (FIXTURES / name).read_bytes()


def record(manufacturer, model, impedance, **kw):
    provenance = {k: {"source_id": "x", "url": "https://x.test", "retrieved_on": "2026-11-17", "method": "page"}
                  for k in PARAMS}
    return make_record(manufacturer=manufacturer, model=model, parameters={**PARAMS, "nominal_impedance_ohm": impedance},
                       field_provenance=provenance, impedance_ohm=impedance, **kw)


def index():
    return DriverIndex([record("BC", "12NDL88", 8.0), record("BC", "12PE32", 8.0), record("18Sound", "6NMB420", 8.0),
                        record("Celestion", "FTX1225", 8.0), record("Beyma", "12P80Fe", 8.0),
                        record("Dayton", "PA310-8", 8.0), record("Eminence", "Beta-12A", 8.0),
                        record("Eminence", "Kappa Pro-15A", 8.0)])


# -- helpers -------------------------------------------------------------------

def test_impedance_prefers_the_labelled_value_over_re():
    assert impedance_from_text("impedance R = 8 Ohm, DC resistance Re = 5,3 Ohm") == 8.0
    assert impedance_from_text("Re = 5,3 Ohm") is None
    assert impedance_from_text("Power: 1000W RMS @ 8 Ohm", allow_unlabelled=True) == 8.0
    assert impedance_from_text("dual 4 ohm and 8 ohm", allow_unlabelled=True) is None
    assert impedance_from_text("Impedanz: 4 Ohm") == 4.0


def test_model_suffix_money_and_brand_helpers():
    assert split_model_impedance("12NW76 8 Ohms") == ("12NW76", 8.0)
    assert split_model_impedance("Delta-12B - 16 Ohm") == ("Delta-12B", 16.0)
    assert split_model_impedance("Eighteen Sound 6 NMB 420, 16 Ohm") == ("Eighteen Sound 6 NMB 420", 16.0)
    assert split_model_impedance("12PE32") == ("12PE32", None)
    assert split_model_impedance("BL10-100X (8 Ohm)") == ("BL10-100X", 8.0)
    assert parse_money("205,30") == 205.30 and parse_money("1.234,50") == 1234.5
    assert parse_money("1,234.50") == 1234.5 and parse_money("CA$760.79") == 760.79 and parse_money("0") is None
    assert brand_and_rest("Eighteen Sound 6 NMB 420") == ("18Sound", "6 NMB 420")
    assert brand_and_rest("B&C 12PE32") == ("BC", "12PE32")
    assert brand_and_rest("Acme 12X") == (None, "Acme 12X")


def test_sitemap_reader_handles_gzip():
    import gzip
    xml = b'<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9"><url><loc>https://s.test/a</loc></url></urlset>'
    assert sitemap_locs(gzip.compress(xml)) == ["https://s.test/a"] == sitemap_locs(xml)


def test_rfc9309_rules_apply_wildcards_and_longest_match():
    text = ("User-agent: *\nDisallow: *?filter=\nDisallow: /intl/\nAllow: /intl/en$\nDisallow: /*sort_by*\n\n"
            "User-agent: GPTBot\nDisallow: /\n")
    rules = robots_rules(text, "horn-simulation-research")
    assert rules_allow(rules, "https://s.test/b-c_loudspeaker_12_inch.html")
    assert not rules_allow(rules, "https://s.test/cat.html?filter=x")
    assert not rules_allow(rules, "https://s.test/intl/x")
    assert rules_allow(rules, "https://s.test/intl/en")
    assert not rules_allow(rules, "https://s.test/eminence/?sort_by=price")
    named = robots_rules("User-agent: *\nDisallow: /\n\nUser-agent: horn-simulation-research\nAllow: /\n",
                         "horn-simulation-research")
    assert rules_allow(named, "https://s.test/x")


def test_title_fallback_finds_multiword_models_but_not_revisions():
    idx = index()
    assert idx.model_from_title("Eminence", "Eminence Kappa Pro-15A 15 inch 8 Ohm") == "KAPPAPRO15A"
    assert idx.model_from_title("Eminence", "Eminence Beta-12A 8 Ohm") == "BETA12A"
    assert idx.model_from_title("Eminence", "Eminence Beta-12A-2 12 inch") is None
    assert idx.model_from_title("Eminence", "Eminence Beta-12A 12 inch") == "BETA12A"


def test_title_fallback_sends_a_possible_model_suffix_to_review():
    idx = DriverIndex([record("Eminence", "Legend BP1525", 8.0), record("Celestion", "Truvox 0615", 8.0)])
    suffixed = Listing(url="https://s.test/a", title="Eminence Legend BP 1525 A", manufacturer="Eminence",
                       model="Legend BP 1525 A", price=160.0, currency="GBP", impedance_ohm=8.0)
    driver_id, confidence, note = idx.match(suffixed)
    assert (driver_id, confidence) == (None, None) and "suffix" in note
    for title in ("Eminence Legend BP 1525 8 Ohm", "Eminence Legend BP 1525 Speaker", "Eminence Legend BP 1525"):
        plain = Listing(url="https://s.test/b", title=title, manufacturer="Eminence", model="Legend BP 1525 x",
                        price=160.0, currency="GBP", impedance_ohm=8.0)
        assert idx.match(plain)[1] == "exact", title
    described = Listing(url="https://s.test/c", title="Celestion Truvox 0615 Mid/bass driver, 8 Ohm",
                        manufacturer="Celestion", model="Truvox", price=44.0, currency="GBP", impedance_ohm=8.0)
    assert idx.match(described)[1] == "exact"


@pytest.mark.parametrize("title, url", [
    ("B&C 12NDL88 8 Ohms B-Stock", "https://s.test/bc_12ndl88_8_ohms_b_stock.htm"),
    ("B&C 12NDL88 8 Ohm", "https://s.test/bc_12ndl88_8_ohms_b_stock.htm"),
    ("B&C 12NDL88 8 Ohm B-Ware", "https://s.test/x"),
    ("B&C 12NDL88 8 Ohm (used)", "https://s.test/x"),
    ("B&C 12NDL88 8 Ohm demo", "https://s.test/x"),
    ("B&C 12NDL88 8 Ohm refurbished", "https://s.test/x"),
])
def test_b_stock_and_used_items_never_match(title, url):
    listing = Listing(url=url, title=title, manufacturer="B&C", model="12NDL88", price=150.0, currency="GBP",
                      impedance_ohm=8.0)
    assert index().match(listing)[:2] == (None, None)
    new = Listing(url="https://s.test/bc_12ndl88_8_ohms.htm", title="B&C 12NDL88 8 Ohm", manufacturer="B&C",
                  model="12NDL88", price=215.0, currency="GBP", impedance_ohm=8.0)
    assert index().match(new)[:2] == ("bc-12ndl88", "exact")


def test_parts_and_kits_never_match_the_driver():
    kit = Listing(url="https://s.test", title="Celestion Repair Kit for FTX1225", manufacturer="Celestion",
                  model="Repair", price=90.0, currency="EUR", impedance_ohm=8.0)
    assert index().match(kit) == (None, None, "a part or kit for a driver, not the driver")


# -- Thomann -----------------------------------------------------------------------

def test_thomann_category_page_reads_each_product_block():
    shop = ThomannUK(FakeFetcher("thomann-uk"))
    url = "https://www.thomann.co.uk/b-c_loudspeaker_12_inch.html"
    listings = shop.parse_listings(page(url, fixture("thomann_uk_category.html")), ProductRef(url))
    assert [(l.model, l.impedance_ohm, l.price, l.in_stock) for l in listings] == [
        ("12NDL88", 8.0, 215.0, True), ("12PE32", None, 135.0, False), ("12TBX100", None, 205.0, True)]
    assert listings[0].url == "https://www.thomann.co.uk/bc_12ndl88_8_ohms.htm"
    assert listings[0].sku == "641689" and listings[0].currency == "GBP" and listings[0].vat_included is True
    assert ThomannUK.last_page(fixture("thomann_uk_category.html").decode()) == 1


def test_thomann_product_page_and_soft_404():
    shop = ThomannUK(FakeFetcher("thomann-uk"))
    url = "https://www.thomann.co.uk/bc_12pe32.htm"
    listing = shop.parse_listing(page(url, fixture("thomann_uk_product.html")), ProductRef(url))
    assert (listing.manufacturer, listing.model, listing.price, listing.impedance_ohm) == ("B&C", "12PE32", 135.0, 8.0)
    assert listing.title == "B&C 12PE32"
    unknown = "https://www.thomann.co.uk/bc_99xx.htm"
    assert shop.parse_listing(page(unknown, fixture("thomann_uk_unknown_item.html")), ProductRef(unknown)) is None


def test_thomann_run_fetches_product_page_only_for_known_models_without_impedance(tmp_path):
    category = "https://www.thomann.co.uk/b-c_loudspeaker_12_inch.html"
    product = "https://www.thomann.co.uk/bc_12pe32.htm"
    fetcher = FakeFetcher("thomann-uk", {category: (200, fixture("thomann_uk_category.html")),
                                         product: (200, fixture("thomann_uk_product.html"))})

    class OneCategory(ThomannUK):
        def discover(self, limit=None):
            yield ProductRef(category)

    review = tmp_path / "review.csv"
    summary = OneCategory(fetcher).run(index(), tmp_path / "prices", review)
    assert [u for u, _ in fetcher.calls] == [category, product]   # 12TBX100 is not in the index: no fetch
    offers = {o.driver_id: o for o in load_offers(tmp_path / "prices")}
    assert set(offers) == {"bc-12ndl88", "bc-12pe32"} and summary.written == 2
    assert offers["bc-12pe32"].price == 135.0 and offers["bc-12pe32"].region == "UK"
    rows = list(csv.DictReader(review.open()))
    assert [r["model"] for r in rows] == ["12TBX100"] and rows[0]["note"] == "model not in the database"


def test_thomann_de_uses_german_pages_not_intl():
    shop = ThomannDE(FakeFetcher("thomann-de"))
    assert shop.product_url("bc_12pe32.htm?type=manufacturer") == "https://www.thomann.de/de/bc_12pe32.htm"
    assert (shop.region, shop.currency) == ("EU", "EUR")


# -- other shops -----------------------------------------------------------------------

def test_lautsprechershop_blocks_give_variants_with_their_own_price():
    shop = Lautsprechershop(FakeFetcher("lautsprechershop"))
    url = "https://www.lautsprechershop.de/pa/18sound_woofer_en.htm"
    listings = shop.parse_listings(page(url, fixture("lautsprechershop_category.html")), ProductRef(url))
    assert [(l.manufacturer, l.model, l.impedance_ohm, l.price, l.sku) for l in listings] == [
        ("18Sound", "6 NMB 420", 8.0, 205.30, "18s-6nmb420"), ("18Sound", "6 NMB 420", 16.0, 183.90, "18s-6nmb420-16")]
    idx = index()
    assert idx.match(listings[0])[:2] == ("18sound-6nmb420", "exact")
    assert idx.match(listings[1])[1] is None        # 16 ohm variant is not in the database


def test_bax_reads_the_main_product_not_related_items():
    shop = BaxShopUK(FakeFetcher("bax-shop-uk"))
    url = "https://www.bax-shop.co.uk/12-inch-speakers/celestion-ftx1225"
    listing = shop.parse_listing(page(url, fixture("bax_product.html")), ProductRef(url))
    assert (listing.manufacturer, listing.model, listing.price, listing.impedance_ohm, listing.in_stock) == (
        "Celestion", "FTX1225", 249.0, 8.0, False)
    assert index().match(listing)[:2] == ("celestion-ftx1225", "exact")


def test_solen_reads_json_listing_and_never_matches_pair_prices(tmp_path):
    url = "https://www.solen.ca/en/products?brand=beyma&perPage=24&page=1"
    fetcher = FakeFetcher("solen", {url: (200, fixture("solen_beyma_page1.json"))})
    shop = Solen(fetcher)
    listings = shop.parse_listings(page(url, fixture("solen_beyma_page1.json")), ProductRef(url, manufacturer="Beyma"))
    woofer = listings[1]
    assert (woofer.model, woofer.price, woofer.impedance_ohm, woofer.in_stock) == ("12P80Fe", 279.5, 8.0, False)
    assert index().match(woofer)[:2] == ("beyma-12p80fe", "exact")
    assert listings[0].impedance_ohm is None and "per pair" in listings[0].title
    shop.fetch(ProductRef(url))
    assert fetcher.calls[-1] == (url, "application/json")


def test_dayton_msrp_and_madisound_product_pages():
    url = "https://www.daytonaudio.com/product/76/pa310-8-12-pro-woofer-8-ohm"
    listing = DaytonMSRP(FakeFetcher("dayton-msrp")).parse_listing(page(url, fixture("dayton_product.html")), ProductRef(url))
    assert (listing.model, listing.price, listing.impedance_ohm, listing.in_stock) == ("PA310-8", 142.99, 8.0, True)
    assert index().match(listing)[:2] == ("dayton-pa310-8", "exact")
    shop = DaytonMSRP(FakeFetcher("dayton-msrp"))
    assert shop.to_offer(listing, "dayton-pa310-8", "exact", "2026-11-17").price_kind == "list"
    url = "https://www.madisoundspeakerstore.com/special-sale/eminence-alphalite-6b-16-ohm-6.5-pro-midbass/"
    listing = Madisound(FakeFetcher("madisound")).parse_listing(page(url, fixture("madisound_product.html")), ProductRef(url))
    assert (listing.manufacturer, listing.model, listing.price, listing.impedance_ohm) == ("Eminence", "AlphaLite-6B", 19.0, 16.0)


def test_price_coverage_counts_exact_offers_per_region_and_source():
    from horn_drivers.prices import PriceOffer
    from horn_drivers.sources.cli import price_coverage
    records = [record("BC", "12PE32", 8.0, driver_type="cone"), record("BC", "12NDL88", 8.0, driver_type="cone"),
               record("Dayton", "PA310-8", 8.0, driver_type="cone")]
    offer = lambda driver, source, region, currency, confidence="exact": PriceOffer(
        driver, source, "https://s.test", None, 10.0, currency, True, region, None, "2026-11-17", confidence)
    offers = [offer("bc-12pe32", "thomann-uk", "UK", "GBP"), offer("bc-12pe32", "thomann-de", "EU", "EUR"),
              offer("bc-12ndl88", "thomann-uk", "UK", "GBP", "fuzzy"), offer("dayton-pa310-8", "dayton-msrp", "US", "USD")]
    result = price_coverage(records, offers, ["B&C"], "cone")
    assert (result["drivers"], result["priced"], result["in_two_or_more_regions"]) == (2, 1, 1)
    assert result["per_region"] == {"EU": 1, "UK": 1} and result["per_source"] == {"thomann-de": 1, "thomann-uk": 1}
