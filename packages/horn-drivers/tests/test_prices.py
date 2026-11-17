"""Price snapshots: storage, latest non-stale offer, ECB conversion, validation."""

import datetime as dt
import json

import pytest

from horn_drivers.prices import (PriceOffer, append_offers, latest_offers, load_fx, parse_ecb_xml, price_snapshot,
                                 read_offers, validate_prices, write_fx)

DAY = dt.date(2026, 11, 17)
ECB_XML = """<?xml version="1.0" encoding="UTF-8"?>
<gesmes:Envelope xmlns:gesmes="http://www.gesmes.org/xml/2002-08-01" xmlns="http://www.ecb.int/vocabulary/2002-08-01/eurofxref">
<gesmes:subject>Reference rates</gesmes:subject>
<Cube><Cube time="2026-11-16"><Cube currency="USD" rate="1.0800"/><Cube currency="GBP" rate="0.8400"/></Cube></Cube>
</gesmes:Envelope>"""


def offer(driver="bc-8ndl51", source="thomann-uk", price=100.0, currency="GBP", region="UK",
          day="2026-11-10", match="exact", in_stock=True):
    return PriceOffer(driver_id=driver, source_id=source, url=f"https://{source}.test/{driver}", sku="123",
                      price=price, currency=currency, vat_included=True, region=region, in_stock=in_stock,
                      retrieved_on=day, match_confidence=match)


def test_offer_fields_are_checked():
    with pytest.raises(ValueError, match="price"):
        offer(price=-1)
    with pytest.raises(ValueError, match="currency"):
        offer(currency="pounds")
    with pytest.raises(ValueError, match="region"):
        offer(region="Mars")
    with pytest.raises(ValueError, match="match_confidence"):
        offer(match="maybe")
    with pytest.raises(ValueError, match="unknown offer fields"):
        PriceOffer.from_json({**offer().to_json(), "description": "copied text"})


def test_append_writes_jsonl_per_source_and_day_without_repeats(tmp_path):
    path = append_offers(tmp_path, "thomann-uk", [offer(), offer()], "2026-11-10")
    append_offers(tmp_path, "thomann-uk", [offer()], "2026-11-10")
    assert path == tmp_path / "thomann-uk" / "2026-11-10.jsonl"
    assert len(path.read_text().splitlines()) == 1
    assert json.loads(path.read_text())["currency"] == "GBP"
    with pytest.raises(ValueError):
        append_offers(tmp_path, "other", [offer()])
    with pytest.raises(ValueError):
        append_offers(tmp_path, "fx", [])


def test_latest_offer_per_source_wins_and_stale_is_dropped():
    offers = [offer(price=120, day="2026-11-01"), offer(price=110, day="2026-11-10"),
              offer(source="bax-uk", price=90, day="2026-09-01"),            # stale
              offer(source="madisound", region="US", currency="USD", price=99, day="2026-11-12"),
              offer(source="fuzzy-shop", price=50, match="fuzzy")]
    latest = latest_offers(offers, as_of=DAY)
    assert [o.price for o in latest[("bc-8ndl51", "UK")]] == [110]
    assert [o.price for o in latest[("bc-8ndl51", "US")]] == [99]
    # A stale newest offer is not replaced by an older one from the same source.
    assert ("x", "UK") not in latest_offers([offer(driver="x", day="2026-09-01"), offer(driver="x", day="2026-08-01")], as_of=DAY)
    # Offers dated after the ranking date are ignored.
    assert latest_offers([offer(day="2026-11-20")], as_of=DAY) == {}


def test_ecb_rates_parse_store_and_convert(tmp_path):
    rate_date, rates = parse_ecb_xml(ECB_XML)
    assert rate_date == "2026-11-16" and rates == {"USD": 1.08, "GBP": 0.84}
    write_fx(tmp_path, rate_date, rates, "2026-11-17")
    fx = load_fx(tmp_path, DAY)
    assert fx.rate_date == "2026-11-16"
    assert fx.convert(108, "USD", "GBP") == pytest.approx(84.0)
    assert fx.convert(84, "GBP", "EUR") == pytest.approx(100.0)
    assert load_fx(tmp_path, dt.date(2026, 11, 1)) is None


def test_price_snapshot_picks_cheapest_in_region_currency(tmp_path):
    append_offers(tmp_path, "thomann-uk", [offer(price=110), offer(driver="cheap", price=40)], "2026-11-10")
    append_offers(tmp_path, "eu-shop", [offer(source="eu-shop", price=100, currency="EUR")], "2026-11-12")
    append_offers(tmp_path, "oos-shop", [offer(source="oos-shop", price=10, in_stock=False)], "2026-11-12")
    write_fx(tmp_path, *parse_ecb_xml(ECB_XML), "2026-11-17")
    snap = price_snapshot(tmp_path, ["bc-8ndl51", "cheap", "none"], "UK", as_of=DAY, max_price=90)
    best = snap["offers"]["bc-8ndl51"]
    # 100 EUR = 84 GBP beats 110 GBP; the 10 GBP offer is out of stock.
    assert (best["price"], best["currency"], best["listed_currency"], best["source_id"]) == (84.0, "GBP", "EUR", "eu-shop")
    assert best["fx_rate"] == pytest.approx(0.84)
    assert snap["fx"]["rate_date"] == "2026-11-16"
    assert snap["unpriced"] == ["none"] and snap["over_max_price"] == []
    assert price_snapshot(tmp_path, ["bc-8ndl51"], "UK", as_of=DAY, max_price=50)["over_max_price"] == ["bc-8ndl51"]
    with pytest.raises(ValueError):
        price_snapshot(tmp_path, [], "Atlantis")


def test_without_rates_foreign_offers_are_skipped(tmp_path):
    append_offers(tmp_path, "eu-shop", [offer(source="eu-shop", currency="EUR")], "2026-11-12")
    snap = price_snapshot(tmp_path, ["bc-8ndl51"], "UK", as_of=DAY)
    assert snap["unpriced"] == ["bc-8ndl51"] and snap["fx"] is None


def test_validator_warns_on_stale_unknown_and_fuzzy_and_errors_on_bad_lines(tmp_path):
    append_offers(tmp_path, "thomann-uk", [offer(day="2026-09-01")], "2026-09-01")
    append_offers(tmp_path, "thomann-uk", [offer(driver="ghost", match="fuzzy", day="2026-11-10")], "2026-11-10")
    (tmp_path / "thomann-uk" / "2026-11-11.jsonl").write_text('{"driver_id": "x"}\n')
    errors, warnings = validate_prices(tmp_path, known_driver_ids={"bc-8ndl51"}, as_of=DAY)
    assert len(errors) == 1 and "missing offer fields" in errors[0]
    assert any("older than 30 days" in w for w in warnings)
    assert any("unknown driver ghost" in w for w in warnings)
    assert any("fuzzy" in w for w in warnings)
    offers, bad = read_offers(tmp_path)
    assert len(offers) == 2 and len(bad) == 1


def test_validator_cli_checks_prices(tmp_path, monkeypatch, capsys):
    from horn_drivers import validator
    db = tmp_path / "db"
    (db / "BC").mkdir(parents=True)
    (db / "BC" / "x.json").write_text(json.dumps({"driver_id": "x", "manufacturer": "BC", "model_name": "X", "parameters": {
        "fs_hz": 61.0, "re_ohm": 5.1, "bl_tm": 12.3, "sd_m2": 0.022, "mms_kg": 0.022, "le_h": 0.0004, "qms": 6.0, "qes": 0.36}}))
    prices = tmp_path / "prices"
    append_offers(prices, "thomann-uk", [offer(driver="x", day="2026-01-01")], "2026-01-01")
    monkeypatch.setattr("sys.argv", ["horn-validate-drivers", "--db", str(db), "--prices", str(prices)])
    with pytest.raises(SystemExit) as exit_info:
        validator.main()
    assert exit_info.value.code == 0
    assert "older than 30 days" in capsys.readouterr().out


def test_list_price_is_used_only_without_a_shop_offer(tmp_path):
    cheap_list = PriceOffer(**{**offer(driver="d1", source="maker-msrp", price=50).to_json(), "price_kind": "list"})
    append_offers(tmp_path, "maker-msrp", [cheap_list, PriceOffer(**{**cheap_list.to_json(), "driver_id": "d2"})],
                  "2026-11-10")
    append_offers(tmp_path, "thomann-uk", [offer(driver="d1", price=80)], "2026-11-10")
    snap = price_snapshot(tmp_path, ["d1", "d2"], "UK", as_of=DAY)
    assert (snap["offers"]["d1"]["price"], snap["offers"]["d1"]["price_kind"]) == (80, "offer")
    assert (snap["offers"]["d2"]["price"], snap["offers"]["d2"]["price_kind"]) == (50, "list")
    with pytest.raises(ValueError, match="price_kind"):
        PriceOffer(**{**offer().to_json(), "price_kind": "sale"})


def test_line_without_price_kind_loads_as_offer():
    data = offer().to_json()
    del data["price_kind"]
    assert PriceOffer.from_json(data).price_kind == "offer"
