"""Prices in the ranking: max_price filter, cheapest feasible marker and report output."""

import json
from pathlib import Path

import pytest

from horn_drivers.loader import load_driver
from horn_drivers.prices import PriceOffer, append_offers
from horn_analysis.auto_report import generate_auto_report
from horn_analysis.pricing import annotate_ranking, cheapest_feasible, describe_pricing, format_price
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "examples/800-1600-modal"


def offer(price, currency="GBP", source="thomann-uk", url="https://shop.test/a"):
    return {"price": price, "currency": "GBP", "listed_price": price, "listed_currency": currency, "fx_rate": None,
            "vat_included": True, "in_stock": True, "source_id": source, "url": url, "sku": "1",
            "retrieved_on": "2026-11-10", "match_confidence": "exact"}


PRICING = {"region": "UK", "currency": "GBP", "as_of": "2026-11-17", "stale_after_days": 30, "max_price": None,
           "fx": None, "offers": {"a": offer(200.0), "b": offer(80.0)}, "unpriced": ["c"],
           "snapshot_dates": ["2026-11-10"], "sources": ["thomann-uk"]}


def rows():
    return [{"driver_id": "a", "horn_label": "h1", "composite_score": 0.9, "model_feasible": True},
            {"driver_id": "c", "horn_label": "h1", "composite_score": 0.8, "model_feasible": True},
            {"driver_id": "b", "horn_label": "h2", "composite_score": 0.5, "model_feasible": True}]


def test_price_is_attached_and_cheapest_feasible_found_beyond_the_shown_rows():
    feasible = rows()
    shown = feasible[:2]                       # the cheapest row is outside the top N
    summary = annotate_ranking(feasible, shown, PRICING)
    assert shown[0]["price"]["price"] == 200.0 and shown[1]["price"] is None
    assert summary["cheapest_feasible"]["driver_id"] == "b"
    assert summary["cheapest_feasible"]["price"] == 80.0
    assert not any(r["cheapest_feasible"] for r in shown)
    # Composite scores are untouched: price is not scored.
    assert [r["composite_score"] for r in shown] == [0.9, 0.8]
    assert annotate_ranking(rows(), rows(), None) is None


def test_previous_cheapest_survives_refinement():
    previous = {"driver_id": "z", "horn_label": "h0", "price": 10.0, "composite_score": 0.1, "currency": "GBP"}
    assert cheapest_feasible(rows(), previous)["driver_id"] == "z"
    feasible = rows()
    feasible[0]["model_feasible"] = False
    feasible[2]["model_feasible"] = False
    for row in feasible:
        row["price"] = PRICING["offers"].get(row["driver_id"])
    assert cheapest_feasible(feasible) is None


def test_price_text_shows_conversion_and_summary_lines():
    converted = {**offer(84.0, currency="EUR"), "listed_price": 100.0}
    assert format_price(converted) == "84.00 GBP (listed 100.00 EUR)"
    assert format_price(None) == "no current offer"
    summary = annotate_ranking(rows(), rows(), {**PRICING, "max_price": 150.0, "excluded_unpriced": 3,
                                                "excluded_over_max_price": 2})
    text = "\n".join(describe_pricing(summary))
    assert "region UK" in text and "2026-11-10" in text and "not part of the score" in text
    assert "Max price 150 GBP: excluded 2 over the limit and 3 without a current offer" in text
    assert "Cheapest feasible: h2 + b at 80.00 GBP" in text


def _write_db(tmp_path):
    db = tmp_path / "db" / "Acme"
    db.mkdir(parents=True)
    base = {"fs_hz": 300.0, "re_ohm": 6.0, "bl_tm": 8.0, "sd_m2": 0.002, "mms_kg": 0.003, "le_h": 0.0002,
            "qms": 5.0, "qes": 0.4}
    for name in ("cheap", "dear", "unpriced"):
        (db / f"{name}.json").write_text(json.dumps({"driver_id": name, "manufacturer": "Acme", "model_name": name,
                                                     "driver_type": "cone", "parameters": base}))
    prices = tmp_path / "prices"
    make = lambda d, p: PriceOffer(driver_id=d, source_id="thomann-uk", url=f"https://shop.test/{d}", sku=None,
                                   price=p, currency="GBP", vat_included=True, region="UK", in_stock=True,
                                   retrieved_on="2026-11-10", match_confidence="exact")
    append_offers(prices, "thomann-uk", [make("cheap", 50.0), make("dear", 500.0)], "2026-11-10")
    return tmp_path / "db", prices


def _prescreen(monkeypatch, tmp_path, *extra):
    from horn_analysis import prescreen
    db, prices = _write_db(tmp_path)
    out = tmp_path / "prescreen.json"
    monkeypatch.setattr("sys.argv", ["prescreen", "--drivers-db", str(db), "--target-f-low", "800",
                                     "--target-f-high", "1600", "--output", str(out), "--prices-dir", str(prices),
                                     "--price-as-of", "2026-11-17", *extra])
    prescreen.main()
    return json.loads(out.read_text())


def test_prescreen_max_price_keeps_only_drivers_priced_under_the_limit(monkeypatch, tmp_path):
    screen = _prescreen(monkeypatch, tmp_path, "--price-region", "UK", "--max-price", "100")
    assert screen["drivers"] == ["cheap"]
    pricing = screen["pricing"]
    assert (pricing["excluded_over_max_price"], pricing["excluded_unpriced"]) == (1, 1)
    assert set(pricing["offers"]) == {"cheap"} and pricing["offers"]["cheap"]["price"] == 50.0


def test_prescreen_region_without_limit_keeps_every_driver(monkeypatch, tmp_path):
    screen = _prescreen(monkeypatch, tmp_path, "--price-region", "UK")
    assert sorted(screen["drivers"]) == ["cheap", "dear", "unpriced"]
    assert screen["pricing"]["unpriced"] == ["unpriced"]


def test_prescreen_without_region_has_no_pricing(monkeypatch, tmp_path):
    assert "pricing" not in _prescreen(monkeypatch, tmp_path)
    with pytest.raises(SystemExit):
        _prescreen(monkeypatch, tmp_path / "x", "--max-price", "100")


def test_report_shows_price_sensitivity_and_cheapest_feasible(tmp_path):
    preserved = json.loads((EXAMPLE / "ranking.json").read_text())[0]
    driver = load_driver(str(ROOT / "data/drivers-curated"), preserved["driver_id"])
    target = TargetSpec(800, 1600)
    label = preserved["horn_label"]
    ranked = rank_horn_drivers(str(EXAMPLE / "horn-impedance.csv"), label, preserved["throat_radius"], [driver],
                               target, top_n=1, length=preserved["length"])
    for row in ranked:
        row.update({k: preserved[k] for k in ("profile", "throat_radius", "mouth_radius", "length")})
        row["model_feasible"] = True   # exercise the price path whatever the current verdict
    pricing = {**PRICING, "offers": {driver.driver_id: offer(123.0)}}
    summary = annotate_ranking(ranked, ranked, pricing)
    generate_auto_report(ranked, {label: str(EXAMPLE / "horn-impedance.csv")}, {driver.driver_id: driver}, None,
                         target, str(tmp_path), top_n=5, pricing=summary)
    text = (tmp_path / "auto_summary.txt").read_text()
    assert "Prices: region UK, GBP" in text
    assert "Price: 123.00 GBP (thomann-uk, 2026-11-10, https://shop.test/a); cheapest feasible" in text
    assert "Driver rating:" in text
    page = (tmp_path / "auto_report.html").read_text()
    assert "Price (as listed, shipping excluded)" in page and "cheapest feasible" in page
    assert "Sensitivity (dB, 1 W, 1 m)" in page
    ranking = json.loads((tmp_path / "auto_ranking.json").read_text())
    assert ranking[0]["price"]["price"] == 123.0 and ranking[0]["cheapest_feasible"] is True


def test_report_without_pricing_has_no_price_column(tmp_path):
    from horn_analysis.html_report import _render_pricing, _render_rankings_rows
    assert _render_pricing(None) == ""
    assert "Price" not in _render_rankings_rows([{"driver_id": "a", "horn_label": "h"}], {}, show_price=False)
