"""Driver prices in the ranking: an optional filter and a cheapest-feasible marker.

Price never enters the composite score. When a price region is chosen, the
prescreen stores one offer per driver (``horn_drivers.prices.price_snapshot``)
in ``prescreen_result.json`` under ``pricing``. Ranking stages copy that offer
onto each row as ``price`` and mark the cheapest model-feasible design.
"""

from typing import Iterable, List, Optional


def attach_prices(rows: Iterable[dict], pricing: Optional[dict]) -> None:
    """Set ``row["price"]`` to the driver's offer (or None) on every row."""
    if not pricing:
        return
    offers = pricing.get("offers", {})
    for row in rows:
        row["price"] = offers.get(row.get("driver_id"))


def _summary(row: dict) -> dict:
    price = row["price"]
    return {"driver_id": row.get("driver_id"), "horn_label": row.get("horn_label"),
            "manufacturer": row.get("manufacturer"), "model_name": row.get("model_name"),
            "composite_score": row.get("composite_score"), "price": price["price"],
            "currency": price["currency"], "source_id": price["source_id"], "url": price["url"],
            "retrieved_on": price["retrieved_on"]}


def cheapest_feasible(rows: Iterable[dict], previous: Optional[dict] = None) -> Optional[dict]:
    """The cheapest model-feasible priced row (ties: higher score), as a summary.

    *previous* is an earlier cheapest summary (for example from before the
    refinement added rows); it competes on the same terms.
    """
    candidates = [_summary(r) for r in rows if r.get("model_feasible") and r.get("price")]
    if previous:
        candidates.append(previous)
    if not candidates:
        return None
    return min(candidates, key=lambda c: (c["price"], -(c.get("composite_score") or 0.0), c["driver_id"], c["horn_label"] or ""))


def annotate_ranking(feasible: List[dict], shown: List[dict], pricing: Optional[dict],
                     previous: Optional[dict] = None) -> Optional[dict]:
    """Attach prices, mark the cheapest feasible design and return the run's pricing summary.

    *feasible* is every model-feasible row (before the top-N cut), *shown*
    the rows that will be written. Returns None when no price region was chosen.
    """
    if not pricing:
        return None
    attach_prices(feasible, pricing)
    attach_prices(shown, pricing)
    cheapest = cheapest_feasible(feasible, previous)
    for row in shown:
        row["cheapest_feasible"] = bool(cheapest and row.get("driver_id") == cheapest["driver_id"]
                                        and row.get("horn_label") == cheapest["horn_label"])
    return {
        "region": pricing["region"], "currency": pricing["currency"], "as_of": pricing["as_of"],
        "max_price": pricing.get("max_price"), "stale_after_days": pricing.get("stale_after_days"),
        "fx": pricing.get("fx"), "snapshot_dates": pricing.get("snapshot_dates", []),
        "sources": pricing.get("sources", []), "priced_drivers": len(pricing.get("offers", {})),
        "excluded_unpriced": pricing.get("excluded_unpriced", 0),
        "excluded_over_max_price": pricing.get("excluded_over_max_price", 0),
        "cheapest_feasible": cheapest,
        "note": "Price is not part of the score. Prices are as listed by the shop, shipping excluded.",
    }


def format_price(price: Optional[dict]) -> str:
    """'84.00 GBP' style text, with the listed price when it was converted."""
    if not price:
        return "no current offer"
    text = f"{price['price']:.2f} {price['currency']}"
    if price.get("listed_currency") and price["listed_currency"] != price["currency"]:
        text += f" (listed {price['listed_price']:.2f} {price['listed_currency']})"
    return text


def describe_pricing(summary: Optional[dict]) -> List[str]:
    """Plain-text lines for the summary file."""
    if not summary:
        return []
    vat = "VAT as listed by each shop"
    lines = [f"Prices: region {summary['region']}, {summary['currency']}, {vat}, shipping excluded; "
             f"price is not part of the score.",
             f"  Snapshot dates: {', '.join(summary['snapshot_dates']) or 'none'}; sources: "
             f"{', '.join(summary['sources']) or 'none'}; offers older than {summary.get('stale_after_days')} days "
             f"on {summary['as_of']} are not used."]
    if summary.get("fx"):
        lines.append(f"  Converted with ECB reference rates of {summary['fx']['rate_date']}.")
    if summary.get("max_price") is not None:
        lines.append(f"  Max price {summary['max_price']:g} {summary['currency']}: excluded "
                     f"{summary['excluded_over_max_price']} over the limit and {summary['excluded_unpriced']} "
                     "without a current offer.")
    cheapest = summary.get("cheapest_feasible")
    if cheapest:
        lines.append(f"  Cheapest feasible: {cheapest['horn_label']} + {cheapest['driver_id']} at "
                     f"{cheapest['price']:.2f} {cheapest['currency']} ({cheapest['source_id']}, "
                     f"{cheapest['retrieved_on']}), score {cheapest['composite_score']:.3f}.")
    else:
        lines.append("  Cheapest feasible: no feasible design has a current price.")
    return lines
