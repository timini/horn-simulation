"""Driver price snapshots, kept outside the driver records.

Layout::

    data/prices/<source_id>/<YYYY-MM-DD>.jsonl   one offer per line
    data/prices/fx/ecb-<YYYY-MM-DD>.json         ECB euro reference rates

An offer stores the price exactly as listed: currency, whether VAT is
included, the region the shop sells to, stock state, retrieval date and how
sure the match to ``driver_id`` is. Shipping is not included. Conversion
happens only at ranking time (``price_snapshot``), with the stored ECB rate
file, and the rate date is recorded next to the converted price.

Only ``exact`` and ``manual`` matches are used for ranking; ``fuzzy``
matches are kept for review. An offer older than ``STALE_AFTER_DAYS`` on the
ranking date is not used, and the validator warns about it.
"""

import datetime as _dt
import json
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

STALE_AFTER_DAYS = 30
MATCH_CONFIDENCE = ("exact", "manual", "fuzzy")
RANKABLE_MATCHES = ("exact", "manual")
# "offer": a shop's selling price. "list": a maker's list price (MSRP), used
# for a driver only when no shop offer for it exists in the region.
PRICE_KINDS = ("offer", "list")
FX_DIR = "fx"
ECB_DAILY_URL = "https://www.ecb.europa.eu/stats/eurofxref/eurofxref-daily.xml"
# Region code -> the currency a buyer there pays in. Prices for a region are
# compared and filtered in this currency.
REGION_CURRENCY = {
    "UK": "GBP", "EU": "EUR", "US": "USD", "CA": "CAD", "AU": "AUD",
    "CH": "CHF", "JP": "JPY", "NO": "NOK", "SE": "SEK", "DK": "DKK", "PL": "PLN",
}
OFFER_FIELDS = ("driver_id", "source_id", "url", "sku", "price", "currency", "vat_included",
                "region", "in_stock", "retrieved_on", "match_confidence", "price_kind")
# Fields an older snapshot line may lack, with the value it then takes.
OFFER_DEFAULTS = {"price_kind": "offer"}


@dataclass(frozen=True)
class PriceOffer:
    driver_id: str
    source_id: str
    url: str
    sku: Optional[str]
    price: float
    currency: str
    vat_included: Optional[bool]
    region: str
    in_stock: Optional[bool]
    retrieved_on: str            # YYYY-MM-DD
    match_confidence: str        # exact | manual | fuzzy
    price_kind: str = "offer"    # offer | list

    def __post_init__(self):
        errors = offer_errors(asdict(self))
        if errors:
            raise ValueError(f"invalid price offer for {self.driver_id!r}: " + "; ".join(errors))

    def to_json(self) -> dict:
        return asdict(self)

    @classmethod
    def from_json(cls, data: dict) -> "PriceOffer":
        unknown = set(data) - set(OFFER_FIELDS)
        if unknown:
            raise ValueError(f"unknown offer fields: {', '.join(sorted(unknown))}")
        missing = [name for name in OFFER_FIELDS if name not in data and name not in OFFER_DEFAULTS]
        if missing:
            raise ValueError(f"missing offer fields: {', '.join(missing)}")
        return cls(**{name: data.get(name, OFFER_DEFAULTS.get(name)) for name in OFFER_FIELDS})

    @property
    def date(self) -> _dt.date:
        return _dt.date.fromisoformat(self.retrieved_on)


def offer_errors(data: dict) -> List[str]:
    """Field-level problems with one offer dict (empty when valid)."""
    errors = []
    for name in ("driver_id", "source_id", "url", "currency", "region", "retrieved_on", "match_confidence"):
        if not isinstance(data.get(name), str) or not data.get(name):
            errors.append(f"{name} must be a nonempty string")
    price = data.get("price")
    if isinstance(price, bool) or not isinstance(price, (int, float)) or not math.isfinite(price) or price <= 0:
        errors.append("price must be a finite positive number")
    if isinstance(data.get("currency"), str) and not re.fullmatch(r"[A-Z]{3}", data["currency"]):
        errors.append("currency must be an ISO 4217 code such as GBP")
    if isinstance(data.get("region"), str) and data["region"] not in REGION_CURRENCY:
        errors.append(f"region must be one of {', '.join(REGION_CURRENCY)}")
    if isinstance(data.get("url"), str) and not data["url"].startswith(("http://", "https://")):
        errors.append("url must be http(s)")
    for name in ("vat_included", "in_stock"):
        if data.get(name) is not None and not isinstance(data.get(name), bool):
            errors.append(f"{name} must be true, false or null")
    if data.get("sku") is not None and not isinstance(data.get("sku"), str):
        errors.append("sku must be a string or null")
    if data.get("match_confidence") not in MATCH_CONFIDENCE:
        errors.append(f"match_confidence must be one of {', '.join(MATCH_CONFIDENCE)}")
    if data.get("price_kind", "offer") not in PRICE_KINDS:
        errors.append(f"price_kind must be one of {', '.join(PRICE_KINDS)}")
    try:
        _dt.date.fromisoformat(str(data.get("retrieved_on")))
    except ValueError:
        errors.append("retrieved_on must be YYYY-MM-DD")
    return errors


# -- files -------------------------------------------------------------------

def snapshot_path(prices_dir: Path, source_id: str, day: str) -> Path:
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", source_id) or source_id == FX_DIR:
        raise ValueError(f"source_id {source_id!r} is not a safe folder name")
    _dt.date.fromisoformat(day)
    return Path(prices_dir) / source_id / f"{day}.jsonl"


def append_offers(prices_dir: Path, source_id: str, offers: Iterable[PriceOffer], day: Optional[str] = None) -> Path:
    """Append offers to ``<prices_dir>/<source_id>/<day>.jsonl``, skipping exact repeats."""
    offers = list(offers)
    if any(offer.source_id != source_id for offer in offers):
        raise ValueError("every offer must come from the source it is filed under")
    day = day or (max(offer.retrieved_on for offer in offers) if offers else _dt.date.today().isoformat())
    path = snapshot_path(prices_dir, source_id, day)
    path.parent.mkdir(parents=True, exist_ok=True)
    existing = set(path.read_text().splitlines()) if path.exists() else set()
    with path.open("a") as handle:
        for offer in offers:
            line = json.dumps(offer.to_json(), sort_keys=True)
            if line not in existing:
                handle.write(line + "\n")
                existing.add(line)
    return path


def _snapshot_files(prices_dir: Path):
    root = Path(prices_dir)
    if not root.is_dir():
        return []
    return sorted(p for p in root.glob("*/*.jsonl") if p.parent.name != FX_DIR)


def read_offers(prices_dir: Path) -> Tuple[List[PriceOffer], List[str]]:
    """All offers under *prices_dir*, plus a list of problems with bad lines."""
    offers, errors = [], []
    for path in _snapshot_files(prices_dir):
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if not line.strip():
                continue
            where = f"{path.parent.name}/{path.name}:{number}"
            try:
                offer = PriceOffer.from_json(json.loads(line))
            except (ValueError, TypeError) as error:
                errors.append(f"{where}: {error}")
                continue
            if offer.source_id != path.parent.name:
                errors.append(f"{where}: source_id {offer.source_id!r} filed under {path.parent.name!r}")
                continue
            offers.append(offer)
    return offers, errors


def load_offers(prices_dir: Path) -> List[PriceOffer]:
    """All valid offers; raises ValueError when any line is invalid."""
    offers, errors = read_offers(prices_dir)
    if errors:
        raise ValueError("invalid price snapshot lines:\n  " + "\n  ".join(errors))
    return offers


def is_stale(offer: PriceOffer, as_of: _dt.date, max_age_days: int = STALE_AFTER_DAYS) -> bool:
    return (as_of - offer.date).days > max_age_days


def latest_offers(offers: Iterable[PriceOffer], *, as_of: Optional[_dt.date] = None,
                  max_age_days: int = STALE_AFTER_DAYS,
                  accepted: Tuple[str, ...] = RANKABLE_MATCHES) -> Dict[Tuple[str, str], List[PriceOffer]]:
    """Latest usable offer per (driver_id, region), one per source.

    For each driver, region and source the most recent offer is taken (ties
    on the same day: the lowest price). It is dropped when it is stale on
    *as_of*, dated after *as_of*, or its match is not in *accepted*; an
    older offer from the same source never stands in for a stale newer one.
    """
    as_of = as_of or _dt.date.today()
    newest: Dict[Tuple[str, str, str], PriceOffer] = {}
    for offer in offers:
        if offer.date > as_of:
            continue
        key = (offer.driver_id, offer.region, offer.source_id)
        current = newest.get(key)
        if current is None or (offer.date, -offer.price) > (current.date, -current.price):
            newest[key] = offer
    result: Dict[Tuple[str, str], List[PriceOffer]] = {}
    for (driver_id, region, _), offer in sorted(newest.items()):
        if offer.match_confidence in accepted and not is_stale(offer, as_of, max_age_days):
            result.setdefault((driver_id, region), []).append(offer)
    return result


# -- exchange rates ------------------------------------------------------------

@dataclass(frozen=True)
class FxRates:
    """Euro reference rates: 1 EUR = rates[code] units of that currency."""
    rate_date: str
    rates: Dict[str, float]
    source_url: str = ECB_DAILY_URL
    retrieved_on: Optional[str] = None
    path: Optional[str] = None

    def convert(self, amount: float, from_currency: str, to_currency: str) -> float:
        if from_currency == to_currency:
            return float(amount)
        table = {"EUR": 1.0, **self.rates}
        if from_currency not in table or to_currency not in table:
            raise KeyError(f"no ECB rate for {from_currency} or {to_currency} on {self.rate_date}")
        return float(amount) / table[from_currency] * table[to_currency]

    def rate(self, from_currency: str, to_currency: str) -> float:
        return self.convert(1.0, from_currency, to_currency)


def parse_ecb_xml(text: str) -> Tuple[str, Dict[str, float]]:
    """(rate date, {currency: rate}) from the ECB daily reference rate XML."""
    root = ET.fromstring(text)
    dated = [el for el in root.iter() if el.tag.endswith("Cube") and el.get("time")]
    if len(dated) != 1:
        raise ValueError("expected one dated Cube element in the ECB file")
    rates = {el.get("currency"): float(el.get("rate")) for el in dated[0]
             if el.tag.endswith("Cube") and el.get("currency")}
    if not rates:
        raise ValueError("ECB file holds no rates")
    _dt.date.fromisoformat(dated[0].get("time"))
    return dated[0].get("time"), rates


def write_fx(prices_dir: Path, rate_date: str, rates: Dict[str, float], retrieved_on: str,
             source_url: str = ECB_DAILY_URL) -> Path:
    path = Path(prices_dir) / FX_DIR / f"ecb-{rate_date}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"source": "European Central Bank euro foreign exchange reference rates",
                                "source_url": source_url, "base": "EUR", "rate_date": rate_date,
                                "retrieved_on": retrieved_on, "rates": dict(sorted(rates.items()))},
                               indent=2) + "\n")
    return path


def load_fx(prices_dir: Path, as_of: Optional[_dt.date] = None) -> Optional[FxRates]:
    """The newest stored ECB rate file dated on or before *as_of*, or None."""
    as_of = as_of or _dt.date.today()
    best = None
    for path in sorted((Path(prices_dir) / FX_DIR).glob("ecb-*.json")):
        data = json.loads(path.read_text())
        if _dt.date.fromisoformat(data["rate_date"]) <= as_of:
            best = FxRates(rate_date=data["rate_date"], rates={k: float(v) for k, v in data["rates"].items()},
                           source_url=data.get("source_url", ECB_DAILY_URL),
                           retrieved_on=data.get("retrieved_on"), path=str(path))
    return best


def fetch_ecb_rates(fetcher, prices_dir: Path, url: str = ECB_DAILY_URL) -> Path:
    """Fetch the ECB daily XML once through a PoliteFetcher and store the rates."""
    result = fetcher.get(url, max_age_days=1.0)
    if not result.ok:
        raise RuntimeError(f"ECB rates not available: HTTP {result.status}")
    rate_date, rates = parse_ecb_xml(result.text)
    return write_fx(prices_dir, rate_date, rates, result.retrieved_on, url)


# -- ranking-time view ---------------------------------------------------------

def price_snapshot(prices_dir: Path, driver_ids: Iterable[str], region: str, *,
                   as_of: Optional[_dt.date] = None, max_price: Optional[float] = None) -> dict:
    """Cheapest usable offer per driver in *region*, in the region's currency.

    Returns a JSON-ready dict::

        {"region", "currency", "as_of", "stale_after_days", "max_price",
         "fx": {"source_url", "rate_date"} | None,
         "offers": {driver_id: {"price", "currency", "listed_price", "listed_currency",
                                "fx_rate", "vat_included", "in_stock", "source_id",
                                "url", "sku", "retrieved_on", "match_confidence", "price_kind"}},
         "unpriced": [driver_id, ...], "over_max_price": [driver_id, ...],
         "snapshot_dates": [YYYY-MM-DD, ...], "sources": [source_id, ...],
         "invalid_lines_skipped": int}

    A shop offer is preferred over a maker's list price, and an in-stock
    offer over a cheaper one that is out of stock (stock unknown counts as
    in stock). Offers in another currency are
    converted with the stored ECB rate; without a rate file they are skipped
    and the driver may end up unpriced.
    """
    if region not in REGION_CURRENCY:
        raise ValueError(f"price region must be one of {', '.join(REGION_CURRENCY)}")
    if max_price is not None and not (math.isfinite(max_price) and max_price > 0):
        raise ValueError("max_price must be finite and positive")
    as_of = as_of or _dt.date.today()
    currency = REGION_CURRENCY[region]
    fx = load_fx(prices_dir, as_of)
    all_offers, bad_lines = read_offers(prices_dir)
    latest = latest_offers(all_offers, as_of=as_of)
    offers, unpriced, over = {}, [], []
    for driver_id in sorted(set(driver_ids)):
        candidates = []
        for offer in latest.get((driver_id, region), []):
            try:
                converted = offer.price if offer.currency == currency else fx.convert(offer.price, offer.currency, currency)
            except (AttributeError, KeyError):
                continue
            candidates.append((offer.price_kind == "list", offer.in_stock is False, converted, offer.source_id, offer))
        if not candidates:
            unpriced.append(driver_id)
            continue
        _, _, converted, _, offer = min(candidates)
        offers[driver_id] = {
            "price": round(converted, 2), "currency": currency,
            "listed_price": offer.price, "listed_currency": offer.currency,
            "fx_rate": None if offer.currency == currency else fx.rate(offer.currency, currency),
            "vat_included": offer.vat_included, "in_stock": offer.in_stock,
            "source_id": offer.source_id, "url": offer.url, "sku": offer.sku,
            "retrieved_on": offer.retrieved_on, "match_confidence": offer.match_confidence,
            "price_kind": offer.price_kind,
        }
        if max_price is not None and converted > max_price:
            over.append(driver_id)
    used_fx = any(entry["fx_rate"] is not None for entry in offers.values())
    return {
        "region": region, "currency": currency, "as_of": as_of.isoformat(),
        "stale_after_days": STALE_AFTER_DAYS, "max_price": max_price,
        "fx": {"source_url": fx.source_url, "rate_date": fx.rate_date} if fx and used_fx else None,
        "offers": offers, "unpriced": unpriced, "over_max_price": over,
        "snapshot_dates": sorted({entry["retrieved_on"] for entry in offers.values()}),
        "sources": sorted({entry["source_id"] for entry in offers.values()}),
        "invalid_lines_skipped": len(bad_lines),
    }


def validate_prices(prices_dir: Path, known_driver_ids: Optional[Iterable[str]] = None,
                    as_of: Optional[_dt.date] = None) -> Tuple[List[str], List[str]]:
    """(errors, warnings) for every snapshot file.

    Errors: unreadable or invalid lines. Warnings: offers whose latest
    snapshot is older than STALE_AFTER_DAYS, offers naming drivers not in
    the database, and fuzzy matches still awaiting review.
    """
    as_of = as_of or _dt.date.today()
    offers, errors = read_offers(prices_dir)
    warnings: List[str] = []
    known = set(known_driver_ids) if known_driver_ids is not None else None
    newest: Dict[Tuple[str, str, str], PriceOffer] = {}
    for offer in offers:
        key = (offer.driver_id, offer.region, offer.source_id)
        if key not in newest or offer.date > newest[key].date:
            newest[key] = offer
    for (driver_id, region, source_id), offer in sorted(newest.items()):
        if is_stale(offer, as_of):
            warnings.append(f"{source_id} {region} price for {driver_id} is from {offer.retrieved_on}, "
                            f"older than {STALE_AFTER_DAYS} days")
        if offer.date > as_of:
            warnings.append(f"{source_id} {region} price for {driver_id} is dated in the future ({offer.retrieved_on})")
        if known is not None and driver_id not in known:
            warnings.append(f"{source_id} price names unknown driver {driver_id}")
        if offer.match_confidence == "fuzzy":
            warnings.append(f"{source_id} {region} price for {driver_id} is a fuzzy match awaiting review")
    fx_dir = Path(prices_dir) / FX_DIR
    for path in sorted(fx_dir.glob("*.json")) if fx_dir.is_dir() else []:
        try:
            data = json.loads(path.read_text())
            _dt.date.fromisoformat(data["rate_date"])
            if not data["rates"] or not all(float(v) > 0 for v in data["rates"].values()):
                raise ValueError("rates must be positive")
        except (ValueError, KeyError, TypeError) as error:
            errors.append(f"{FX_DIR}/{path.name}: {error}")
    return errors, warnings
