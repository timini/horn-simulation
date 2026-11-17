"""Source adapter interface.

An adapter turns one site into facts. It never touches the network itself:
``discover`` yields product references, the base class fetches each through
the shared ``PoliteFetcher`` (robots.txt, pacing, cache, stop on blocks), and
the adapter's parse method reads the fetched page. Parse methods are pure, so
they are tested offline against cached pages and rerun after a parser fix
without refetching.

Two kinds:

``DriverSource`` (manufacturer or aggregator)
    implements ``discover`` and ``parse_driver``; ``parse_driver`` returns a
    ``DriverResult`` whose ``record`` is a schema version 2 record built with
    ``horn_drivers.records.make_record`` (SI units, field_provenance for each
    field), or ``None`` with a ``reason``.

``PriceSource`` (retailer)
    implements ``discover`` and ``parse_listing``; ``parse_listing`` returns a
    ``Listing`` (facts read from the page) or ``None`` for a soft 404 or a
    page without one clear product. The base class matches listings to
    driver records with ``DriverIndex`` and builds ``PriceOffer`` lines.

Store only facts: numbers, units, URL, date, SKU. No descriptions or images.
"""

import csv
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar, Dict, Iterable, Iterator, List, Optional, Tuple

from horn_drivers.fetch import DisallowedByRobots, FetchResult, PoliteFetcher, ResponseTooLarge, SourceBlocked
from horn_drivers.identity import canonical_manufacturer, normalise_model, record_identity
from horn_drivers.prices import PriceOffer, append_offers
from horn_drivers.records import check_structure, write_incomplete, write_record


@dataclass(frozen=True)
class ProductRef:
    """A product page found by discovery."""
    url: str
    manufacturer: Optional[str] = None     # as the source names it
    model: Optional[str] = None
    sku: Optional[str] = None
    hints: Tuple[Tuple[str, str], ...] = ()  # adapter-specific context, e.g. (("category", "woofer"),)

    def hint(self, key: str, default: Optional[str] = None) -> Optional[str]:
        return dict(self.hints).get(key, default)


@dataclass
class DriverResult:
    ref: ProductRef
    record: Optional[dict]
    reason: Optional[str] = None           # why record is None, or why it is incomplete
    incomplete: bool = False               # True: write to drivers-incomplete with reason


@dataclass(frozen=True)
class Listing:
    """Facts read from one retailer product page."""
    url: str
    title: str                             # page or product-block title, used for the match check
    manufacturer: str
    model: str
    price: float
    currency: str
    sku: Optional[str] = None
    impedance_ohm: Optional[float] = None
    in_stock: Optional[bool] = None
    vat_included: Optional[bool] = None    # None: use the source's registry value


PARSE_ERROR_LIMIT = 2


@dataclass
class RunSummary:
    source_id: str
    pages: int = 0
    written: int = 0
    incomplete: int = 0
    skipped: List[Tuple[str, str]] = field(default_factory=list)   # (url, reason)
    review: int = 0
    parse_errors: int = 0                  # pages whose parser raised; also listed in skipped
    blocked_reason: Optional[str] = None
    outputs: List[str] = field(default_factory=list)

    def too_many_parse_errors(self) -> bool:
        """True when parser exceptions look like an adapter bug, not odd pages:
        more than PARSE_ERROR_LIMIT pages and more than 5 % of fetched pages."""
        return self.parse_errors > PARSE_ERROR_LIMIT and self.parse_errors > 0.05 * max(self.pages, 1)

    def to_dict(self) -> dict:
        return {"source_id": self.source_id, "pages": self.pages, "written": self.written,
                "incomplete": self.incomplete, "skipped": len(self.skipped), "parse_errors": self.parse_errors,
                "review": self.review, "blocked_reason": self.blocked_reason, "outputs": self.outputs}


class SourceAdapter(ABC):
    """Base for every source. Subclasses set source_id and kind and register
    with ``@horn_drivers.sources.register_adapter``."""

    source_id: ClassVar[str] = ""
    kind: ClassVar[str] = ""               # manufacturer | aggregator | retailer

    def __init__(self, fetcher: PoliteFetcher, entry=None):
        if fetcher.source_id != self.source_id:
            raise ValueError(f"fetcher is for {fetcher.source_id}, adapter is {self.source_id}")
        self.fetcher = fetcher
        self.entry = entry

    @abstractmethod
    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        """Yield product references, at most *limit*. May fetch index pages or
        sitemaps through ``self.fetcher.get``."""

    def fetch(self, ref: ProductRef) -> FetchResult:
        return self.fetcher.get(ref.url)

    def _pages(self, limit: Optional[int], summary: RunSummary) -> Iterator[Tuple[ProductRef, FetchResult]]:
        """Fetched pages for discovered refs; stops cleanly when the source is blocked."""
        try:
            for ref in self.discover(limit):
                try:
                    page = self.fetch(ref)
                except DisallowedByRobots as error:
                    summary.skipped.append((ref.url, str(error)))
                    continue
                except ResponseTooLarge as error:
                    summary.skipped.append((ref.url, f"response too large: {error}"))
                    continue
                summary.pages += 1
                if not page.ok:
                    summary.skipped.append((ref.url, f"HTTP {page.status}"))
                    continue
                yield ref, page
        except SourceBlocked as error:
            summary.blocked_reason = str(error)


class DriverSource(SourceAdapter):
    """A source of Thiele/Small records (manufacturer or aggregator)."""

    @abstractmethod
    def parse_driver(self, page: FetchResult, ref: ProductRef) -> DriverResult:
        """Build a version 2 record from one fetched page. No network access."""

    def run(self, db_dir: Path, incomplete_dir: Path, limit: Optional[int] = None) -> RunSummary:
        """Discover, fetch, parse, load-check and write records.

        A record the loader rejects, or that fails the schema check, goes to
        *incomplete_dir* with the reason, as the scraper does.
        """
        from horn_drivers.loader import _driver_from_dict
        summary = RunSummary(self.source_id)
        for ref, page in self._pages(limit, summary):
            try:
                result = self.parse_driver(page, ref)
            except (KeyError, ValueError, TypeError, AttributeError) as error:
                # A page the parser does not understand is skipped with the reason, not fatal.
                summary.skipped.append((ref.url, f"parse error: {type(error).__name__}: {error}"))
                summary.parse_errors += 1
                continue
            if result.record is None:
                summary.skipped.append((ref.url, result.reason or "no record"))
                continue
            reason = result.reason if result.incomplete else None
            if reason is None:
                errors, _ = check_structure(result.record)
                if errors:
                    reason = "schema: " + "; ".join(errors)
            if reason is None:
                try:
                    _driver_from_dict(result.record)
                except (KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
                    reason = f"loader rejects record: {error}"
            if reason is None:
                summary.outputs.append(str(write_record(db_dir, result.record)))
                summary.written += 1
            else:
                summary.outputs.append(str(write_incomplete(incomplete_dir, result.record, reason)))
                summary.incomplete += 1
        return summary


class DriverIndex:
    """Exact-match index from (canonical manufacturer, normalised model, impedance) to driver_id."""

    def __init__(self, records: Iterable[dict]):
        self.exact: Dict[Tuple[str, str, Optional[float]], str] = {}
        self.by_model: Dict[Tuple[str, str], List[Tuple[Optional[float], str]]] = {}
        for record in records:
            identity = record_identity(record)
            manufacturer = canonical_manufacturer(identity.get("manufacturer"))
            if not manufacturer or not identity.get("model"):
                continue
            impedance = (identity.get("variant") or {}).get("impedance_ohm")
            impedance = float(impedance) if impedance else None
            names = [identity["model"], *identity.get("aliases", [])]
            for name in names:
                model = normalise_model(name, manufacturer)
                self.exact[(manufacturer, model, impedance)] = record["driver_id"]
                self.by_model.setdefault((manufacturer, model), []).append((impedance, record["driver_id"]))

    def match(self, listing: Listing) -> Tuple[Optional[str], Optional[str], str]:
        """(driver_id, confidence, note). Confidence is "exact", "fuzzy" or None.

        Exact: canonical manufacturer, normalised model and impedance agree,
        and the title names both the manufacturer and the model. A listing
        without impedance whose model has a single variant is "fuzzy" and
        goes to review. Anything else is unmatched.
        """
        manufacturer = canonical_manufacturer(listing.manufacturer)
        if manufacturer is None:
            return None, None, f"unknown manufacturer {listing.manufacturer!r}"
        model = normalise_model(listing.model, manufacturer)
        title = normalise_model(listing.title)
        if model not in title:
            return None, None, "title does not name the model"
        from horn_drivers.identity import manufacturer_aliases
        if not any(normalise_model(alias) in title for alias in manufacturer_aliases(manufacturer)):
            return None, None, "title does not name the manufacturer"
        impedance = float(listing.impedance_ohm) if listing.impedance_ohm else None
        if impedance is not None and (manufacturer, model, impedance) in self.exact:
            return self.exact[(manufacturer, model, impedance)], "exact", ""
        variants = self.by_model.get((manufacturer, model), [])
        if len(variants) == 1 and (impedance is None or variants[0][0] is None):
            return variants[0][1], "fuzzy", "impedance not confirmed on both sides"
        if variants:
            return None, None, "several variants and no impedance match"
        return None, None, "model not in the database"


REVIEW_FIELDS = ("source_id", "url", "title", "manufacturer", "model", "impedance_ohm", "sku",
                 "price", "currency", "candidate_driver_id", "confidence", "note", "retrieved_on")


class PriceSource(SourceAdapter):
    """A retailer. Class attributes give the defaults its offers carry."""

    kind: ClassVar[str] = "retailer"
    region: ClassVar[str] = ""             # key of horn_drivers.prices.REGION_CURRENCY
    currency: ClassVar[str] = ""
    vat_included: ClassVar[Optional[bool]] = None

    @abstractmethod
    def parse_listing(self, page: FetchResult, ref: ProductRef) -> Optional[Listing]:
        """Facts from one product page, or None for a soft 404 / no single product."""

    def to_offer(self, listing: Listing, driver_id: str, confidence: str, retrieved_on: str) -> PriceOffer:
        return PriceOffer(driver_id=driver_id, source_id=self.source_id, url=listing.url, sku=listing.sku,
                          price=float(listing.price), currency=listing.currency,
                          vat_included=self.vat_included if listing.vat_included is None else listing.vat_included,
                          region=self.region, in_stock=listing.in_stock, retrieved_on=retrieved_on,
                          match_confidence=confidence)

    def run(self, index: DriverIndex, prices_dir: Path, review_csv: Path,
            limit: Optional[int] = None) -> RunSummary:
        """Discover, fetch, parse and match; write exact offers to the day's
        snapshot and everything else to the review CSV (fuzzy offers are
        written too, marked fuzzy, so the ranking ignores them until reviewed)."""
        summary = RunSummary(self.source_id)
        offers: List[PriceOffer] = []
        review_rows = []
        for ref, page in self._pages(limit, summary):
            listing = self.parse_listing(page, ref)
            if listing is None:
                summary.skipped.append((ref.url, "no single product on page"))
                continue
            driver_id, confidence, note = index.match(listing)
            if confidence is not None:
                offers.append(self.to_offer(listing, driver_id, confidence, page.retrieved_on))
            if confidence != "exact":
                review_rows.append({"source_id": self.source_id, "url": listing.url, "title": listing.title,
                                    "manufacturer": listing.manufacturer, "model": listing.model,
                                    "impedance_ohm": listing.impedance_ohm, "sku": listing.sku,
                                    "price": listing.price, "currency": listing.currency,
                                    "candidate_driver_id": driver_id, "confidence": confidence,
                                    "note": note, "retrieved_on": page.retrieved_on})
        if offers:
            for day in sorted({offer.retrieved_on for offer in offers}):
                path = append_offers(prices_dir, self.source_id, [o for o in offers if o.retrieved_on == day], day)
                summary.outputs.append(str(path))
            summary.written = sum(1 for o in offers if o.match_confidence == "exact")
        if review_rows:
            review_csv = Path(review_csv)
            review_csv.parent.mkdir(parents=True, exist_ok=True)
            new_file = not review_csv.exists()
            with review_csv.open("a", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=REVIEW_FIELDS)
                if new_file:
                    writer.writeheader()
                writer.writerows(review_rows)
            summary.outputs.append(str(review_csv))
        summary.review = len(review_rows)
        return summary
