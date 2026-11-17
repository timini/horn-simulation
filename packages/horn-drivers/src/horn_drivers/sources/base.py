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
import re
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


@dataclass
class RunSummary:
    source_id: str
    pages: int = 0
    written: int = 0
    incomplete: int = 0
    skipped: List[Tuple[str, str]] = field(default_factory=list)   # (url, reason)
    review: int = 0
    blocked_reason: Optional[str] = None
    outputs: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"source_id": self.source_id, "pages": self.pages, "written": self.written,
                "incomplete": self.incomplete, "skipped": len(self.skipped), "review": self.review,
                "blocked_reason": self.blocked_reason, "outputs": self.outputs}


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
            result = self.parse_driver(page, ref)
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

    def model_from_title(self, manufacturer: str, title: str) -> Optional[str]:
        """The longest known model of *manufacturer* named in *title*, or None.

        The title is split into letter-and-digit tokens and a model counts as
        named when a run of consecutive tokens joins to its normalised form.
        A run followed by a short revision token ("2", "II", "MK2") is
        ambiguous (Beta-12A versus Beta-12A-2) and gives None unless the
        longer name is itself a known model.
        """
        canonical = canonical_manufacturer(manufacturer)
        known = {model for (maker, model) in self.by_model if maker == canonical}
        if not known:
            return None
        tokens = re.findall(r"[A-Z0-9]+", title.upper().replace("&", " "))
        best: Optional[Tuple[int, str]] = None
        for start in range(len(tokens)):
            joined = ""
            for end in range(start, min(len(tokens), start + 6)):
                joined += tokens[end]
                if joined in known:
                    following = tokens[end + 1] if end + 1 < len(tokens) else ""
                    after = tokens[end + 2] if end + 2 < len(tokens) else ""
                    if (_REVISION.fullmatch(following) and after not in _UNIT_WORDS
                            and (joined + following) not in known):
                        continue
                    if best is None or len(joined) > best[0]:
                        best = (len(joined), joined)
        return best[1] if best else None

    def match(self, listing: Listing) -> Tuple[Optional[str], Optional[str], str]:
        """(driver_id, confidence, note). Confidence is "exact", "fuzzy" or None.

        Exact: canonical manufacturer, normalised model and impedance agree,
        and the title names both the manufacturer and the model. A listing
        without impedance whose model has a single variant is "fuzzy" and
        goes to review. A part or kit for a driver (repair kit, diaphragm)
        is never matched. When the listing's model is not in the index, a
        known model named in the title is used (``model_from_title``).
        Anything else is unmatched.
        """
        manufacturer = canonical_manufacturer(listing.manufacturer)
        if manufacturer is None:
            return None, None, f"unknown manufacturer {listing.manufacturer!r}"
        if _ACCESSORY.search(listing.title) or _ACCESSORY.search(listing.model):
            return None, None, "a part or kit for a driver, not the driver"
        model = normalise_model(listing.model, manufacturer)
        title = normalise_model(listing.title)
        from_title = ""
        if (manufacturer, model) not in self.by_model:
            found = self.model_from_title(manufacturer, listing.title)
            if found:
                model, from_title = found, f"model {found} read from the title; "
        if not model or model not in title:
            return None, None, "title does not name the model"
        from horn_drivers.identity import manufacturer_aliases
        if not any(normalise_model(alias) in title for alias in manufacturer_aliases(manufacturer)):
            return None, None, "title does not name the manufacturer"
        impedance = float(listing.impedance_ohm) if listing.impedance_ohm else None
        if impedance is not None and (manufacturer, model, impedance) in self.exact:
            return self.exact[(manufacturer, model, impedance)], "exact", from_title.rstrip("; ")
        variants = self.by_model.get((manufacturer, model), [])
        if len(variants) == 1 and (impedance is None or variants[0][0] is None):
            return variants[0][1], "fuzzy", from_title + "impedance not confirmed on both sides"
        if variants:
            known = sorted({str(v[0]) for v in variants})
            return None, None, from_title + f"impedance {impedance} not among database variants {', '.join(known)}"
        return None, None, "model not in the database"

    def knows_model(self, listing: Listing) -> bool:
        """True when the listing's model (or one named in its title) is in the index."""
        manufacturer = canonical_manufacturer(listing.manufacturer)
        if manufacturer is None:
            return False
        if (manufacturer, normalise_model(listing.model, manufacturer)) in self.by_model:
            return True
        return self.model_from_title(manufacturer, listing.title) is not None


_REVISION = re.compile(r"\d{1,2}|II|III|IV|MK\d*|V\d")
# A listing for a part of a driver (recone or repair kit, diaphragm, gasket)
# names the driver model too; it is never a price for the driver itself.
_ACCESSORY = re.compile(r"\b(?:repair|recone|reconing|re-cone|diaphragms?|membranes?|ersatzmembran\w*|"
                        r"reparatur\w*|spare|kit|gasket|grilles?|replacement\s+(?:cone|part))\b", re.IGNORECASE)
_UNIT_WORDS = frozenset({"OHM", "OHMS", "INCH", "IN", "W", "WATT", "WATTS", "MM", "CM", "X"})

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

    def parse_listings(self, page: FetchResult, ref: ProductRef) -> List[Listing]:
        """Every product on a fetched page. A category page may hold many,
        each read from its own product block; the default reads one product
        page with ``parse_listing``. No network access."""
        listing = self.parse_listing(page, ref)
        return [listing] if listing is not None else []

    def detail_ref(self, listing: Listing) -> Optional[ProductRef]:
        """A product page to fetch when a category listing lacks the impedance.

        ``run`` fetches it only when the model is in the driver index, and
        uses the product page's own ``parse_listing`` result instead. The
        default is None (no follow-up fetch)."""
        return None

    def _detail(self, detail: ProductRef, listing: Listing, retrieved_on: str,
                summary: RunSummary) -> Tuple[Listing, str]:
        """The product page's listing in place of *listing*, when it can be read."""
        try:
            page = self.fetch(detail)
        except (DisallowedByRobots, ResponseTooLarge) as error:
            summary.skipped.append((detail.url, str(error)))
            return listing, retrieved_on
        summary.pages += 1
        if not page.ok:
            summary.skipped.append((detail.url, f"HTTP {page.status}"))
            return listing, retrieved_on
        found = self.parse_listing(page, detail)
        if found is None:
            summary.skipped.append((detail.url, "product page shows no single product"))
            return listing, retrieved_on
        return found, page.retrieved_on

    def to_offer(self, listing: Listing, driver_id: str, confidence: str, retrieved_on: str) -> PriceOffer:
        return PriceOffer(driver_id=driver_id, source_id=self.source_id, url=listing.url, sku=listing.sku,
                          price=float(listing.price), currency=listing.currency,
                          vat_included=self.vat_included if listing.vat_included is None else listing.vat_included,
                          region=self.region, in_stock=listing.in_stock, retrieved_on=retrieved_on,
                          match_confidence=confidence)

    def run(self, index: DriverIndex, prices_dir: Path, review_csv: Path,
            limit: Optional[int] = None) -> RunSummary:
        """Discover, fetch, parse and match; write exact offers to the day's
        snapshot and everything else (fuzzy or unmatched) to the review CSV
        only. A reviewer can add a checked offer by hand as "manual"."""
        summary = RunSummary(self.source_id)
        offers: List[PriceOffer] = []
        review_rows = []
        seen = set()
        try:
            for ref, page in self._pages(limit, summary):
                listings = self.parse_listings(page, ref)
                if not listings:
                    summary.skipped.append((ref.url, "no product on page"))
                    continue
                for listing in listings:
                    retrieved_on = page.retrieved_on
                    if listing.impedance_ohm is None and index.knows_model(listing):
                        detail = self.detail_ref(listing)
                        if detail is not None:
                            listing, retrieved_on = self._detail(detail, listing, retrieved_on, summary)
                    key = (listing.url, listing.sku, listing.model, listing.impedance_ohm)
                    if key in seen:
                        continue
                    seen.add(key)
                    driver_id, confidence, note = index.match(listing)
                    if confidence == "exact":
                        offers.append(self.to_offer(listing, driver_id, confidence, retrieved_on))
                    else:
                        review_rows.append({"source_id": self.source_id, "url": listing.url, "title": listing.title,
                                            "manufacturer": listing.manufacturer, "model": listing.model,
                                            "impedance_ohm": listing.impedance_ohm, "sku": listing.sku,
                                            "price": listing.price, "currency": listing.currency,
                                            "candidate_driver_id": driver_id, "confidence": confidence,
                                            "note": note, "retrieved_on": retrieved_on})
        except SourceBlocked as error:
            summary.blocked_reason = str(error)
        if offers:
            for day in sorted({offer.retrieved_on for offer in offers}):
                path = append_offers(prices_dir, self.source_id, [o for o in offers if o.retrieved_on == day], day)
                summary.outputs.append(str(path))
            summary.written = len(offers)
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

