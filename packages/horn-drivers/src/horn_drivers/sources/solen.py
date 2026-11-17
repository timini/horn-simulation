"""Solen price adapter (CA, CAD, sales tax not included).

Solen's brand pages build their product grid in the browser from the shop's
own listing endpoint, which answers JSON when the request asks for it
(``Accept: application/json``); robots.txt allows everything. Each page of
that listing gives, per product: name (the model), a one-line summary that
starts with the brand and model and usually states the impedance, the
selling price, the price unit and stock. Only those facts are read.

A product priced per pair is not a per-driver price; it is left out of the
price file and listed for review instead.
"""

import json
from typing import ClassVar, Iterator, List, Optional, Tuple
from urllib.parse import quote

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import Listing, PriceSource, ProductRef
from horn_drivers.sources.registry import register_adapter
from horn_drivers.sources.retail import impedance_from_text, parse_money

BASE = "https://www.solen.ca/en"
JSON = "application/json"


@register_adapter
class Solen(PriceSource):
    source_id = "solen"
    region, currency, vat_included = "CA", "CAD", False
    # Solen brand slug -> manufacturer name. Of the pilot makers Solen stocks only Beyma.
    brands: ClassVar[Tuple[Tuple[str, str], ...]] = (("beyma", "Beyma"),)
    per_page: ClassVar[int] = 24

    def page_url(self, slug: str, page: int) -> str:
        return f"{BASE}/products?brand={quote(slug)}&perPage={self.per_page}&page={page}"

    def fetch(self, ref: ProductRef) -> FetchResult:
        return self.fetcher.get(ref.url, accept=JSON)

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        count = 0
        for slug, maker in self.brands:
            first = self.fetcher.get(self.page_url(slug, 1), accept=JSON)
            last = 1
            if first.ok:
                try:
                    last = int(json.loads(first.text)["products"]["last_page"])
                except (ValueError, KeyError, TypeError):
                    last = 1
            for page in range(1, last + 1):
                if limit is not None and count >= limit:
                    return
                count += 1
                yield ProductRef(self.page_url(slug, page), manufacturer=maker)

    def parse_listings(self, page: FetchResult, ref: ProductRef) -> List[Listing]:
        try:
            products = json.loads(page.text)["products"]["data"]
        except (ValueError, KeyError, TypeError):
            return []
        listings = []
        for product in products:
            listing = self.product_listing(product, ref.manufacturer or "")
            if listing is not None:
                listings.append(listing)
        return listings

    def product_listing(self, product: dict, maker: str) -> Optional[Listing]:
        try:
            model, slug = product["name"], product["slug"]
            price = parse_money((product.get("selling_price") or product["price"])["amount"])
        except (KeyError, TypeError):
            return None
        summary = product.get("short_description") or ""
        if price is None or not model or not summary:
            return None
        unit = (product.get("unit") or "each").lower()
        title = summary if unit == "each" else f"{summary} (priced per {unit})"
        in_stock = product.get("in_stock") if isinstance(product.get("in_stock"), bool) else None
        listing = Listing(url=f"{BASE}/products/{quote(slug)}", title=title, manufacturer=maker, model=model,
                          price=price, currency="CAD", sku=str(product.get("id")) if product.get("id") else None,
                          impedance_ohm=impedance_from_text(summary, allow_unlabelled=True), in_stock=in_stock)
        if unit != "each":
            # A pair price is not a driver price: send it to review with no impedance, so it cannot match exactly.
            listing = Listing(**{**listing.__dict__, "impedance_ohm": None, "model": f"{model} ({unit})"})
        return listing

    def parse_listing(self, page: FetchResult, ref: ProductRef) -> Optional[Listing]:
        listings = self.parse_listings(page, ref)
        return listings[0] if len(listings) == 1 else None
