"""Lautsprechershop.de price adapter (EU, EUR, VAT included).

The shop lists each maker's drivers on one static category page per
product group (``pa/18sound_woofer_en.htm``, ``pa/eminence_en.htm``). Every
product sits in its own block: a heading with the product name and an
anchor named after the order number, then a ``<preis artikel="...">``
element whose first amount is the gross price. Discovery reads the English
sitemap and keeps the category pages of the pilot makers. Only the name,
order number, gross price and stated impedance are read.
"""

import re
from typing import ClassVar, Iterator, List, Optional, Tuple

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import Listing, PriceSource, ProductRef
from horn_drivers.sources.registry import register_adapter
from horn_drivers.sources.retail import (brand_and_rest, impedance_from_text, parse_money, sitemap_locs,
                                         split_model_impedance, strip_tags)

SITEMAP = "https://www.lautsprechershop.de/sitemap_en.xml"
_HEADING = re.compile(r"<h2[^>]*>(.*?)<a\s+name=\"([^\"]+)\"[^>]*>", re.S | re.I)
_PRICE = re.compile(r"<preis\s+artikel=\"([^\"]+)\"[^>]*>.*?EUR\s*([\d.,]+)", re.S | re.I)


@register_adapter
class Lautsprechershop(PriceSource):
    source_id = "lautsprechershop"
    region, currency, vat_included = "EU", "EUR", True
    groups: ClassVar[Tuple[str, ...]] = ("18sound", "beyma", "celestion", "eminence", "bc")

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        pattern = re.compile(r"^https://www\.lautsprechershop\.de/(?:pa|chassis)/(?:%s)(?:_[a-z]+)*_en\.htm$"
                             % "|".join(map(re.escape, self.groups)))
        index = self.fetcher.get(SITEMAP)
        if not index.ok:
            return
        urls = sorted({u for u in sitemap_locs(index.content) if pattern.match(u)})
        for url in urls[:limit] if limit is not None else urls:
            yield ProductRef(url, hints=(("kind", "category"),))

    def parse_listings(self, page: FetchResult, ref: ProductRef) -> List[Listing]:
        text = page.content.decode("utf-8", errors="replace")
        if "�" in text:
            text = page.content.decode("latin-1")
        headings = list(_HEADING.finditer(text))
        listings = []
        for number, heading in enumerate(headings):
            end = headings[number + 1].start() if number + 1 < len(headings) else len(text)
            block = text[heading.start():end]
            listing = self.block_listing(page.url, strip_tags(heading.group(1)), heading.group(2), block)
            if listing is not None:
                listings.append(listing)
        return listings

    def block_listing(self, page_url: str, name: str, anchor: str, block: str) -> Optional[Listing]:
        price_match = _PRICE.search(block)
        if not price_match or price_match.group(1) != anchor:
            return None             # the block must price its own article
        price = parse_money(price_match.group(2))
        manufacturer, rest = brand_and_rest(name)
        if price is None or manufacturer is None or not rest:
            return None
        model, impedance = split_model_impedance(rest)
        if impedance is None:
            impedance = impedance_from_text(strip_tags(block))
        return Listing(url=f"{page_url}#{anchor}", title=name, manufacturer=manufacturer, model=model,
                       price=price, currency="EUR", sku=anchor, impedance_ohm=impedance)

    def parse_listing(self, page: FetchResult, ref: ProductRef) -> Optional[Listing]:
        listings = self.parse_listings(page, ref)
        return listings[0] if len(listings) == 1 else None
