"""Dayton Audio list prices (US, USD, MSRP).

Dayton Audio's own product pages state the MSRP; the products are sold
through Parts Express, whose pages need JavaScript and are not fetched. The
MSRP is a manufacturer list price, so it is a ceiling rather than a shop
offer. Discovery walks the pro audio drivers category (``?pagenum=N``).
A product page gives the model as the first word of its heading and the
impedance in the heading text.
"""

import re
from typing import ClassVar, Iterator, Optional, Tuple
from urllib.parse import urljoin

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import Listing, PriceSource, ProductRef
from horn_drivers.sources.registry import register_adapter
from horn_drivers.sources.retail import first_token, html_title, impedance_from_text, parse_money, strip_tags

BASE = "https://www.daytonaudio.com"
_PRODUCT_LINK = re.compile(r'href="(/product/\d+/[a-z0-9-]+)"')
_PAGENUM = re.compile(r"pagenum=(\d+)")
_HEADING = re.compile(r'<h1[^>]*itemprop="name"[^>]*>(.*?)</h1>', re.S | re.I)
_MSRP = re.compile(r'class="variant-price"[^>]*>\s*MSRP:\s*\$\s*([\d,]+(?:\.\d{2})?)', re.I)
_AVAILABILITY = re.compile(r'itemprop="availability"\s+href="https?://schema\.org/(\w+)"', re.I)


@register_adapter
class DaytonMSRP(PriceSource):
    source_id = "dayton-msrp"
    region, currency, vat_included = "US", "USD", False
    price_kind = "list"
    categories: ClassVar[Tuple[str, ...]] = ("/category/161/pro-audio-drivers",)

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        seen = set()
        for category in self.categories:
            first = self.fetcher.get(urljoin(BASE, category))
            if not first.ok:
                continue
            pages = sorted({int(n) for n in _PAGENUM.findall(first.text)} | {1})
            for number in pages:
                page = first if number == 1 else self.fetcher.get(urljoin(BASE, f"{category}?pagenum={number}"))
                if not page.ok:
                    continue
                for link in _PRODUCT_LINK.findall(page.text):
                    url = urljoin(BASE, link)
                    if url in seen:
                        continue
                    if limit is not None and len(seen) >= limit:
                        return
                    seen.add(url)
                    yield ProductRef(url)

    def parse_listing(self, page: FetchResult, ref: ProductRef) -> Optional[Listing]:
        text = page.text
        heading, msrp = _HEADING.search(text), _MSRP.search(text)
        if not heading or not msrp:
            return None
        name = strip_tags(heading.group(1))
        price = parse_money(msrp.group(1))
        if price is None or not name:
            return None
        availability = _AVAILABILITY.search(text)
        in_stock = None if not availability else availability.group(1).lower() == "instock"
        title = html_title(text) or f"Dayton Audio {name}"
        return Listing(url=page.url, title=title, manufacturer="Dayton Audio", model=first_token(name),
                       price=price, currency="USD", impedance_ohm=impedance_from_text(name, allow_unlabelled=True),
                       in_stock=in_stock)
