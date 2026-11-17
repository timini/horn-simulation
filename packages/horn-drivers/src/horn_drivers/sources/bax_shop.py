"""Bax-shop UK price adapter (UK, GBP, VAT included).

Discovery reads the pro audio sitemap and keeps loudspeaker product pages
(``/<n>-inch-speakers/<brand>-<model>``) of the pilot makers. A product
page names its main product in a ``data-product-visited`` attribute (id,
SKU, brand, name, stock count) and prices it in its tracking block
``window.gtm.products.item<id>``; other products on the page (related
items) are ignored. The model is the first word after the brand in the
product name; the matcher falls back to models named anywhere in the name.
"""

import html
import json
import re
from typing import ClassVar, Iterator, Optional, Tuple

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import Listing, PriceSource, ProductRef
from horn_drivers.sources.registry import register_adapter
from horn_drivers.sources.retail import (first_token, html_title, impedance_from_text, json_after, parse_money,
                                         sitemap_locs, split_brand)

SITEMAP = "https://sitemap.bax-shop.co.uk/en_gb/sitemap-pro-audio.xml"


@register_adapter
class BaxShopUK(PriceSource):
    source_id = "bax-shop-uk"
    region, currency, vat_included = "UK", "GBP", True
    brand_slugs: ClassVar[Tuple[str, ...]] = ("18-sound", "b-c", "beyma", "celestion", "eminence")

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        pattern = re.compile(r"^https://www\.bax-shop\.co\.uk/[0-9]+-inch-speakers/(?:%s)-[a-z0-9-]+$"
                             % "|".join(map(re.escape, self.brand_slugs)))
        index = self.fetcher.get(SITEMAP)
        if not index.ok:
            return
        urls = sorted({u for u in sitemap_locs(index.content) if pattern.match(u)})
        for url in urls[:limit] if limit is not None else urls:
            yield ProductRef(url)

    def parse_listing(self, page: FetchResult, ref: ProductRef) -> Optional[Listing]:
        text = page.text
        match = re.search(r'data-product-visited="([^"]*)"', text)
        if not match:
            return None
        try:
            visited = json.loads(html.unescape(match.group(1)))
            product_id, brand, name = str(visited["id"]), visited["brand"], visited["name"]
        except (ValueError, KeyError, TypeError):
            return None
        item = json_after(text, f"window.gtm.products.item{product_id} = ")
        price = parse_money(item.get("price")) if isinstance(item, dict) else None
        if price is None:
            price = parse_money(((visited.get("price") or {}).get("actual")))
        rest = split_brand(name, brand) if brand else None
        if price is None or not rest:
            return None
        amount = (visited.get("stock") or {}).get("amount")
        try:
            in_stock = int(amount) > 0 if amount is not None else None
        except (TypeError, ValueError):
            in_stock = None
        title = html_title(text) or name
        return Listing(url=page.url, title=title if name.lower() in title.lower() else name, manufacturer=brand,
                       model=first_token(rest), price=price, currency="GBP", sku=visited.get("sku"),
                       impedance_ohm=impedance_from_text(name, allow_unlabelled=True), in_stock=in_stock)
