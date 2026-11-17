"""Madisound price adapter (US, USD, sales tax not included).

Discovery reads the shop's sitemap and keeps product pages whose address
names a pilot maker (speaker kits excluded). A product page carries a
schema.org Product in JSON-LD with brand, name and an Offer (price,
currency, availability); only those facts are read. The model is the first
word after the brand in the product name, and the impedance comes from the
name ("AlphaLite-6B 16 ohm").

Madisound stocks few pro drivers of the pilot makers (on 2026-11-17 the
sitemap listed one Eminence driver and no B&C, 18 Sound, Celestion cone or
Beyma product pages), so this source adds little to the pilot.
"""

import re
from typing import ClassVar, Iterator, Optional, Tuple

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import Listing, PriceSource, ProductRef
from horn_drivers.sources.registry import register_adapter
from horn_drivers.sources.retail import first_token, impedance_from_text, jsonld_products, parse_money, sitemap_locs, split_brand

SITEMAP = "https://www.madisoundspeakerstore.com/sitemap.xml"


@register_adapter
class Madisound(PriceSource):
    source_id = "madisound"
    region, currency, vat_included = "US", "USD", False
    brand_words: ClassVar[Tuple[str, ...]] = ("eminence", "beyma", "celestion", "18-sound", "eighteen-sound", "b-c")

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        pattern = re.compile(r"^https://www\.madisoundspeakerstore\.com/[^/]+/(?:%s)-[^/]+/$"
                             % "|".join(map(re.escape, self.brand_words)))
        index = self.fetcher.get(SITEMAP)
        if not index.ok:
            return
        urls = sorted({u for u in sitemap_locs(index.content) if pattern.match(u) and "kit" not in u})
        for url in urls[:limit] if limit is not None else urls:
            yield ProductRef(url)

    def parse_listing(self, page: FetchResult, ref: ProductRef) -> Optional[Listing]:
        products = list(jsonld_products(page.text))
        if len(products) != 1:
            return None
        product = products[0]
        offer = product.get("offers")
        if isinstance(offer, list):
            offer = offer[0] if len(offer) == 1 else None
        brand = product.get("brand")
        brand = brand.get("name") if isinstance(brand, dict) else (brand or product.get("manufacturer"))
        name = product.get("name") or ""
        if not isinstance(offer, dict) or not brand or not name:
            return None
        price = parse_money(offer.get("price"))
        rest = split_brand(name, brand)
        if price is None or not rest or (offer.get("priceCurrency") or "USD") != "USD":
            return None
        availability = str(offer.get("availability") or "")
        in_stock = None if not availability else availability.endswith("InStock")
        return Listing(url=page.url, title=name, manufacturer=brand, model=first_token(rest), price=price,
                       currency="USD", sku=product.get("sku"),
                       impedance_ohm=impedance_from_text(name, allow_unlabelled=True), in_stock=in_stock)
