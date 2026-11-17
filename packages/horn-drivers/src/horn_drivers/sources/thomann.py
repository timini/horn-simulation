"""Thomann UK (GBP) and Thomann DE (EUR) price adapters.

Discovery reads the shop's sitemap index and the first sitemap file that
lists brand-and-size loudspeaker categories (for example
``b-c_loudspeaker_12_inch.html`` on the UK site,
``b-c_12_zoll_lautsprecher.html`` on the German one). A category page
carries a JSON block with one entry per product: manufacturer, model,
article number, price with currency and VAT flag, stock state and a short
spec list. The adapter reads only those facts, plus the impedance from the
model name or the spec list.

When a category entry gives no impedance and the model is in the driver
index, the product page is fetched; its tracking block names brand, model
and price, and its meta description usually gives the impedance. A product
page without that block (Thomann answers unknown product URLs with HTTP 200
and the title "Unknown Item") is a soft 404 and yields nothing.
"""

import re
from typing import ClassVar, Iterator, List, Optional, Tuple
from urllib.parse import urljoin

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import Listing, PriceSource, ProductRef
from horn_drivers.sources.registry import register_adapter
from horn_drivers.sources.retail import (html_title, impedance_from_text, json_after, meta_content, sitemap_locs,
                                         split_model_impedance)

# Thomann brand slugs for the pilot manufacturers.
PILOT_BRANDS = ("b-c", "beyma", "celestion", "eighteensound", "eminence")
ARTICLES_MARKER = '"articleListsSettings":{"articles":'
PAGING_MARKER = '"pagingSettings":'


class ThomannBase(PriceSource):
    base: ClassVar[str] = ""               # shop root, for example https://www.thomann.co.uk/
    category_pattern: ClassVar[str] = ""   # regex with {brands}
    brands: ClassVar[Tuple[str, ...]] = PILOT_BRANDS

    # -- discovery -----------------------------------------------------------

    def category_urls(self) -> List[str]:
        pattern = re.compile(self.category_pattern.format(brands="|".join(map(re.escape, self.brands))))
        index = self.fetcher.get(urljoin(self.base, "/sitemap.xml"))
        if not index.ok:
            return []
        for sitemap in sitemap_locs(index.content):
            page = self.fetcher.get(sitemap)
            if not page.ok:
                continue
            found = sorted({url for url in sitemap_locs(page.content) if pattern.match(url)})
            if found:
                return found          # categories sit together in one sitemap file
        return []

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        count = 0
        for url in self.category_urls():
            pages = [url]
            first = self.fetcher.get(url)
            last = self.last_page(first.text) if first.ok else 1
            pages += [f"{url}?pg={n}" for n in range(2, last + 1)]
            for page_url in pages:
                if limit is not None and count >= limit:
                    return
                count += 1
                yield ProductRef(page_url, hints=(("kind", "category"),))

    @staticmethod
    def last_page(text: str) -> int:
        paging = json_after(text, PAGING_MARKER)
        try:
            return max(1, int(paging.get("lastPage") or 1)) if isinstance(paging, dict) else 1
        except (TypeError, ValueError):
            return 1

    # -- parsing ---------------------------------------------------------------

    def product_url(self, relative: str) -> str:
        return urljoin(self.base, relative.split("?", 1)[0])

    def parse_listings(self, page: FetchResult, ref: ProductRef) -> List[Listing]:
        articles = json_after(page.text, ARTICLES_MARKER)
        if not isinstance(articles, list):
            listing = self.parse_listing(page, ref)
            return [listing] if listing else []
        listings = []
        for article in articles:
            listing = self.article_listing(article)
            if listing is not None:
                listings.append(listing)
        return listings

    def article_listing(self, article: dict) -> Optional[Listing]:
        try:
            primary = article["price"]["primary"]
            price = float(primary["rawPrice"])
            currency = primary["currency"]["key"]
            manufacturer, raw_model = article["manufacturer"], article["model"]
            relative = article["relativeLink"]
        except (KeyError, TypeError, ValueError):
            return None
        if price <= 0 or not raw_model:
            return None
        model, impedance = split_model_impedance(raw_model)
        if impedance is None:
            impedance = impedance_from_text(((article.get("texts") or {}).get("description")) or "")
        token = (article.get("availability") or {}).get("token")
        return Listing(url=self.product_url(relative), title=f"{manufacturer} {raw_model}",
                       manufacturer=manufacturer, model=model, price=price, currency=currency,
                       sku=str(article.get("number") or "") or None, impedance_ohm=impedance,
                       in_stock=(token == "in-stock") if token else None,
                       vat_included=primary.get("includesVat"))

    def parse_listing(self, page: FetchResult, ref: ProductRef) -> Optional[Listing]:
        """One product page; None for a soft 404 or a page without the product block."""
        text = page.text
        title = html_title(text)
        if "unknown item" in title.lower():
            return None
        view = text.find('"event":"view_item"')
        items = json_after(text[view:], '"ecommerce":{"items":') if view >= 0 else None
        if not isinstance(items, list) or len(items) != 1:
            return None
        item = items[0]
        try:
            price = float(item["price"])
            manufacturer, raw_model, currency = item["item_brand"], item["item_name"], item["currency"]
        except (KeyError, TypeError, ValueError):
            return None
        if price <= 0:
            return None
        model, impedance = split_model_impedance(raw_model)
        description = meta_content(text, "description") or ""
        if impedance is None:
            impedance = impedance_from_text(description, allow_unlabelled=True)
        page_title = meta_content(text, "og:title") or title
        return Listing(url=page.url, title=page_title, manufacturer=manufacturer, model=model, price=price,
                       currency=currency, sku=str(item.get("item_id") or "") or None, impedance_ohm=impedance)

    def detail_ref(self, listing: Listing) -> Optional[ProductRef]:
        return ProductRef(listing.url, listing.manufacturer, listing.model, listing.sku,
                          hints=(("kind", "product"),))


@register_adapter
class ThomannUK(ThomannBase):
    source_id = "thomann-uk"
    region, currency, vat_included = "UK", "GBP", True
    base = "https://www.thomann.co.uk/"
    category_pattern = r"^https://www\.thomann\.co\.uk/({brands})_loudspeaker_[0-9_]+_inch\.html$"


@register_adapter
class ThomannDE(ThomannBase):
    source_id = "thomann-de"
    region, currency, vat_included = "EU", "EUR", True
    # English pages on thomann.de live under /intl/, which robots.txt disallows; the German pages are used.
    base = "https://www.thomann.de/de/"
    category_pattern = r"^https://www\.thomann\.de/de/({brands})_[0-9-]+_zoll_lautsprecher\.html$"
