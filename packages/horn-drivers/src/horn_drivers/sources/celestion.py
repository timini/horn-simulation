"""Celestion (celestion.com): pro LF cone drivers from the product pages.

Product pages are WordPress/WooCommerce HTML. Each value sits in a
``div.product-detail-spec-col-line`` with a label div and a value div; the
footnotes in ``#footnote`` state how power and sensitivity were measured.
robots.txt sets ``Crawl-delay: 10``, which the fetcher honours.

Discovery pages through the "LF Loudspeakers" category (12 products per
page). The size range is checked on the product page, since the listing does
not state it.
"""

import math
import re
from typing import Iterator, Optional

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import DriverResult, DriverSource, ProductRef
from horn_drivers.sources.common import (REQUIRED_TS, SpecSheet, diameter_label, number_range, sensitivity_1w,
                                         cell_text, split_value_unit, stated_drive, states_one_watt, to_si)
from horn_drivers.sources.registry import register_adapter

BASE = "https://celestion.com"
CATEGORY_URL = f"{BASE}/product-category/pro/lf-drivers/"
PER_PAGE = 12
MIN_INCHES, MAX_INCHES = 5.0, 15.0
_PRODUCT = re.compile(r'href="(https://celestion\.com/product/[a-z0-9\-]+/)"', re.I)
_COUNT = re.compile(r"of\s+(\d+)\s+results", re.I)

LABELS = {
    "fs": "fs_hz", "re": "re_ohm", "qes": "qes", "qms": "qms", "qts": "qts", "sd": "sd_m2", "xmax": "xmax_m",
    "mms": "mms_kg", "bl": "bl_tm", "bi": "bl_tm", "le": "le_h", "le (at 1khz)": "le_h", "cms": "cms_m_per_n",
    "rms": "rms_kg_per_s", "rated impedance": "nominal_impedance_ohm", "power rating": "power_w",
    "continuous power rating": "peak_power_w",
}


@register_adapter
class CelestionSource(DriverSource):
    source_id = "celestion"
    kind = "manufacturer"

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        first = self.fetcher.get(CATEGORY_URL)
        if not first.ok:
            raise RuntimeError(f"Celestion LF category returned HTTP {first.status}")
        match = _COUNT.search(first.text)
        pages = math.ceil(int(match.group(1)) / PER_PAGE) if match else 1
        texts = [first.text]
        for number in range(2, pages + 1):
            page = self.fetcher.get(f"{CATEGORY_URL}page/{number}/")
            if page.ok:
                texts.append(page.text)
        seen = set()
        for text in texts:
            for ref in self.refs_from_listing(text):
                if ref.url in seen:
                    continue
                if limit is not None and len(seen) >= limit:
                    return
                seen.add(ref.url)
                yield ref

    @staticmethod
    def refs_from_listing(html: str) -> Iterator[ProductRef]:
        for url in dict.fromkeys(_PRODUCT.findall(html)):
            yield ProductRef(url=url, manufacturer="Celestion", model=url.rstrip("/").rsplit("/", 1)[-1].upper())

    def parse_driver(self, page: FetchResult, ref: ProductRef) -> DriverResult:
        return DriverResult(ref, *self.record_from_html(page.text, page.url, page.retrieved_on))

    def record_from_html(self, html: str, url: str, retrieved_on: str):
        from bs4 import BeautifulSoup
        soup = BeautifulSoup(html, "html.parser")
        heading = soup.find(id="product-detail-title") or soup.find("h1")
        model = heading.get_text(" ", strip=True) if heading else ""
        if not model:
            return None, "page has no model heading", False
        values = {}
        for line in soup.select("div.product-detail-spec-col-line"):
            cells = line.find_all("div", recursive=False)
            if len(cells) >= 2:
                values.setdefault(cells[0].get_text(" ", strip=True).lower(), cell_text(cells[1]))
        footnotes = {}
        for paragraph in soup.select("#footnote p"):
            strong = paragraph.find("strong")
            if strong:
                footnotes[strong.get_text(strip=True).rstrip(":").lower()] = paragraph.get_text(" ", strip=True)
        inches = None
        if "nominal diameter" in values:
            found = re.search(r"([\d.]+)\s*in", values["nominal diameter"])
            inches = float(found.group(1)) if found else None
        if inches is not None and not MIN_INCHES <= inches <= MAX_INCHES:
            return None, f"{inches:g} in is outside the {MIN_INCHES:g} to {MAX_INCHES:g} in pilot range", False
        sheet = SpecSheet(self.source_id, url, retrieved_on, "Celestion", model)
        for label, parameter in LABELS.items():
            if label not in values or parameter in sheet.parameters:
                continue
            number, unit = split_value_unit(values[label])
            if number is None or (number == 0 and parameter != "le_h"):
                continue
            sheet.put_stated(parameter, number, unit)
        if "sensitivity" in values:
            spl, _ = split_value_unit(values["sensitivity"])
            note = footnotes.get("sensitivity", "")
            drive = stated_drive(note)
            impedance = sheet.parameters.get("nominal_impedance_ohm")
            if spl is not None and states_one_watt(note):
                sheet.put("sensitivity_db_1w1m", spl)
            elif spl is not None and drive and impedance and abs(drive[1] - impedance) <= 0.01:
                sheet.put("sensitivity_db_1w1m", sensitivity_1w(spl, drive[0], drive[1]), "derived")
            elif spl is not None:
                sheet.notes.append(f"Sensitivity {spl:g} dB stated without a usable drive level; not stored.")
        if "frequency range" in values:
            sheet.band(*number_range(values["frequency range"]))
        missing = sheet.missing(REQUIRED_TS)
        if len(missing) == len(REQUIRED_TS):
            return None, "no Thiele/Small parameters on page", False
        record = sheet.record(driver_type="cone", nominal_diameter=diameter_label(inches))
        if missing:
            return record, f"missing {', '.join(missing)} on the manufacturer page", True
        return record, None, False
