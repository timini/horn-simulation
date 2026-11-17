"""Beyma (beyma.com): LF and midrange cone drivers from the product pages.

Product pages are WordPress HTML. Values sit in ``div.block-product-features
div.item`` as a title div and a description div with units ("21.6 N / A",
"163 um / N", "0.053 m2"). Sensitivity is stated as "1 W @ 1 m @ ZN".

Discovery reads the "Low and mid frequency loudspeakers" category, which
links every current model. The size range is checked on the product page.
Beyma also publishes the same values as an XLS catalogue; the pages are used
because they state units.
"""

import re
from typing import Iterator, Optional

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import DriverResult, DriverSource, ProductRef
from horn_drivers.sources.common import (REQUIRED_TS, SpecSheet, cell_text, diameter_label, number_range,
                                         sensitivity_1w, split_value_unit, stated_drive, states_one_watt, to_si)
from horn_drivers.sources.registry import register_adapter

BASE = "https://www.beyma.com"
CATEGORY_URL = f"{BASE}/en/products/c/low-mid-frequency/"
MIN_INCHES, MAX_INCHES = 5.0, 15.0
_PRODUCT = re.compile(r"""href=['"](https://www\.beyma\.com/en/products/c/low-mid-frequency/([A-Z0-9]+)/[a-z0-9\-]+/)['"]""")
_TITLE = re.compile(r"^\s*(?:LOUDSPEAKER|ALTAVOZ)\s+(.+?)\s+(\d+(?:[.,]\d+)?)\s*(?:Ohms?|Oh|OH|Ω)?\s*$", re.I)

LABELS = {
    "fs": "fs_hz", "re": "re_ohm", "qes": "qes", "qms": "qms", "qts": "qts", "sd": "sd_m2", "xmax": "xmax_m",
    "moving mass": "mms_kg", "mms": "mms_kg", "bl factor": "bl_tm", "le @1 khz": "le_h", "le": "le_h",
    "cms": "cms_m_per_n", "rms": "rms_kg_per_s", "nominal impedance": "nominal_impedance_ohm",
    "power capacity": "power_w", "program power": "peak_power_w",
}


def model_from_title(title: str):
    """("12MI100", 8.0) from "LOUDSPEAKER 12MI100 8 Ohm"; (title, None) when the pattern differs."""
    match = _TITLE.match(title or "")
    if not match:
        return (title or "").strip(), None
    return match.group(1).strip(), float(match.group(2).replace(",", "."))


@register_adapter
class BeymaSource(DriverSource):
    source_id = "beyma"
    kind = "manufacturer"

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        page = self.fetcher.get(CATEGORY_URL)
        if not page.ok:
            raise RuntimeError(f"Beyma category returned HTTP {page.status}")
        for count, ref in enumerate(self.refs_from_listing(page.text)):
            if limit is not None and count >= limit:
                return
            yield ref

    @staticmethod
    def refs_from_listing(html: str) -> Iterator[ProductRef]:
        seen = set()
        for url, sku in _PRODUCT.findall(html):
            if sku in seen:
                continue
            seen.add(sku)
            yield ProductRef(url=url, manufacturer="Beyma", sku=sku)

    def parse_driver(self, page: FetchResult, ref: ProductRef) -> DriverResult:
        return DriverResult(ref, *self.record_from_html(page.text, page.url, page.retrieved_on, ref.sku))

    def record_from_html(self, html: str, url: str, retrieved_on: str, sku: Optional[str] = None):
        from bs4 import BeautifulSoup
        soup = BeautifulSoup(html, "html.parser")
        title = soup.find("title")
        model, title_ohms = model_from_title(title.get_text(" ", strip=True) if title else "")
        if not model:
            return None, "page has no model title", False
        values = {}
        for item in soup.select("div.block-product-features div.item"):
            label, value = item.find("div", class_="title"), item.find("div", class_="description")
            if label is not None and value is not None:
                values.setdefault(label.get_text(" ", strip=True).lower(), cell_text(value))
        if not values:
            return None, "page has no specification table", False
        inches = None
        if "nominal diameter" in values:
            found = re.search(r"([\d.]+)\s*in\b", values["nominal diameter"])
            inches = float(found.group(1)) if found else None
        if inches is not None and not MIN_INCHES <= inches <= MAX_INCHES:
            return None, f"{inches:g} in is outside the {MIN_INCHES:g} to {MAX_INCHES:g} in pilot range", False
        sheet = SpecSheet(self.source_id, url, retrieved_on, "Beyma", model)
        for label, parameter in LABELS.items():
            if label not in values or parameter in sheet.parameters:
                continue
            # Beyma writes compound units with spaces ("N / A", "kg / s"); join them first.
            number, unit = split_value_unit(re.sub(r"\s*/\s*", "/", values[label]))
            if number is None or (number == 0 and parameter != "le_h"):
                continue
            sheet.put_stated(parameter, number, unit)
        if "nominal_impedance_ohm" not in sheet.parameters and title_ohms:
            sheet.put("nominal_impedance_ohm", title_ohms)
        if "sensitivity" in values:
            text = values["sensitivity"]
            spl, _ = split_value_unit(text)
            drive = stated_drive(text)
            impedance = sheet.parameters.get("nominal_impedance_ohm")
            if spl is not None and states_one_watt(text):
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
        record = sheet.record(driver_type="cone", nominal_diameter=diameter_label(inches),
                              aliases=[sku] if sku else ())
        if missing:
            return record, f"missing {', '.join(missing)} on the manufacturer page", True
        return record, None, False
