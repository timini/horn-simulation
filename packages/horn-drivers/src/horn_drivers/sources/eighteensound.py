"""Eighteen Sound (eighteensound.it): LF cone drivers from the product pages.

Product pages are static HTML. Each value sits in ``<li><span>label</span>
<sup class="note">n</sup><b>value unit</b></li>``; footnote *n* is the
n-th item of ``ul.catalogSpec`` and states the sensitivity drive voltage.

Discovery reads the LF driver listing, which links every current model and
impedance variant, and keeps the pilot size range. The datasheet PDF next to
each page is not read yet.
"""

import re
from typing import Iterator, Optional
from urllib.parse import unquote, urljoin

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import DriverResult, DriverSource, ProductRef
from horn_drivers.sources.common import (REQUIRED_TS, SpecSheet, cell_text, diameter_label, number_range, sensitivity_1w,
                                         split_value_unit, stated_drive, to_si)
from horn_drivers.sources.registry import register_adapter

BASE = "https://www.eighteensound.it"
LISTING_URL = f"{BASE}/en/products/lf-driver"
MIN_INCHES, MAX_INCHES = 5.0, 15.0
_LINK = re.compile(r'href="(/en/products/lf-driver/(\d+)-(\d+)/(\d+)/([^"/]+))"')

LABELS = {
    "resonance frequency": "fs_hz", "fs": "fs_hz", "re": "re_ohm", "qes": "qes", "qms": "qms", "qts": "qts",
    "sd": "sd_m2", "xmax": "xmax_m", "mms": "mms_kg", "bl": "bl_tm", "le": "le_h", "cms": "cms_m_per_n",
    "nominal impedance": "nominal_impedance_ohm", "nominal power handling": "power_w",
    "continuous power handling": "peak_power_w",
}


def _soup(html: str):
    from bs4 import BeautifulSoup
    return BeautifulSoup(html, "html.parser")


@register_adapter
class EighteenSoundSource(DriverSource):
    source_id = "eighteensound"
    kind = "manufacturer"

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        page = self.fetcher.get(LISTING_URL)
        if not page.ok:
            raise RuntimeError(f"18 Sound listing returned HTTP {page.status}")
        yield from self.refs_from_listing(page.text, limit)

    @staticmethod
    def refs_from_listing(html: str, limit: Optional[int] = None) -> Iterator[ProductRef]:
        seen = set()
        count = 0
        for path, whole, frac, ohms, slug in _LINK.findall(html):
            inches = float(f"{whole}.{frac}")
            key = path.lower().replace("%20", "").replace(" ", "")
            if key in seen or not MIN_INCHES <= inches <= MAX_INCHES:
                continue
            seen.add(key)
            if limit is not None and count >= limit:
                return
            count += 1
            yield ProductRef(url=urljoin(BASE, path), manufacturer="18 Sound", model=unquote(slug).strip(),
                             hints=(("inches", f"{inches:g}"), ("ohms", ohms)))

    def parse_driver(self, page: FetchResult, ref: ProductRef) -> DriverResult:
        return DriverResult(ref, *self.record_from_html(page.text, page.url, page.retrieved_on,
                                                        ref.hint("inches")))

    def record_from_html(self, html: str, url: str, retrieved_on: str, inches: Optional[str] = None):
        soup = _soup(html)
        heading = soup.find("h1")
        model = heading.get_text(" ", strip=True) if heading else ""
        if not model:
            return None, "page has no model heading", False
        notes = [li.get_text(" ", strip=True) for ul in soup.select("ul.catalogSpec") for li in ul.find_all("li")]
        values = {}
        for item in soup.select("div.productDetailsContent li"):
            label, value = item.find("span"), item.find("b")
            if label is None or value is None:
                continue
            note = item.find("sup", class_="note")
            values[label.get_text(" ", strip=True).lower()] = (cell_text(value),
                                                               note.get_text(strip=True) if note else "")
        sheet = SpecSheet(self.source_id, url, retrieved_on, "18Sound", model)
        for label, parameter in LABELS.items():
            if label not in values or parameter in sheet.parameters:
                continue
            number, unit = split_value_unit(values[label][0])
            if number is None or (number == 0 and parameter != "le_h"):
                continue
            sheet.put_stated(parameter, number, unit)
        if "sensitivity" in values:
            spl, _ = split_value_unit(values["sensitivity"][0])
            index = values["sensitivity"][1]
            note = notes[int(index) - 1] if index.isdigit() and 0 < int(index) <= len(notes) else ""
            drive = stated_drive(note) or next(filter(None, map(stated_drive, notes)), None)
            impedance = sheet.parameters.get("nominal_impedance_ohm")
            if spl is not None and drive and impedance and abs(drive[1] - impedance) > 0.01:
                sheet.notes.append(f"Sensitivity {spl:g} dB: the drive note names {drive[1]:g} ohm for a "
                                   f"{impedance:g} ohm model, so the 1 W figure is unknown; not stored.")
            elif spl is not None and drive:
                sheet.put("sensitivity_db_1w1m", sensitivity_1w(spl, drive[0], drive[1]), "derived")
            elif spl is not None:
                sheet.notes.append(f"Sensitivity {spl:g} dB stated without a drive level; not stored.")
        if "frequency range" in values:
            sheet.band(*number_range(values["frequency range"][0]))
        missing = sheet.missing(REQUIRED_TS)
        if len(missing) == len(REQUIRED_TS):
            return None, "no Thiele/Small parameters on page", False
        size = float(inches) if inches else None
        record = sheet.record(driver_type="cone", nominal_diameter=diameter_label(size))
        if missing:
            return record, f"missing {', '.join(missing)} on the manufacturer page", True
        return record, None, False
