"""B&C Speakers (bcspeakers.com): LF cone drivers from the product pages.

The site is a Remix app. Every page embeds its loader data as JSON
(``window.__remixContext``), so values are read from structured fields with
their units instead of from the rendered text.

Discovery reads the LF driver listing (one page holds every current model
and impedance variant) and yields one product page per variant within the
pilot size range. Archived models are not listed there.
"""

import json
import re
from typing import Iterator, Optional, Tuple
from urllib.parse import quote

from horn_drivers.fetch import FetchResult
from horn_drivers.sources.base import DriverResult, DriverSource, ProductRef
from horn_drivers.sources.common import (REQUIRED_TS, SpecSheet, diameter_label, sensitivity_1w, stated_drive,
                                         to_si)
from horn_drivers.sources.registry import register_adapter

BASE = "https://www.bcspeakers.com"
LISTING_URL = f"{BASE}/en/products/lf-driver"
MIN_INCHES, MAX_INCHES = 5.0, 15.0

_CONTEXT = re.compile(r"window\.__remixContext\s*=\s*(\{.*?\});\s*</script>", re.S)

# page field -> (parameter, preferred unit key for dict values)
PARAMS = {
    "fs": ("fs_hz", None), "re": ("re_ohm", None), "qes": ("qes", None), "qms": ("qms", None),
    "qts": ("qts", None), "sd": ("sd_m2", "cm2"), "xmax": ("xmax_m", "mm"), "mms": ("mms_kg", None),
    "bl": ("bl_tm", None), "le": ("le_h", None), "cms": ("cms_m_per_n", None),
}


def remix_data(html: str) -> dict:
    """The route loader data (everything but ``root``) from a page."""
    match = _CONTEXT.search(html)
    if not match:
        raise ValueError("page has no __remixContext data")
    loader = json.loads(match.group(1))["state"]["loaderData"]
    routes = [value for key, value in loader.items() if key != "root"]
    if len(routes) != 1 or not isinstance(routes[0], dict):
        raise ValueError("unexpected loader data layout")
    return routes[0]


def field_value(spec: Optional[dict], unit_key: Optional[str] = None) -> Tuple[Optional[float], str]:
    """(value, unit) of one spec field; dict values are read at *unit_key*."""
    if not isinstance(spec, dict):
        return None, ""
    value, units = spec.get("value"), spec.get("units") or {}
    if isinstance(value, dict):
        key = unit_key or next(iter(value), None)
        value = value.get(key)
        unit = units.get(key, key or "")
    else:
        unit = units.get("default", "") if isinstance(units, dict) else ""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None, unit
    return float(value), unit


def product_url(inches: float, ohms: float, name: str) -> str:
    return f"{LISTING_URL}/{inches:g}/{ohms:g}/{quote(name, safe='')}"


@register_adapter
class BCSource(DriverSource):
    source_id = "bc"
    kind = "manufacturer"

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        page = self.fetcher.get(LISTING_URL)
        if not page.ok:
            raise RuntimeError(f"B&C listing returned HTTP {page.status}")
        yield from self.refs_from_listing(page.text, limit)

    @staticmethod
    def refs_from_listing(html: str, limit: Optional[int] = None) -> Iterator[ProductRef]:
        products = remix_data(html).get("products") or {}
        count = 0
        for size_key in sorted(products, key=float):
            if not MIN_INCHES <= float(size_key) <= MAX_INCHES:
                continue
            for family in products[size_key]:
                for variant in family.get("variants", []):
                    if variant.get("archived") or variant.get("noLongerInProduction"):
                        continue
                    ohms = variant.get("nominalImpedance")
                    if not variant.get("name") or not ohms:
                        continue
                    if limit is not None and count >= limit:
                        return
                    count += 1
                    yield ProductRef(url=product_url(float(size_key), float(ohms), variant["name"]),
                                     manufacturer="B&C", model=variant["name"],
                                     hints=(("inches", f"{float(size_key):g}"), ("ohms", f"{float(ohms):g}")))

    def parse_driver(self, page: FetchResult, ref: ProductRef) -> DriverResult:
        try:
            product = remix_data(page.text)["product"]
        except (ValueError, KeyError) as error:
            return DriverResult(ref, None, f"no product data: {error}")
        return DriverResult(ref, *self.record_from_product(product, page.url, page.retrieved_on))

    def record_from_product(self, product: dict, url: str, retrieved_on: str):
        """(record, reason, incomplete) for one product variant."""
        name = (product.get("name") or "").strip()
        if not name:
            return None, "product has no name", False
        sheet = SpecSheet(self.source_id, url, retrieved_on, "BC", name)
        params = product.get("params") or {}
        for page_field, (parameter, unit_key) in PARAMS.items():
            value, unit = field_value(params.get(page_field), unit_key)
            if value is None or (value == 0 and parameter not in ("le_h",)):
                continue
            sheet.put(parameter, to_si(value, unit))
        specs = product.get("lfSpecs") or {}
        impedance, _ = field_value(specs.get("nominalImpedance"))
        sheet.put("nominal_impedance_ohm", impedance or product.get("nominalImpedance"))
        power, _ = field_value(specs.get("nominalPowerHandling"))
        sheet.put("power_w", power)
        program, _ = field_value(specs.get("continuousPowerHandling"))
        sheet.put("peak_power_w", program)
        spl, _ = field_value(specs.get("sensitivity"))
        if spl is not None:
            drive = stated_drive((specs.get("sensitivity") or {}).get("notes", ""))
            if drive and impedance and abs(drive[1] - impedance) > 0.01:
                sheet.notes.append(f"Sensitivity {spl:g} dB: the drive note names {drive[1]:g} ohm for a "
                                   f"{impedance:g} ohm model, so the 1 W figure is unknown; not stored.")
            elif drive:
                sheet.put("sensitivity_db_1w1m", sensitivity_1w(spl, drive[0], drive[1]), "derived")
            else:
                sheet.notes.append(f"Sensitivity {spl:g} dB stated without a drive level; not stored.")
        band = (specs.get("frequencyRange") or {}).get("value")
        if isinstance(band, dict):
            sheet.band(band.get("min"), band.get("max"))
        inches = product.get("inches")
        missing = sheet.missing(REQUIRED_TS)
        if len(missing) == len(REQUIRED_TS):
            return None, "no Thiele/Small parameters on page", False
        record = sheet.record(driver_type="cone", nominal_diameter=diameter_label(inches))
        if missing:
            return record, f"missing {', '.join(missing)} on the manufacturer page", True
        return record, None, False
