"""Parsing helpers shared by the retailer adapters.

Everything here works on text already fetched; nothing touches the network.
The helpers read facts only: a price, a currency, an impedance, a model
name. They never keep descriptions.
"""

import gzip
import html
import json
import re
import xml.etree.ElementTree as ET
from typing import Iterator, List, Optional, Tuple

from horn_drivers.identity import canonical_manufacturer, manufacturer_aliases

_OHM = r"(\d+(?:[.,]\d+)?)\s*-?\s*(?:ohms?|Ω|Ohm)\b"
_LABELLED_OHM = re.compile(r"(?:impedance|impedanz|imp\.)(?:\s*[a-z]{1,2}\b)?\s*[:=]?\s*" + _OHM, re.IGNORECASE)
_ANY_OHM = re.compile(r"(?<![\w.,])" + _OHM, re.IGNORECASE)
_MODEL_OHM_SUFFIX = re.compile(r"[\s,/-]*\(?\s*-?\s*(\d+(?:[.,]\d+)?)\s*(?:ohms?|Ω)\s*\)?\s*$", re.IGNORECASE)


def _number(text: str) -> float:
    return float(text.replace(",", "."))


def _unique(values: List[float]) -> Optional[float]:
    distinct = sorted(set(values))
    return distinct[0] if len(distinct) == 1 else None


def impedance_from_text(text: str, allow_unlabelled: bool = False) -> Optional[float]:
    """Nominal impedance stated in *text*, or None when absent or ambiguous.

    A value labelled "Impedance"/"Impedanz" wins. Without a label, a single
    distinct "<n> Ohm" value is used only when *allow_unlabelled* is set,
    because spec lists also quote the DC resistance Re in ohms.
    """
    if not text:
        return None
    labelled = [_number(m.group(1)) for m in _LABELLED_OHM.finditer(text)]
    if labelled:
        return _unique(labelled)
    if allow_unlabelled:
        return _unique([_number(m.group(1)) for m in _ANY_OHM.finditer(text)])
    return None


def split_model_impedance(model: str) -> Tuple[str, Optional[float]]:
    """("12NW76", 8.0) from "12NW76 8 Ohms"; the model unchanged when no suffix."""
    match = _MODEL_OHM_SUFFIX.search(model)
    if not match or match.start() == 0:
        return model.strip(), None
    return model[:match.start()].strip(" ,-/"), _number(match.group(1))


def parse_money(text: str) -> Optional[float]:
    """A price from "205,30", "1.234,50", "1,234.50", "CA$760.79" or "142.99"."""
    if text is None:
        return None
    digits = re.sub(r"[^\d.,]", "", str(text))
    if not digits:
        return None
    if "," in digits and "." in digits:
        if digits.rfind(",") > digits.rfind("."):
            digits = digits.replace(".", "").replace(",", ".")
        else:
            digits = digits.replace(",", "")
    elif "," in digits:
        head, _, tail = digits.rpartition(",")
        digits = f"{head.replace(',', '')}.{tail}" if len(tail) in (1, 2) else digits.replace(",", "")
    try:
        value = float(digits)
    except ValueError:
        return None
    return value if value > 0 else None


def split_brand(name: str, manufacturer: str) -> Optional[str]:
    """The text after a leading brand name, or None when *name* does not start with it."""
    canonical = canonical_manufacturer(manufacturer)
    if canonical is None:
        return None
    for alias in sorted(manufacturer_aliases(canonical), key=len, reverse=True):
        pattern = r"^\s*" + r"[\s\-]*".join(re.escape(part) for part in alias.split()) + r"(?:®|\b)[\s,:-]*"
        match = re.match(pattern, name, flags=re.IGNORECASE)
        if match:
            return name[match.end():].strip()
    return None


def first_token(text: str) -> str:
    """The first whitespace-delimited token, without trailing punctuation."""
    parts = text.split()
    return parts[0].strip(",;:") if parts else ""


def sitemap_locs(content: bytes) -> List[str]:
    """<loc> values from a sitemap or sitemap index, gzip or plain."""
    if content[:2] == b"\x1f\x8b":
        content = gzip.decompress(content)
    try:
        root = ET.fromstring(content)
    except ET.ParseError:
        return [html.unescape(m) for m in re.findall(r"<loc>\s*([^<\s]+)\s*</loc>", content.decode("utf-8", "replace"))]
    return [el.text.strip() for el in root.iter() if el.tag.endswith("loc") and el.text]


def json_after(text: str, marker: str):
    """Decode the JSON value that starts right after *marker* in *text*, or None."""
    start = text.find(marker)
    if start < 0:
        return None
    try:
        value, _ = json.JSONDecoder().raw_decode(text, start + len(marker))
    except ValueError:
        return None
    return value


def jsonld_products(text: str) -> Iterator[dict]:
    """Every schema.org Product object in the page's JSON-LD blocks."""
    for block in re.findall(r"<script[^>]*application/ld\+json[^>]*>(.*?)</script>", text, flags=re.S | re.I):
        try:
            data = json.loads(block, strict=False)
        except ValueError:
            continue
        stack = [data]
        while stack:
            item = stack.pop()
            if isinstance(item, list):
                stack.extend(item)
            elif isinstance(item, dict):
                kind = item.get("@type")
                if kind == "Product" or (isinstance(kind, list) and "Product" in kind):
                    yield item
                stack.extend(v for v in item.values() if isinstance(v, (dict, list)))


def html_title(text: str) -> str:
    match = re.search(r"<title[^>]*>(.*?)</title>", text, flags=re.S | re.I)
    return html.unescape(re.sub(r"\s+", " ", match.group(1))).strip() if match else ""


def meta_content(text: str, name: str) -> Optional[str]:
    """content of <meta name=... > or <meta property=...>."""
    for pattern in (rf'<meta[^>]+(?:name|property)="{re.escape(name)}"[^>]*content="([^"]*)"',
                    rf'<meta[^>]+content="([^"]*)"[^>]*(?:name|property)="{re.escape(name)}"'):
        match = re.search(pattern, text, flags=re.I)
        if match:
            return html.unescape(match.group(1))
    return None


def strip_tags(fragment: str) -> str:
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", fragment))).strip()


def brand_and_rest(name: str) -> Tuple[Optional[str], str]:
    """(canonical manufacturer, remaining text) for a name that starts with a known brand."""
    from horn_drivers.identity import MANUFACTURERS
    best: Tuple[int, Optional[str], str] = (0, None, name)
    for canonical in MANUFACTURERS:
        rest = split_brand(name, canonical)
        if rest is not None and len(name) - len(rest) > best[0]:
            best = (len(name) - len(rest), canonical, rest)
    return best[1], best[2]
