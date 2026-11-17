"""Canonical driver identity: manufacturer names, model normalisation, match keys.

The canonical manufacturer name is the folder name used under ``data/drivers``
(for example ``BC`` for B&C Speakers), so new records land next to the old
ones. Aliases cover the spellings that retailers and datasheets use.

A match key is (canonical manufacturer, normalised model, impedance). Source
adapters and the retailer matcher use it to join a listing to a driver
record. Anything that does not match exactly goes to review, not to a guess.
"""

import re
from typing import Dict, List, Optional, Tuple

# canonical name -> (display name, aliases). Aliases are compared after
# _fold(), so case, spaces and punctuation do not matter.
MANUFACTURERS: Dict[str, Tuple[str, Tuple[str, ...]]] = {
    "18Sound": ("18 Sound", ("Eighteen Sound", "Eighteensound", "18 Sound", "18-Sound")),
    "BC": ("B&C Speakers", ("B&C", "B & C", "B and C", "BC Speakers", "B&C Speakers")),
    "Beyma": ("Beyma", ()),
    "BMS": ("BMS", ("BMS Speakers",)),
    "Celestion": ("Celestion", ()),
    "Ciare": ("Ciare", ()),
    "Dayton": ("Dayton Audio", ("Dayton Audio",)),
    "Eminence": ("Eminence", ("Eminence Speaker",)),
    "Faital": ("Faital Pro", ("FaitalPRO", "Faital Pro", "Faital")),
    "Fostex": ("Fostex", ()),
    "JBL": ("JBL Professional", ("JBL Pro", "JBL Professional")),
    "LaVoce": ("LaVoce", ("La Voce",)),
    "Monacor": ("Monacor", ()),
    "Peerless": ("Peerless by Tymphany", ("Tymphany", "Peerless by Tymphany")),
    "PRV": ("PRV Audio", ("PRV Audio",)),
    "RCF": ("RCF", ()),
    "SB Acoustics": ("SB Acoustics", ("SBAcoustics", "SB Audience")),
    "Scan-Speak": ("Scan-Speak", ("ScanSpeak", "Scan Speak")),
    "SEAS": ("SEAS", ()),
    "Selenium": ("Selenium", ("JBL Selenium",)),
    "SICA": ("SICA", ("Sica Altoparlanti",)),
    "Visaton": ("Visaton", ()),
}


def _fold(text: str) -> str:
    """Lower case with everything but letters and digits removed; & reads as 'and'."""
    return re.sub(r"[^a-z0-9]", "", text.lower().replace("&", "and"))


_ALIAS_INDEX: Dict[str, str] = {}
for _canonical, (_display, _aliases) in MANUFACTURERS.items():
    for _name in (_canonical, _display, *_aliases):
        _ALIAS_INDEX[_fold(_name)] = _canonical


def canonical_manufacturer(name: Optional[str]) -> Optional[str]:
    """Canonical manufacturer for *name*, or None when it is not on the list."""
    if not name:
        return None
    return _ALIAS_INDEX.get(_fold(name))


def is_canonical_manufacturer(name: Optional[str]) -> bool:
    return name in MANUFACTURERS


def manufacturer_aliases(canonical: str) -> List[str]:
    """Every known spelling of a canonical manufacturer, canonical first."""
    display, aliases = MANUFACTURERS[canonical]
    return list(dict.fromkeys((canonical, display, *aliases)))


def normalise_model(model: str, manufacturer: Optional[str] = None) -> str:
    """Model name reduced for matching: upper case, no spaces or punctuation.

    A leading manufacturer alias is stripped, so "Eminence Beta-12A" and
    "BETA12A" give the same result for manufacturer "Eminence".
    """
    text = model.strip()
    canonical = canonical_manufacturer(manufacturer) if manufacturer else None
    if canonical:
        for alias in sorted(manufacturer_aliases(canonical), key=len, reverse=True):
            pattern = r"^\s*" + r"[\s\-]*".join(re.escape(part) for part in alias.split()) + r"[\s\-]+"
            text = re.sub(pattern, "", text, flags=re.IGNORECASE)
    return re.sub(r"[^A-Z0-9]", "", text.upper())


def match_key(manufacturer: str, model: str, impedance_ohm: Optional[float]) -> Tuple[str, str, Optional[float]]:
    """(canonical manufacturer, normalised model, impedance) for exact matching.

    Raises ValueError when the manufacturer is not on the canonical list, so
    an unknown brand is never silently matched under a wrong name.
    """
    canonical = canonical_manufacturer(manufacturer)
    if canonical is None:
        raise ValueError(f"Unknown manufacturer {manufacturer!r}; add it to horn_drivers.identity.MANUFACTURERS")
    return canonical, normalise_model(model, canonical), (float(impedance_ohm) if impedance_ohm else None)


def record_identity(record: dict) -> dict:
    """The ``identity`` block of a record, built from v1 fields when absent."""
    identity = record.get("identity")
    if isinstance(identity, dict):
        return identity
    params = record.get("parameters") or {}
    return {
        "manufacturer": canonical_manufacturer(record.get("manufacturer")) or record.get("manufacturer"),
        "model": record.get("model_name"),
        "variant": {"impedance_ohm": params.get("nominal_impedance_ohm")},
        "aliases": [],
    }
