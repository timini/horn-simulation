"""Helpers shared by the manufacturer adapters: numbers, units, record building.

Manufacturer pages state values with display units ("0.37 mH", "346.0 cm2",
"14.5 Txm"). ``to_si`` turns a value and its unit into the SI value stored
under ``parameters``. ``SpecSheet`` collects the values read from one page
and builds a version 2 record with one provenance entry per field.

Sensitivity is stored as dB SPL for 1 W at 1 m. Makers state it at a drive
voltage (2.83 V for 8 ohm, 4 V for 16 ohm, 2 V for 4 ohm). ``sensitivity_1w``
converts a stated voltage into the 1 W figure for the nominal impedance.
"""

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from horn_drivers.records import make_record, provenance_entry

_NUMBER = re.compile(r"[-+]?\d+(?:[.,]\d+)?")

# unit -> (parameter family, factor to SI). Units are compared after _unit_key().
_UNIT_FACTORS: Dict[str, float] = {
    "hz": 1.0, "khz": 1000.0,
    "ohm": 1.0, "ohms": 1.0, "ω": 1.0, "Ω": 1.0,
    "tm": 1.0, "txm": 1.0, "t.m": 1.0, "n/a": 1.0,
    "cm2": 1e-4, "cm²": 1e-4, "m2": 1.0, "m²": 1.0, "mm2": 1e-6,
    "g": 1e-3, "gr": 1e-3, "kg": 1.0,
    "mh": 1e-3, "uh": 1e-6, "µh": 1e-6, "μh": 1e-6, "h": 1.0,
    "mm": 1e-3, "m": 1.0, "cm": 1e-2,
    "w": 1.0, "db": 1.0, "": 1.0,
    "mm/n": 1e-3, "m/n": 1.0, "um/n": 1e-6, "µm/n": 1e-6, "μm/n": 1e-6,
    "kg/s": 1.0, "ns/m": 1.0, "n.s/m": 1.0,
    "dm3": 1e-3, "dm³": 1e-3, "l": 1e-3, "litres": 1e-3, "liters": 1e-3,
}


# Units each stored parameter may be stated in. A value whose unit belongs to
# another quantity ("Rms 1.804 mm/N") is a page error and is not stored.
_ALLOWED_UNITS: Dict[str, Tuple[str, ...]] = {
    "fs_hz": ("hz", "khz"),
    "re_ohm": ("ohm", "ohms", "ω", "Ω"),
    "nominal_impedance_ohm": ("ohm", "ohms", "ω", "Ω"),
    "bl_tm": ("tm", "txm", "t.m", "n/a"),
    "sd_m2": ("cm2", "cm²", "m2", "m²", "mm2"),
    "mms_kg": ("g", "gr", "kg"),
    "le_h": ("mh", "uh", "µh", "μh", "h"),
    "xmax_m": ("mm", "m", "cm"),
    "power_w": ("w",),
    "peak_power_w": ("w",),
    "cms_m_per_n": ("mm/n", "m/n", "um/n", "µm/n", "μm/n"),
    "rms_kg_per_s": ("kg/s", "ns/m", "n.s/m"),
    "vas_m3": ("dm3", "dm³", "l", "litres", "liters"),
}

# A supplied Cms or Rms must agree with Fs, Mms and Qms within this fraction;
# the simulation uses a supplied value as given.
CONSISTENCY_TOLERANCE = 0.10


def unit_allowed(name: str, unit: str) -> bool:
    """True when *unit* is a unit *name* may be stated in (unknown names pass)."""
    allowed = _ALLOWED_UNITS.get(name)
    return allowed is None or _unit_key(unit) in allowed


def _unit_key(unit: str) -> str:
    return unit.strip().lower().replace("·", ".").replace(" ", "")


def first_number(text) -> Optional[float]:
    """First number in *text* ("350 W" gives 350.0), or None. Accepts a number."""
    if text is None:
        return None
    if isinstance(text, (int, float)):
        return float(text)
    match = _NUMBER.search(str(text))
    return float(match.group(0).replace(",", ".")) if match else None


def number_range(text) -> Tuple[Optional[float], Optional[float]]:
    """(low, high) from "65 - 5000 Hz" or "65 Hz - 3000 Hz"; (None, None) when absent."""
    if text is None:
        return None, None
    numbers = [float(n.replace(",", ".")) for n in _NUMBER.findall(str(text).replace("\u2013", "-"))]
    numbers = [abs(n) for n in numbers]
    if len(numbers) >= 2:
        return numbers[0], numbers[1]
    return None, None


def to_si(value: Optional[float], unit: str) -> Optional[float]:
    """*value* in *unit* converted to SI; raises KeyError for an unknown unit."""
    if value is None:
        return None
    key = _unit_key(unit)
    if key not in _UNIT_FACTORS:
        raise KeyError(f"unknown unit {unit!r}")
    return value * _UNIT_FACTORS[key]


def split_value_unit(text: str) -> Tuple[Optional[float], str]:
    """("0.37", "mH") style split of a display value such as "0.37 mH (x in)".

    The imperial repeat is dropped, whether in brackets or after a slash
    followed by a number ("44.51g/1.57oz"); "mm/N" stays a unit. The
    unit is the first word after the number.
    """
    if text is None:
        return None, ""
    text = re.split(r"\s*/\s*(?=\d)", re.sub(r"\(.*?\)", "", str(text)))[0].strip()
    match = re.match(r"\s*([-+]?\d+(?:[.,]\d+)?)\s*([^\s(]*)", text)
    if not match:
        return None, ""
    return float(match.group(1).replace(",", ".")), match.group(2)


_VOLTAGE = re.compile(r"(\d+(?:[.,]\d+)?)\s*V(?:rms)?\b[^.]*?\bfor\s+(\d+(?:[.,]\d+)?)\s*ohm", re.IGNORECASE)
_ONE_WATT = re.compile(r"(?<![\d.])1\s*W\b", re.IGNORECASE)


def states_one_watt(note: str) -> bool:
    """True when a sensitivity note says the figure is for 1 W (e.g. "at 1W, 1m")."""
    return bool(note and _ONE_WATT.search(note))


def stated_drive(note: str) -> Optional[Tuple[float, float]]:
    """(volts, ohms) from a note such as "Applied RMS Voltage is set to 4V for 16 ohm"."""
    if not note:
        return None
    match = _VOLTAGE.search(note)
    if not match:
        return None
    return float(match.group(1).replace(",", ".")), float(match.group(2).replace(",", "."))


def sensitivity_1w(spl_db: float, volts: float, impedance_ohm: float) -> float:
    """SPL at 1 W into the nominal impedance, from SPL measured at *volts*."""
    return round(spl_db - 20.0 * math.log10(volts / math.sqrt(impedance_ohm)), 1)


@dataclass
class SpecSheet:
    """Values read from one page, each with its provenance method."""

    source_id: str
    url: str
    retrieved_on: str
    manufacturer: str
    model: str
    parameters: Dict[str, float] = field(default_factory=dict)
    methods: Dict[str, str] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)
    usable_f_low_hz: Optional[float] = None
    usable_f_high_hz: Optional[float] = None

    def put(self, name: str, value: Optional[float], method: str = "page") -> None:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return
        self.parameters[name] = float(f"{float(value):.6g}")
        self.methods[name] = method

    def put_stated(self, name: str, number: float, unit: str) -> None:
        """Store a page value in SI. A dimensioned value shown without a unit is
        left out with a note rather than guessed."""
        if not unit and name not in UNITLESS:
            self.notes.append(f"{name} shown as {number:g} without a unit; not stored.")
            return
        if unit and not unit_allowed(name, unit):
            self.notes.append(f"{name} shown as {number:g} {unit}, a unit of another quantity; not stored.")
            return
        self.put(name, to_si(number, unit))

    def drop_inconsistent(self) -> None:
        """Drop a supplied Cms or Rms that disagrees with Fs, Mms and Qms.

        Cms must give Fs = 1/(2 pi sqrt(Mms Cms)) and Rms must equal
        2 pi Fs Mms / Qms, each within CONSISTENCY_TOLERANCE. A value that
        does not is left out with a note; the loader derives it instead.
        """
        p = self.parameters
        fs, mms, qms = p.get("fs_hz"), p.get("mms_kg"), p.get("qms")
        cms = p.get("cms_m_per_n")
        if cms and fs and mms:
            implied = 1.0 / (2 * math.pi * math.sqrt(mms * cms))
            if abs(implied - fs) / fs > CONSISTENCY_TOLERANCE:
                self._drop("cms_m_per_n", f"page Cms {cms:.4g} m/N gives Fs {implied:.1f} Hz against the "
                                          f"stated {fs:g} Hz")
        rms = p.get("rms_kg_per_s")
        if rms and fs and mms and qms:
            expected = 2 * math.pi * fs * mms / qms
            if abs(rms - expected) / expected > CONSISTENCY_TOLERANCE:
                self._drop("rms_kg_per_s", f"page Rms {rms:.4g} kg/s against {expected:.3g} kg/s from "
                                           f"Fs, Mms and Qms")

    def _drop(self, name: str, why: str) -> None:
        self.parameters.pop(name, None)
        self.methods.pop(name, None)
        self.notes.append(f"{why}; differs by more than {100 * CONSISTENCY_TOLERANCE:g} %, so not stored.")

    def band(self, low: Optional[float], high: Optional[float]) -> None:
        if low and high and 0 < low < high:
            self.usable_f_low_hz, self.usable_f_high_hz = float(low), float(high)

    def missing(self, names) -> List[str]:
        return [name for name in names if self.parameters.get(name) is None]

    def record(self, *, driver_type: str = "cone", nominal_diameter: Optional[str] = None,
               aliases=(), coil: Optional[str] = None) -> dict:
        self.drop_inconsistent()
        provenance = {name: provenance_entry(self.source_id, self.url, self.retrieved_on, method)
                      for name, method in self.methods.items()}
        if self.usable_f_low_hz is not None:
            for name in ("usable_f_low_hz", "usable_f_high_hz"):
                provenance[name] = provenance_entry(self.source_id, self.url, self.retrieved_on, "page")
        return make_record(manufacturer=self.manufacturer, model=self.model, parameters=self.parameters,
                           field_provenance=provenance, driver_type=driver_type,
                           impedance_ohm=self.parameters.get("nominal_impedance_ohm"), coil=coil,
                           aliases=aliases, nominal_diameter=nominal_diameter,
                           usable_f_low_hz=self.usable_f_low_hz, usable_f_high_hz=self.usable_f_high_hz,
                           id_includes_impedance=True, notes=self.notes)


UNITLESS = ("qms", "qes", "qts")
REQUIRED_TS = ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg", "le_h")


def diameter_label(inches: Optional[float]) -> Optional[str]:
    return f"{inches:g}in" if inches else None


def cell_text(tag) -> str:
    """Text of an HTML cell with <sup>/<sub> joined to the unit ("cm<sup>2</sup>" reads "cm2")."""
    for small in tag.find_all(["sup", "sub"]):
        small.replace_with(small.get_text())
    return " ".join(tag.get_text("").split())
