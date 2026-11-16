"""Driver database loader for horn simulation.

Loads loudspeaker driver Thiele-Small parameters from a database,
converting stored units to SI and deriving missing parameters.

Supports two database layouts:
  - **Directory** (v3): ``data/drivers/{Manufacturer}/{driver-id}.json``
    One JSON file per driver.  This is the preferred format.
  - **Single-file** (v1/v2): a single ``drivers.json`` containing all
    drivers in an array or dict-of-dicts.  Retained for backward
    compatibility.

All public functions accept a ``db_path`` that may point to either a
directory or a JSON file.
"""

import json
import math
import warnings
from pathlib import Path
from typing import List, Optional

from horn_core.parameters import DriverParameters


# Relative mismatch between a supplied Qts and the Qts implied by supplied
# Qms and Qes (Qms*Qes/(Qms+Qes)) above which a record is inconsistent.
# Datasheets round Qs to two or three figures, which explains a few percent.
Q_IDENTITY_TOLERANCE = 0.10
# Qms above this is physically unusual and usually a symptom of deriving it
# from nearly equal Qes and Qts.
QMS_PLAUSIBLE_MAX = 30.0

# Values far outside anything a loudspeaker driver can have. They almost
# always mean the number was entered in the wrong unit.
_GROSS_LIMITS = {
    "sd_m2": (0.5, "Sd is in m^2; use sd_cm2 for square centimetres"),
    "mms_kg": (2.0, "Mms is in kg; use mms_g for grams"),
    "le_h": (0.05, "Le is in H; use le_mh for millihenries"),
    "xmax_m": (0.1, "Xmax is in m; use xmax_mm for millimetres"),
    "exit_area_m2": (0.1, "exit area is in m^2; use exit_area_cm2 or exit_diameter_mm"),
}

RHO_AIR = 1.204   # kg/m^3, 20 C
C_AIR = 343.0     # m/s


def q_identity_error(qms: float, qes: float, qts: float) -> float:
    """Relative error of the supplied Qts against Qms*Qes/(Qms+Qes)."""
    implied = qms * qes / (qms + qes)
    return abs(implied - qts) / qts


def _first(params: dict, *names):
    """First alias that is present (not None), so explicit zeros still reach validation."""
    for name in names:
        value = params.get(name)
        if value is not None:
            return value
    return None


def _scaled(value, factor):
    return value * factor if value is not None else None


def _driver_from_dict(d: dict) -> DriverParameters:
    """Convert a driver dict (from JSON) into a DriverParameters instance.

    Handles both v1 (flat) and v2 (nested ``parameters``) formats, and
    converts common datasheet units on the fly: mH, cm^2, mm, g, litres.
    A record that gives Vas instead of Mms gets Mms from Vas, Sd and Fs.
    """
    # v2 format nests T-S params under "parameters"
    params = d.get("parameters", d)
    driver_id = d.get("driver_id", "unknown")

    # Unit conversions — source may use convenience units
    # Only missing values fall back to aliases. Explicit zero limits must
    # reach validation, and zero inductance is a valid limiting case.
    le_h = params.get("le_h")
    if le_h is None:
        le_h = _mh_to_h(params.get("le_mh"))
    sd_m2 = _first(params, "sd_m2", "sd_sq_meters")
    if sd_m2 is None:
        sd_m2 = _cm2_to_m2(_first(params, "sd_cm2", "sd_sq_cm"))
    xmax_m = params.get("xmax_m")
    if xmax_m is None:
        xmax_m = _mm_to_m(params.get("xmax_mm"))
    exit_area_m2 = params.get("exit_area_m2")
    if exit_area_m2 is None:
        exit_area_m2 = _cm2_to_m2(params.get("exit_area_cm2"))
    if exit_area_m2 is None and params.get("exit_diameter_mm") is not None:
        exit_area_m2 = math.pi * (params["exit_diameter_mm"] * 1e-3 / 2.0) ** 2
    front_chamber_volume_m3 = params.get("front_chamber_volume_m3")
    if front_chamber_volume_m3 is None:
        front_chamber_volume_m3 = _cm3_to_m3(params.get("front_chamber_volume_cm3"))
    front_chamber_area_m2 = params.get("front_chamber_area_m2")
    if front_chamber_area_m2 is None:
        front_chamber_area_m2 = _cm2_to_m2(params.get("front_chamber_area_cm2"))
    rear_chamber_volume_m3 = params.get("rear_chamber_volume_m3")
    if rear_chamber_volume_m3 is None:
        rear_chamber_volume_m3 = _cm3_to_m3(params.get("rear_chamber_volume_cm3"))
    re_ohm = params.get("re_ohm")
    if re_ohm is None:
        re_ohm = params.get("re_ohms", 0.0)
    fs_hz = _first(params, "fs_hz", "fs")
    if fs_hz is None:
        raise KeyError("fs_hz")
    mms_kg = params.get("mms_kg")
    if mms_kg is None:
        mms_kg = _scaled(params.get("mms_g"), 1e-3)
    vas_m3 = _scaled(_first(params, "vas_l", "vas_liters", "vas_litres"), 1e-3)
    if mms_kg is None and vas_m3 is not None and sd_m2 and vas_m3 > 0:
        # Cms = Vas / (rho c^2 Sd^2); Mms = 1 / ((2 pi fs)^2 Cms)
        cms = vas_m3 / (RHO_AIR * C_AIR ** 2 * sd_m2 ** 2)
        mms_kg = 1.0 / ((2.0 * math.pi * fs_hz) ** 2 * cms)
    power_w = _first(params, "power_w", "power_aes_w", "power_aes_watts")

    for name, value in (("sd_m2", sd_m2), ("mms_kg", mms_kg), ("le_h", le_h),
                        ("xmax_m", xmax_m), ("exit_area_m2", exit_area_m2)):
        limit, hint = _GROSS_LIMITS[name]
        if isinstance(value, (int, float)) and value > limit:
            raise ValueError(f"{driver_id}: {name}={value} is not a plausible driver value ({hint})")

    qms, qes, qts = params.get("qms"), params.get("qes"), params.get("qts")
    if all(isinstance(q, (int, float)) and q > 0 for q in (qms, qes, qts)):
        error = q_identity_error(qms, qes, qts)
        if error > Q_IDENTITY_TOLERANCE:
            raise ValueError(
                f"{driver_id}: Q factors are inconsistent: Qms={qms}, Qes={qes} imply "
                f"Qts={qms * qes / (qms + qes):.3f}, record has {qts} ({error:.0%} apart)")

    driver = DriverParameters(
        driver_id=driver_id,
        manufacturer=d.get("manufacturer", ""),
        model_name=d.get("model_name", ""),
        fs_hz=fs_hz,
        re_ohm=re_ohm,
        bl_tm=params.get("bl_tm", 0.0),
        sd_m2=sd_m2 if sd_m2 is not None else 0.0,
        mms_kg=mms_kg if mms_kg is not None else 0.0,
        le_h=le_h if le_h is not None else 0.0,
        qms=qms,
        qes=qes,
        qts=qts,
        exit_area_m2=exit_area_m2,
        driver_type=d.get("driver_type"),
        nominal_diameter=d.get("nominal_diameter"),
        xmax_m=xmax_m,
        nominal_impedance_ohm=params.get("nominal_impedance_ohm"),
        power_w=power_w,
        peak_power_w=params.get("peak_power_w"),
        usable_f_low_hz=d.get("usable_f_low_hz"),
        usable_f_high_hz=d.get("usable_f_high_hz"),
        parameter_source=d.get("parameter_source"),
        interface_model=d.get("interface_model"),
        mmd_kg=params.get("mmd_kg"),
        rear_load_mass_kg=params.get("rear_load_mass_kg"),
        cms_m_per_n=params.get("cms_m_per_n"),
        rms_kg_per_s=params.get("rms_kg_per_s"),
        front_chamber_volume_m3=front_chamber_volume_m3,
        front_chamber_area_m2=front_chamber_area_m2,
        rear_chamber_volume_m3=rear_chamber_volume_m3,
    )
    if qms is None and driver.qms is not None and driver.qms > QMS_PLAUSIBLE_MAX:
        warnings.warn(
            f"{driver_id}: derived Qms={driver.qms:.1f} is implausibly high; "
            "supply Qms or Rms from the datasheet", RuntimeWarning, stacklevel=3)
    return driver


def _load_from_directory(db_dir: Path) -> List[dict]:
    """Load all driver JSON files from a directory tree.

    Expects ``db_dir/{Manufacturer}/{driver-id}.json``.
    """
    drivers: List[dict] = []
    for manufacturer_dir in sorted(db_dir.iterdir()):
        if not manufacturer_dir.is_dir():
            continue
        for driver_file in sorted(manufacturer_dir.glob("*.json")):
            try:
                d = json.loads(driver_file.read_text())
                drivers.append(d)
            except (json.JSONDecodeError, KeyError) as e:
                print(f"WARN: skipping {driver_file}: {e}")
    return drivers


def load_drivers_raw(db_path: str) -> List[dict]:
    """Load all drivers as raw dicts (before conversion to DriverParameters).

    Accepts a directory (v3) or a single JSON file (v1/v2).
    """
    p = Path(db_path)

    # v3: directory of per-driver JSON files
    if p.is_dir():
        return _load_from_directory(p)

    # v1/v2: single JSON file
    raw = json.loads(p.read_text())
    if isinstance(raw, dict) and "drivers" in raw:
        return raw["drivers"]
    # v1 dict-of-dicts
    drivers = []
    for driver_id, data in raw.items():
        data.setdefault("driver_id", driver_id)
        drivers.append(data)
    return drivers


def load_drivers(db_path: str) -> List[DriverParameters]:
    """Load valid drivers, with explicit warnings for rejected database records.

    Supports:
      - v3 directory: ``db_path/{Manufacturer}/{driver-id}.json``
      - v2 single file: ``{ "schema_version": 2, "drivers": [...] }``
      - v1 single file: ``{ "driver_id": { ... }, ... }``
    """
    drivers = []
    for record in load_drivers_raw(db_path):
        try:
            drivers.append(_driver_from_dict(record))
        except (KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
            warnings.warn(f"Rejected driver {record.get('driver_id', 'unknown')}: {error}", RuntimeWarning, stacklevel=2)
    return drivers


def load_driver(db_path: str, driver_id: str) -> DriverParameters:
    """Load a single driver by ID, raising its validation error if invalid."""
    for record in load_drivers_raw(db_path):
        if record.get("driver_id") == driver_id:
            return _driver_from_dict(record)
    raise KeyError(f"Driver '{driver_id}' not found in {db_path}")


# ---------------------------------------------------------------------------
# Internal unit-conversion helpers
# ---------------------------------------------------------------------------

def _mh_to_h(val: Optional[float]) -> Optional[float]:
    return val * 1e-3 if val is not None else None


def _mm_to_m(val: Optional[float]) -> Optional[float]:
    return val * 1e-3 if val is not None else None


def _cm2_to_m2(val: Optional[float]) -> Optional[float]:
    return val * 1e-4 if val is not None else None


def _cm3_to_m3(val: Optional[float]) -> Optional[float]:
    return val * 1e-6 if val is not None else None
