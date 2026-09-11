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
import warnings
from pathlib import Path
from typing import List, Optional

from horn_core.parameters import DriverParameters


def _driver_from_dict(d: dict) -> DriverParameters:
    """Convert a driver dict (from JSON) into a DriverParameters instance.

    Handles both v1 (flat) and v2 (nested ``parameters``) formats, and
    converts common non-SI units (mH → H, cm² → m², mm → m) on the fly.
    """
    if d.get("catalogue_status") == "quarantined":
        raise ValueError("Catalogue record quarantined: " + d.get("quarantine_reason", "source identity unverified"))
    # v2 format nests T-S params under "parameters"
    params = d.get("parameters", d)

    # Unit conversions — source may use convenience units
    # Only missing values fall back to aliases. Explicit zero limits must
    # reach validation, and zero inductance is a valid limiting case.
    le_h = params.get("le_h")
    if le_h is None:
        le_h = _mh_to_h(params.get("le_mh"))
    if le_h is None:
        raise ValueError("Missing voice-coil inductance; unknown is not zero")
    sd_m2 = params.get("sd_m2")
    if sd_m2 is None:
        sd_m2 = params.get("sd_sq_meters")
    xmax_m = params.get("xmax_m")
    if xmax_m is None:
        xmax_m = _mm_to_m(params.get("xmax_mm"))
    exit_area_m2 = params.get("exit_area_m2")
    if exit_area_m2 is None:
        exit_area_m2 = _cm2_to_m2(params.get("exit_area_cm2"))
    re_ohm = params.get("re_ohm")
    if re_ohm is None:
        re_ohm = params.get("re_ohms", 0.0)

    return DriverParameters(
        driver_id=d.get("driver_id", "unknown"),
        manufacturer=d.get("manufacturer", ""),
        model_name=d.get("model_name", ""),
        fs_hz=params["fs_hz"],
        re_ohm=re_ohm,
        bl_tm=params.get("bl_tm", 0.0),
        sd_m2=sd_m2 if sd_m2 is not None else 0.0,
        mms_kg=params.get("mms_kg", 0.0),
        le_h=le_h if le_h is not None else 0.0,
        qms=params.get("qms"),
        qes=params.get("qes"),
        qts=params.get("qts"),
        exit_area_m2=exit_area_m2,
        driver_type=d.get("driver_type"),
        nominal_diameter=d.get("nominal_diameter") if d.get("nominal_diameter_verified") is not False else None,
        overall_diameter_m=d.get("overall_diameter_m"),
        xmax_m=xmax_m,
        nominal_impedance_ohm=params.get("nominal_impedance_ohm"),
        power_w=params.get("power_w"),
        peak_power_w=params.get("peak_power_w"),
        usable_f_low_hz=d.get("usable_f_low_hz"),
        usable_f_high_hz=d.get("usable_f_high_hz"),
        parameter_source=d.get("parameter_source"),
        interface_model=d.get("interface_model"),
        mmd_kg=params.get("mmd_kg"),
        rear_load_mass_kg=params.get("rear_load_mass_kg"),
        cms_m_per_n=params.get("cms_m_per_n"),
        rms_kg_per_s=params.get("rms_kg_per_s"),
    )


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
