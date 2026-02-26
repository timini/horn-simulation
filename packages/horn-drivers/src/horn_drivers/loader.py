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
from pathlib import Path
from typing import List, Optional

from horn_core.parameters import DriverParameters


def _driver_from_dict(d: dict) -> DriverParameters:
    """Convert a driver dict (from JSON) into a DriverParameters instance.

    Handles both v1 (flat) and v2 (nested ``parameters``) formats, and
    converts common non-SI units (mH → H, cm² → m², mm → m) on the fly.
    """
    # v2 format nests T-S params under "parameters"
    params = d.get("parameters", d)

    # Unit conversions — source may use convenience units
    le_h = params.get("le_h") or _mh_to_h(params.get("le_mh"))
    sd_m2 = params.get("sd_m2") or params.get("sd_sq_meters")
    xmax_m = params.get("xmax_m") or _mm_to_m(params.get("xmax_mm"))

    return DriverParameters(
        driver_id=d.get("driver_id", "unknown"),
        manufacturer=d.get("manufacturer", ""),
        model_name=d.get("model_name", ""),
        fs_hz=params["fs_hz"],
        re_ohm=params.get("re_ohm") or params.get("re_ohms", 0.0),
        bl_tm=params.get("bl_tm", 0.0),
        sd_m2=sd_m2 or 0.0,
        mms_kg=params.get("mms_kg", 0.0),
        le_h=le_h or 0.0,
        qms=params.get("qms"),
        qes=params.get("qes"),
        qts=params.get("qts"),
        driver_type=d.get("driver_type"),
        nominal_diameter=d.get("nominal_diameter"),
        xmax_m=xmax_m,
        nominal_impedance_ohm=params.get("nominal_impedance_ohm"),
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
    """Load all drivers from a database (directory or JSON file).

    Supports:
      - v3 directory: ``db_path/{Manufacturer}/{driver-id}.json``
      - v2 single file: ``{ "schema_version": 2, "drivers": [...] }``
      - v1 single file: ``{ "driver_id": { ... }, ... }``
    """
    return [_driver_from_dict(d) for d in load_drivers_raw(db_path)]


def load_driver(db_path: str, driver_id: str) -> DriverParameters:
    """Load a single driver by ID from a database."""
    for driver in load_drivers(db_path):
        if driver.driver_id == driver_id:
            return driver
    raise KeyError(f"Driver '{driver_id}' not found in {db_path}")


# ---------------------------------------------------------------------------
# Internal unit-conversion helpers
# ---------------------------------------------------------------------------

def _mh_to_h(val: Optional[float]) -> Optional[float]:
    return val * 1e-3 if val is not None else None


def _mm_to_m(val: Optional[float]) -> Optional[float]:
    return val * 1e-3 if val is not None else None
