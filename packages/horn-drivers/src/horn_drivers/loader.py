"""Driver database loader for horn simulation.

Loads compression driver Thiele-Small parameters from a JSON database,
converting stored units to SI and deriving missing parameters.
"""

import json
from pathlib import Path
from typing import List, Optional

from horn_core.parameters import DriverParameters


def _driver_from_dict(d: dict) -> DriverParameters:
    """Convert a driver dict (from JSON) into a DriverParameters instance.

    Handles both v1 (flat) and v2 (nested ``parameters``) formats, and
    converts common non-SI units (mH -> H, cm^2 -> m^2, mm -> m) on the fly.
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
        xmax_m=xmax_m,
        nominal_impedance_ohm=params.get("nominal_impedance_ohm"),
    )


def load_drivers(db_path: str) -> List[DriverParameters]:
    """Load all drivers from a JSON database file.

    Supports v1 (dict-of-dicts keyed by driver_id) and v2 (object with
    ``schema_version`` and ``drivers`` array) formats.
    """
    raw = json.loads(Path(db_path).read_text())

    # v2 format: { "schema_version": 2, "drivers": [...] }
    if isinstance(raw, dict) and "drivers" in raw:
        return [_driver_from_dict(d) for d in raw["drivers"]]

    # v1 format: { "driver_id": { ... }, ... }
    drivers = []
    for driver_id, data in raw.items():
        data.setdefault("driver_id", driver_id)
        drivers.append(_driver_from_dict(data))
    return drivers


def load_driver(db_path: str, driver_id: str) -> DriverParameters:
    """Load a single driver by ID from a JSON database file."""
    for driver in load_drivers(db_path):
        if driver.driver_id == driver_id:
            return driver
    raise KeyError(f"Driver '{driver_id}' not found in {db_path}")


def _mh_to_h(val: Optional[float]) -> Optional[float]:
    return val * 1e-3 if val is not None else None


def _mm_to_m(val: Optional[float]) -> Optional[float]:
    return val * 1e-3 if val is not None else None
