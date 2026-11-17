"""Driver record schema version 2: identity, per-field provenance, writers.

Version 1 records (the scraped and curated ones) keep a single
``parameter_source`` URL and a free-form ``provenance`` block. Version 2 adds:

``identity``
    ``{"manufacturer": <canonical name>, "model": str,
    "variant": {"impedance_ohm": float | None, "coil": str | None},
    "aliases": [str, ...]}``. Aliases hold retailer SKUs and spellings.

``field_provenance``
    ``{<field>: {"source_id", "url", "retrieved_on", "method"}}`` for every
    stored parameter (and the band fields). ``method`` is one of
    ``page``, ``pdf``, ``derived``, ``manual``. An entry may carry
    ``conflicts``: a list of ``{"source_id", "url", "value"}`` that lost the
    per-field conflict rule.

The top-level ``manufacturer`` and ``model_name`` stay, so every reader of
version 1 records works unchanged. Prices never live in a driver record; see
``horn_drivers.prices``.
"""

import datetime as _dt
import json
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from horn_drivers.identity import canonical_manufacturer, is_canonical_manufacturer

RECORD_SCHEMA_VERSION_V2 = 2
SUPPORTED_RECORD_SCHEMA_VERSIONS = (1, 2)
PROVENANCE_METHODS = ("page", "pdf", "derived", "manual")
TOP_LEVEL_FIELDS = ("usable_f_low_hz", "usable_f_high_hz")
ESSENTIAL_PROVENANCE_FIELDS = ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg", "le_h")


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def driver_id_for(manufacturer: str, model: str, impedance_ohm: Optional[float] = None) -> str:
    """Stable id in the scraper's scheme (canonical manufacturer plus model).

    Pass *impedance_ohm* when the model name is shared by several impedance
    variants; it is appended as ``-<n>ohm``.
    """
    canonical = canonical_manufacturer(manufacturer) or manufacturer
    base = _slug(f"{canonical}-{model}")
    if impedance_ohm:
        base += "-" + _slug(f"{impedance_ohm:g}") + "ohm"
    return base


def provenance_entry(source_id: str, url: str, retrieved_on: str, method: str, **extra) -> dict:
    if method not in PROVENANCE_METHODS:
        raise ValueError(f"method must be one of {PROVENANCE_METHODS}, not {method!r}")
    _dt.date.fromisoformat(retrieved_on)
    return {"source_id": source_id, "url": url, "retrieved_on": retrieved_on, "method": method, **extra}


def make_record(*, manufacturer: str, model: str, parameters: Dict[str, float],
                field_provenance: Dict[str, dict], driver_type: str = "unknown",
                impedance_ohm: Optional[float] = None, coil: Optional[str] = None,
                aliases: Iterable[str] = (), nominal_diameter: Optional[str] = None,
                usable_f_low_hz: Optional[float] = None, usable_f_high_hz: Optional[float] = None,
                driver_id: Optional[str] = None, id_includes_impedance: bool = False,
                notes: Iterable[str] = ()) -> dict:
    """Build a version 2 driver record. Parameters must already be in SI units.

    Raises ValueError for an unknown manufacturer or a structurally invalid
    record; physics checks are left to the loader and validator.
    """
    canonical = canonical_manufacturer(manufacturer)
    if canonical is None:
        raise ValueError(f"Unknown manufacturer {manufacturer!r}")
    if impedance_ohm is None:
        impedance_ohm = parameters.get("nominal_impedance_ohm")
    record = {
        "driver_id": driver_id or driver_id_for(canonical, model, impedance_ohm if id_includes_impedance else None),
        "record_schema_version": RECORD_SCHEMA_VERSION_V2,
        "manufacturer": canonical,
        "model_name": model,
        "identity": {"manufacturer": canonical, "model": model,
                     "variant": {"impedance_ohm": impedance_ohm, "coil": coil},
                     "aliases": sorted(set(aliases))},
        "driver_type": driver_type,
        "parameters": dict(parameters),
        "field_provenance": dict(field_provenance),
    }
    if nominal_diameter:
        record["nominal_diameter"] = nominal_diameter
    if usable_f_low_hz is not None and usable_f_high_hz is not None:
        record["usable_f_low_hz"], record["usable_f_high_hz"] = usable_f_low_hz, usable_f_high_hz
    sources = {entry.get("url") for entry in field_provenance.values() if entry.get("url")}
    record["parameter_source"] = sorted(sources)[0] if len(sources) == 1 else None
    if notes:
        record["provenance"] = {"notes": list(notes)}
    errors, _ = check_structure(record)
    if errors:
        raise ValueError(f"{record['driver_id']}: " + "; ".join(errors))
    return record


def check_structure(record: dict) -> Tuple[List[str], List[str]]:
    """Schema checks for one record: (errors, warnings).

    Version 1 records pass with no errors. Version 2 records need a canonical
    identity and a complete provenance entry for each stored essential field.
    """
    errors: List[str] = []
    warnings: List[str] = []
    version = record.get("record_schema_version", 1)
    if version not in SUPPORTED_RECORD_SCHEMA_VERSIONS:
        return [f"record_schema_version {version!r} is not supported (expected one of {SUPPORTED_RECORD_SCHEMA_VERSIONS})"], []
    if version == 1:
        return errors, warnings

    identity = record.get("identity")
    if not isinstance(identity, dict):
        errors.append("version 2 record needs an identity block")
    else:
        manufacturer = identity.get("manufacturer")
        if not is_canonical_manufacturer(manufacturer):
            errors.append(f"identity.manufacturer {manufacturer!r} is not a canonical name (horn_drivers.identity)")
        elif record.get("manufacturer") not in (None, manufacturer):
            errors.append(f"manufacturer {record.get('manufacturer')!r} differs from identity.manufacturer {manufacturer!r}")
        if not identity.get("model"):
            errors.append("identity.model is empty")
        if not isinstance(identity.get("aliases", []), list):
            errors.append("identity.aliases must be a list")
        if not isinstance(identity.get("variant", {}), dict):
            errors.append("identity.variant must be an object")

    params = record.get("parameters") or {}
    provenance = record.get("field_provenance")
    if not isinstance(provenance, dict):
        errors.append("version 2 record needs a field_provenance map")
        return errors, warnings
    for name, entry in provenance.items():
        if name not in params and name not in TOP_LEVEL_FIELDS:
            warnings.append(f"field_provenance names {name}, which the record does not store")
        if not isinstance(entry, dict):
            errors.append(f"field_provenance.{name} must be an object")
            continue
        missing = [key for key in ("source_id", "url", "retrieved_on", "method") if not entry.get(key)]
        if missing:
            errors.append(f"field_provenance.{name} lacks {', '.join(missing)}")
            continue
        if entry["method"] not in PROVENANCE_METHODS:
            errors.append(f"field_provenance.{name}.method {entry['method']!r} not in {PROVENANCE_METHODS}")
        try:
            _dt.date.fromisoformat(str(entry["retrieved_on"]))
        except ValueError:
            errors.append(f"field_provenance.{name}.retrieved_on {entry['retrieved_on']!r} is not an ISO date")
        if entry["method"] != "manual" and not str(entry["url"]).startswith(("http://", "https://")):
            errors.append(f"field_provenance.{name}.url is not an http(s) URL")
        for conflict in entry.get("conflicts") or []:
            if isinstance(conflict, dict):
                warnings.append(f"{name}: {entry['source_id']} value {params.get(name)!r} kept; "
                                f"{conflict.get('source_id')} gave {conflict.get('value')!r}")
    for name in ESSENTIAL_PROVENANCE_FIELDS:
        if params.get(name) is not None and name not in provenance:
            errors.append(f"{name} has no field_provenance entry")
    for name in list(params) + [f for f in TOP_LEVEL_FIELDS if record.get(f) is not None]:
        if name not in ESSENTIAL_PROVENANCE_FIELDS and name not in provenance:
            warnings.append(f"{name} has no field_provenance entry")
    return errors, warnings


def parameter_source(record: dict) -> Optional[str]:
    """Single source URL for display: v1 parameter_source, else the Fs provenance URL."""
    if record.get("parameter_source"):
        return record["parameter_source"]
    provenance = record.get("field_provenance") or {}
    for name in ("fs_hz", *sorted(provenance)):
        entry = provenance.get(name)
        if isinstance(entry, dict) and entry.get("url"):
            return entry["url"]
    return None


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def write_record(db_dir: Path, record: dict) -> Path:
    """Write ``db_dir/<manufacturer>/<driver_id>.json`` atomically."""
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]*", record["driver_id"]):
        raise ValueError(f"driver_id {record['driver_id']!r} is not a safe file name")
    path = Path(db_dir) / record["manufacturer"] / f"{record['driver_id']}.json"
    _atomic_write(path, json.dumps(record, indent=4) + "\n")
    return path


def write_incomplete(incomplete_dir: Path, record: dict, reason: str) -> Path:
    """Write a record that failed loading or validation, with the reason, for review."""
    return write_record(Path(incomplete_dir), {**record, "incomplete_reason": reason})
