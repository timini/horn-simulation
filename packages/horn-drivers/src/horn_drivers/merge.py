"""Merge manufacturer records into a driver database without silent overwrites.

Adapters write version 2 records to a staging folder. ``merge_staged`` then
joins each staged record to the existing record for the same driver
(canonical manufacturer, normalised model, impedance) and applies the
per-field conflict rule from the source plan:

- the manufacturer value wins;
- where the old value differs by more than 5 % (10 % for Le and Xmax), the
  old value is kept in ``field_provenance[field]["conflicts"]`` and a row is
  written to the conflict report;
- a field the manufacturer page lacks keeps the old value and its source;
- the matched record keeps its ``driver_id``, so runs and tests that name it
  still work, and git history keeps the old file.

Staged records that match nothing are added. Old records that match nothing
are left as they are. Staged records with an impedance suffix in their id
(``bc-8ndl51-8ohm``) lose the suffix when the model has a single impedance
variant and the plain id is free.
"""

import csv
import datetime as _dt
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from horn_drivers.identity import canonical_manufacturer, normalise_model, record_identity
from horn_drivers.records import driver_id_for, provenance_entry, write_record

COMPARED_FIELDS = ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg", "le_h", "qms", "qes", "qts", "xmax_m",
                   "nominal_impedance_ohm")
LOOSE_FIELDS = {"le_h": 0.10, "xmax_m": 0.10}
DEFAULT_TOLERANCE = 0.05
CARRIED_FIELDS = ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg", "le_h", "qms", "qes", "qts", "xmax_m",
                  "nominal_impedance_ohm", "rms_kg_per_s", "cms_m_per_n")
REPORT_FIELDS = ("driver_id", "manufacturer", "model", "impedance_ohm", "field", "manufacturer_value",
                 "previous_value", "difference_pct", "tolerance_pct", "manufacturer_url", "previous_source",
                 "previous_url")


@dataclass
class MergeSummary:
    added: List[str] = field(default_factory=list)
    updated: List[str] = field(default_factory=list)
    duplicates: List[Tuple[str, str]] = field(default_factory=list)   # (staged id, reason)
    conflicts: List[dict] = field(default_factory=list)
    untouched_previous: int = 0
    outputs: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"added": len(self.added), "updated": len(self.updated), "duplicates": self.duplicates,
                "conflict_rows": len(self.conflicts),
                "drivers_with_conflicts": len({row["driver_id"] for row in self.conflicts}),
                "untouched_previous": self.untouched_previous, "outputs": self.outputs}


def _impedance(record: dict) -> Optional[float]:
    value = (record_identity(record).get("variant") or {}).get("impedance_ohm")
    if value is None:
        value = (record.get("parameters") or {}).get("nominal_impedance_ohm")
    return float(value) if value else None


def _model_key(record: dict) -> str:
    identity = record_identity(record)
    return normalise_model(identity.get("model") or record.get("model_name") or "", identity.get("manufacturer"))


def relative_difference(new: float, old: float) -> float:
    if new == old:
        return 0.0
    scale = max(abs(new), abs(old))
    return abs(new - old) / scale if scale else 0.0


def _previous_source(record: dict) -> Tuple[str, Optional[str], str]:
    """(source_id, url, retrieved_on) describing where an old v1 record came from."""
    url = record.get("parameter_source")
    provenance = record.get("provenance") or {}
    retrieved = provenance.get("retrieved_on") or "1970-01-01"
    if url and "loudspeakerdatabase.com" in url:
        source = "loudspeakerdatabase"
    elif record.get("field_provenance"):
        entry = record["field_provenance"].get("fs_hz") or next(iter(record["field_provenance"].values()))
        return entry.get("source_id", "unknown"), entry.get("url"), entry.get("retrieved_on", retrieved)
    else:
        source = provenance.get("source_kind") or "unknown"
    return source, url, retrieved


def merge_one(new: dict, old: dict) -> Tuple[dict, List[dict]]:
    """Merged record (manufacturer values win) and the conflict rows."""
    merged = json.loads(json.dumps(new))
    merged["driver_id"] = old["driver_id"]
    params = merged["parameters"]
    provenance = merged["field_provenance"]
    old_params = old.get("parameters") or {}
    new_source = next(iter(provenance.values()))["source_id"] if provenance else None
    old_source, old_url, old_date = _previous_source(old)
    same_source = old.get("record_schema_version") == 2 and old_source == new_source
    rows = []
    if same_source:
        # A rerun of the same adapter: keep conflicts found earlier, add none.
        for name, entry in (old.get("field_provenance") or {}).items():
            if entry.get("conflicts") and name in provenance:
                provenance[name]["conflicts"] = entry["conflicts"]
    else:
        for name in COMPARED_FIELDS:
            new_value, old_value = params.get(name), old_params.get(name)
            if new_value is None or old_value is None:
                continue
            difference = relative_difference(float(new_value), float(old_value))
            tolerance = LOOSE_FIELDS.get(name, DEFAULT_TOLERANCE)
            if difference <= tolerance:
                continue
            provenance[name].setdefault("conflicts", []).append(
                {"source_id": old_source, "url": old_url, "value": old_value, "retrieved_on": old_date})
            rows.append({"driver_id": merged["driver_id"], "manufacturer": merged["manufacturer"],
                         "model": merged["model_name"], "impedance_ohm": _impedance(merged), "field": name,
                         "manufacturer_value": new_value, "previous_value": old_value,
                         "difference_pct": round(100 * difference, 1), "tolerance_pct": round(100 * tolerance),
                         "manufacturer_url": provenance[name]["url"], "previous_source": old_source,
                         "previous_url": old_url})
        for name in CARRIED_FIELDS:
            if params.get(name) is None and old_params.get(name) is not None and old_url:
                params[name] = old_params[name]
                provenance[name] = provenance_entry(old_source, old_url, old_date, "page",
                                                    note="not on the manufacturer page; kept from the previous record")
        notes = list((merged.get("provenance") or {}).get("notes", []))
        notes.append(f"Replaces the {old_source} record retrieved {old_date}. Previous values that differ by more "
                     "than 5 % (10 % for Le and Xmax) are kept under field_provenance conflicts.")
        merged["provenance"] = {**(merged.get("provenance") or {}), "notes": notes}
        if "mms_kg" in provenance and "may be the dry moving mass" in json.dumps(old.get("provenance", {})):
            provenance["mms_kg"]["note"] = "manufacturer Mms; the previous value may have been Mmd"
    aliases = set(merged["identity"].get("aliases", []))
    if old.get("model_name") and old["model_name"] != merged["model_name"]:
        aliases.add(old["model_name"])
    merged["identity"]["aliases"] = sorted(aliases)
    return merged, rows


def load_staged(folder: Path, makers: Optional[List[str]] = None) -> List[dict]:
    """Every record under ``folder/<maker>/*.json``, optionally for some makers only."""
    records = []
    for path in sorted(Path(folder).glob("*/*.json")):
        if makers and path.parent.name not in makers:
            continue
        records.append(json.loads(path.read_text()))
    return records


def merge_staged(staged: List[dict], db_dir: Path, report_csv: Optional[Path] = None,
                 dry_run: bool = False) -> MergeSummary:
    """Merge *staged* v2 records into *db_dir*; write the conflict report CSV."""
    db_dir = Path(db_dir)
    summary = MergeSummary()
    by_maker: Dict[str, List[dict]] = {}
    for record in staged:
        by_maker.setdefault(record["manufacturer"], []).append(record)
    for maker, records in sorted(by_maker.items()):
        existing = load_folder_maker(db_dir, maker)
        existing_ids = {record["driver_id"] for record in existing}
        old_by_model: Dict[str, List[dict]] = {}
        for record in existing:
            old_by_model.setdefault(_model_key(record), []).append(record)
        new_by_model: Dict[str, List[dict]] = {}
        for record in records:
            new_by_model.setdefault(_model_key(record), []).append(record)
        used = set()
        written_ids = set()
        for record in sorted(records, key=lambda r: r["driver_id"]):
            key, impedance = _model_key(record), _impedance(record)
            candidates = [r for r in old_by_model.get(key, []) if r["driver_id"] not in used]
            exact = [r for r in candidates if _impedance(r) == impedance]
            old = exact[0] if exact else None
            if old is None and len(candidates) == 1 and len(new_by_model[key]) == 1 \
                    and not any(_impedance(r) == _impedance(candidates[0]) for r in new_by_model[key]):
                old = candidates[0]
            if old is not None:
                merged, rows = merge_one(record, old)
                used.add(old["driver_id"])
                summary.updated.append(merged["driver_id"])
                summary.conflicts.extend(rows)
            else:
                merged = json.loads(json.dumps(record))
                plain = driver_id_for(maker, merged["model_name"])
                if len(new_by_model[key]) == 1 and plain not in existing_ids and plain not in written_ids:
                    merged["driver_id"] = plain
                if merged["driver_id"] in written_ids or (merged["driver_id"] in existing_ids
                                                          and merged["driver_id"] not in used):
                    summary.duplicates.append((record["driver_id"], "id already taken by another record"))
                    continue
                summary.added.append(merged["driver_id"])
            written_ids.add(merged["driver_id"])
            if not dry_run:
                write_record(db_dir, merged)
        summary.untouched_previous += len(existing_ids - used)
    if report_csv is not None and summary.conflicts and not dry_run:
        report_csv = Path(report_csv)
        report_csv.parent.mkdir(parents=True, exist_ok=True)
        with report_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS)
            writer.writeheader()
            writer.writerows(sorted(summary.conflicts, key=lambda r: (r["driver_id"], r["field"])))
        summary.outputs.append(str(report_csv))
    return summary


def load_folder_maker(db_dir: Path, maker: str) -> List[dict]:
    folder = Path(db_dir) / maker
    if not folder.is_dir():
        return []
    return [json.loads(path.read_text()) for path in sorted(folder.glob("*.json"))]


def default_report_path(report_dir: Path, label: str) -> Path:
    return Path(report_dir) / f"{label}-{_dt.date.today().isoformat()}.csv"
