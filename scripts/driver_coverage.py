#!/usr/bin/env python3
"""Write docs/DRIVER_COVERAGE.md: what the driver database and price snapshots hold.

Reads only files in the repository (no network):

- data/drivers/<maker>/*.json and data/drivers-curated/<maker>/*.json
- data/drivers-incomplete/ (record files per maker and compression-drivers.json)
- data/prices/<source>/<date>.jsonl (the newest file per source)
- data/sources.json (source registry: status and notes)
- conflict CSVs under data/drivers-incomplete/conflicts/ and data/sources/*/conflicts.csv

The output depends on the data only, not on the day it is run, so a rerun on
unchanged data gives the same file.

    python3 scripts/driver_coverage.py [--root .] [--out docs/DRIVER_COVERAGE.md]
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.parse import urlparse

FIELDS = [
    ("fs_hz", "Fs"), ("re_ohm", "Re"), ("qes", "Qes"), ("qms", "Qms"), ("qts", "Qts"),
    ("sd_m2", "Sd"), ("mms_kg", "Mms"), ("bl_tm", "Bl"), ("le_h", "Le"), ("xmax_m", "Xmax"),
    ("cms_m_per_n", "Cms"), ("rms_kg_per_s", "Rms"), ("nominal_impedance_ohm", "Nominal impedance"),
    ("power_w", "Power (rated)"), ("peak_power_w", "Power (program or peak)"),
    ("sensitivity_db_1w1m", "Sensitivity (1 W, 1 m)"),
]
# Match levels auto mode ranks on (horn_drivers.prices.RANKABLE_MATCHES).
RANKABLE = ("exact", "manual")


def load_records(folder: Path) -> List[dict]:
    records = []
    if not folder.is_dir():
        return records
    for path in sorted(folder.glob("*/*.json")):
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(record, dict) and "parameters" in record:
            record["_folder"] = path.parent.name
            records.append(record)
    return records


def host_to_source(registry: List[dict]) -> Dict[str, str]:
    hosts = {}
    for entry in registry:
        host = urlparse(entry.get("base_url") or "").netloc.lower()
        if host:
            hosts.setdefault(host.removeprefix("www."), entry["id"])
    return hosts


def source_of_url(url: str, hosts: Dict[str, str]) -> str:
    host = urlparse(url or "").netloc.lower().removeprefix("www.")
    return hosts.get(host, host or "unknown")


def field_source(record: dict, field: str, hosts: Dict[str, str]) -> str:
    """The source id a stored value came from."""
    entry = (record.get("field_provenance") or {}).get(field)
    if entry and entry.get("source_id"):
        return entry["source_id"]
    if record.get("record_schema_version", 1) == 1:
        return source_of_url(record.get("parameter_source", ""), hosts)
    return "unknown"


def primary_source(record: dict, hosts: Dict[str, str]) -> str:
    """The source of the Thiele-Small set (the source of Fs)."""
    return field_source(record, "fs_hz", hosts)


def retrieved_on(record: dict) -> Optional[str]:
    entry = (record.get("field_provenance") or {}).get("fs_hz") or {}
    return entry.get("retrieved_on") or (record.get("provenance") or {}).get("retrieved_on")


def latest_price_files(prices_dir: Path) -> Dict[str, Path]:
    files = {}
    if not prices_dir.is_dir():
        return files
    for folder in sorted(p for p in prices_dir.iterdir() if p.is_dir() and p.name != "fx"):
        snapshots = sorted(folder.glob("*.jsonl"))
        if snapshots:
            files[folder.name] = snapshots[-1]
    return files


def load_offers(path: Path) -> List[dict]:
    offers = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line:
            offers.append(json.loads(line))
    return offers


def count_csv_rows(paths: Iterable[Path]) -> Dict[str, int]:
    counts = {}
    for path in paths:
        with path.open(newline="") as handle:
            counts[str(path)] = sum(1 for _ in csv.DictReader(handle))
    return counts


def table(headers: List[str], rows: List[List]) -> List[str]:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(cell) for cell in row) + " |")
    return lines


def pct(part: int, whole: int) -> str:
    return f"{100.0 * part / whole:.0f} %" if whole else "n/a"


def build_report(root: Path) -> str:
    registry = json.loads((root / "data/sources.json").read_text())["sources"]
    hosts = host_to_source(registry)
    scraped = load_records(root / "data/drivers")
    curated = load_records(root / "data/drivers-curated")
    incomplete = load_records(root / "data/drivers-incomplete")
    compression_file = root / "data/drivers-incomplete/compression-drivers.json"
    compression_legacy = 0
    if compression_file.is_file():
        raw = json.loads(compression_file.read_text())
        items = raw.get("drivers", raw) if isinstance(raw, dict) else raw
        compression_legacy = len(items)
    records = scraped + curated
    total = len(records)
    out: List[str] = []
    dates = sorted(d for d in (retrieved_on(r) for r in records) if d)
    price_files = latest_price_files(root / "data/prices")
    offers_by_source = {sid: load_offers(path) for sid, path in price_files.items()}
    price_dates = sorted(o.get("retrieved_on") for offers in offers_by_source.values() for o in offers
                         if o.get("retrieved_on"))
    as_of = max(dates[-1:] + price_dates[-1:]) if (dates or price_dates) else "unknown"

    out += ["# Driver database coverage", "",
            "Generated by `scripts/driver_coverage.py` (`just driver-coverage`) from the files in the "
            f"repository. Newest data date: {as_of}. Do not edit by hand; rerun the script after a data "
            "or price refresh.", ""]

    # Summary
    versions = Counter(r.get("record_schema_version", 1) for r in records)
    types = Counter(r.get("driver_type", "unknown") for r in records)
    out += ["## Summary", ""]
    out += table(["Item", "Count"], [
        ["Usable records (data/drivers)", len(scraped)],
        ["Curated records (data/drivers-curated)", len(curated)],
        ["Schema version 1 records", versions.get(1, 0)],
        ["Schema version 2 records", versions.get(2, 0)],
        ["Cone drivers", types.get("cone", 0)],
        ["Other or unknown driver type", total - types.get("cone", 0)],
        ["Records held back in data/drivers-incomplete", len(incomplete)],
        ["Legacy entries in drivers-incomplete/compression-drivers.json", compression_legacy],
    ])
    out.append("")

    # Records per maker and source
    by_maker: Dict[str, Counter] = defaultdict(Counter)
    for record in records:
        by_maker[record.get("manufacturer") or record["_folder"]][primary_source(record, hosts)] += 1
    sources = sorted({s for c in by_maker.values() for s in c})
    out += ["## Records per maker and source", "",
            "The source is where the Thiele-Small set (Fs) came from. A record can hold other fields "
            "from other sources; see the next table.", ""]
    rows = []
    for maker in sorted(by_maker, key=str.lower):
        counts = by_maker[maker]
        rows.append([maker] + [counts.get(s, "") for s in sources] + [sum(counts.values())])
    rows.append(["Total"] + [sum(by_maker[m].get(s, 0) for m in by_maker) for s in sources] + [total])
    out += table(["Maker"] + sources + ["Total"], rows)
    out.append("")

    held = Counter(r.get("manufacturer") or r["_folder"] for r in incomplete)
    if held:
        out += ["Records held back in `data/drivers-incomplete` per maker:", ""]
        out += table(["Maker", "Records"], [[m, held[m]] for m in sorted(held, key=str.lower)])
        out.append("")

    # Fields coverage
    out += ["## Field coverage", "",
            "Records that hold each field, and the sources the stored values came from.", ""]
    rows = []
    for field, label in FIELDS:
        holders = [r for r in records if r["parameters"].get(field) is not None]
        cone = sum(1 for r in holders if r.get("driver_type") == "cone")
        by_source = Counter(field_source(r, field, hosts) for r in holders)
        detail = ", ".join(f"{s} {n}" for s, n in sorted(by_source.items(), key=lambda kv: (-kv[1], kv[0])))
        rows.append([label, len(holders), pct(len(holders), total), cone, detail])
    band = [r for r in records if r.get("usable_f_low_hz") is not None and r.get("usable_f_high_hz") is not None]
    band_sources = Counter(field_source(r, "usable_f_low_hz", hosts) for r in band)
    rows.append(["Usable band", len(band), pct(len(band), total),
                 sum(1 for r in band if r.get("driver_type") == "cone"),
                 ", ".join(f"{s} {n}" for s, n in sorted(band_sources.items(), key=lambda kv: (-kv[1], kv[0])))])
    sized = sum(1 for r in records if r.get("nominal_diameter"))
    rows.append(["Nominal size", sized, pct(sized, total), sum(1 for r in records if r.get("nominal_diameter")
                                                               and r.get("driver_type") == "cone"), ""])
    out += table(["Field", "Records", "Share", "Cone records", "Values by source"], rows)
    out.append("")
    both = sum(1 for r in records if r["parameters"].get("power_w") is not None
               and r["parameters"].get("sensitivity_db_1w1m") is not None)
    non_aes = sum(1 for r in records
                  if "power_standard" in ((r.get("field_provenance") or {}).get("power_w") or {}))
    out += [f"Records with both rated power and sensitivity: {both} ({pct(both, total)}). "
            f"Rated power marked as not AES: {non_aes}.", ""]

    # Prices
    out += ["## Prices", "",
            "From the newest snapshot file per source in `data/prices/`. Only offers with "
            "`match_confidence` exact or manual are used by auto mode; a list price (`price_kind` list) is used "
            "only when no shop offer exists in the region.", ""]
    known = {r["driver_id"] for r in records}
    rows = []
    for sid in sorted(offers_by_source):
        offers = offers_by_source[sid]
        exact = [o for o in offers if o.get("match_confidence") in RANKABLE]
        drivers = {o["driver_id"] for o in exact}
        regions = sorted({o.get("region") or "?" for o in offers})
        currencies = sorted({o.get("currency") or "?" for o in offers})
        kinds = sorted({o.get("price_kind", "offer") for o in offers})
        dates_ = sorted({o.get("retrieved_on") for o in offers if o.get("retrieved_on")})
        unknown = len(drivers - known)
        rows.append([sid, ", ".join(regions), ", ".join(currencies), len(offers), len(exact), len(drivers),
                     unknown, ", ".join(kinds), dates_[-1] if dates_ else "", price_files[sid].name])
    out += table(["Source", "Region", "Currency", "Offers", "Exact or manual", "Drivers", "Unknown driver ids",
                  "Price kind", "Newest offer", "File"], rows)
    out.append("")

    region_drivers: Dict[str, set] = defaultdict(set)
    any_driver = set()
    for offers in offers_by_source.values():
        for o in offers:
            if o.get("match_confidence") in RANKABLE and o["driver_id"] in known:
                region_drivers[o.get("region") or "?"].add(o["driver_id"])
                any_driver.add(o["driver_id"])
    regions_per_driver = Counter()
    for region, ids in region_drivers.items():
        for driver_id in ids:
            regions_per_driver[driver_id] += 1
    cone_ids = {r["driver_id"] for r in records if r.get("driver_type") == "cone"}
    rows = [[region, len(ids), len(ids & cone_ids)] for region, ids in sorted(region_drivers.items())]
    rows.append(["Any region", len(any_driver), len(any_driver & cone_ids)])
    rows.append(["Two or more regions", sum(1 for n in regions_per_driver.values() if n >= 2),
                 sum(1 for d, n in regions_per_driver.items() if n >= 2 and d in cone_ids)])
    out += ["Drivers with at least one exact or manual offer, per region:", ""]
    out += table(["Region", "Drivers", "Cone drivers"], rows)
    out.append("")

    maker_of = {r["driver_id"]: r.get("manufacturer") or r["_folder"] for r in records}
    maker_total = Counter(maker_of.values())
    priced = Counter(maker_of[d] for d in any_driver)
    rows = [[m, priced[m], maker_total[m], pct(priced[m], maker_total[m])]
            for m in sorted(priced, key=lambda m: (-priced[m], m.lower()))]
    out += ["Drivers with at least one exact or manual offer, per maker:", ""]
    out += table(["Maker", "Priced", "Records", "Share"], rows)
    out.append("")

    # Freshness
    out += ["## Freshness", "",
            "Date the Thiele-Small set was read (field provenance of Fs, or the record provenance "
            "for schema version 1). Legacy loudspeakerdatabase.com records carry the date of the "
            "commit that wrote them.", ""]
    months = Counter((d or "unknown")[:7] for d in (retrieved_on(r) for r in records))
    out += table(["Month read", "Records"], [[m, months[m]] for m in sorted(months)])
    out.append("")

    # Conflicts
    stored = 0
    conflicted_records = 0
    for record in records:
        n = sum(len(e.get("conflicts") or []) for e in (record.get("field_provenance") or {}).values()
                if isinstance(e, dict))
        stored += n
        conflicted_records += bool(n)
    csv_paths = sorted((root / "data/drivers-incomplete/conflicts").glob("*.csv")) + \
        sorted((root / "data/sources").glob("*/conflicts.csv"))
    csv_counts = count_csv_rows(csv_paths)
    out += ["## Conflicts", "",
            f"Stored conflicts (old values kept under `field_provenance.<field>.conflicts`): {stored} "
            f"on {conflicted_records} records. The validator prints each as a warning.", ""]
    if csv_counts:
        out += table(["Conflict list", "Rows"],
                     [[f"`{Path(p).relative_to(root)}`", n] for p, n in csv_counts.items()])
        out.append("")
    review = sorted((root / "data/sources").glob("*/prefix_review.csv"))
    for path, n in count_csv_rows(review).items():
        out += [f"Name matches waiting for a manual decision in `{Path(path).relative_to(root)}`: {n}.", ""]

    # Sources
    out += ["## Sources blocked or not in use", "",
            "Sources in `data/sources.json` that are not approved, or approved but not run, with the "
            "recorded reason.", ""]
    rows = []
    for entry in registry:
        status = entry.get("approval_status")
        if status == "approved" and not entry.get("status_note"):
            continue
        if status == "approved":
            reason = entry.get("status_note", "")
            if "fail" not in reason.lower() and "not run" not in reason.lower():
                continue
        else:
            reason = entry.get("status_note") or entry.get("approval_note") or ""
        rows.append([entry["id"], entry.get("kind", ""), status, reason.replace("|", "/").replace("\n", " ")])
    out += table(["Source", "Kind", "Status", "Reason"], rows)
    out.append("")
    return "\n".join(out)


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", default=".", help="repository root")
    parser.add_argument("--out", default="docs/DRIVER_COVERAGE.md")
    args = parser.parse_args(argv)
    root = Path(args.root)
    text = build_report(root)
    out = Path(args.out)
    if not out.is_absolute():
        out = root / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
