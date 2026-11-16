"""Add one driver to a database from datasheet values with explicit units.

Every option names its unit, so a datasheet can be typed in as printed:
Sd in cm^2, Mms in grams, Vas in litres, Le in mH, Xmax in mm. The values
are converted to SI, checked by the same loader and validator the pipeline
uses, and written with their source URL and retrieval date. Nothing is
written when a check fails.

Example:
    horn-add-driver --db data/drivers-curated \\
        --manufacturer Acme --model "AC10 8ohm" --type cone \\
        --fs-hz 45 --re-ohm 5.6 --bl-tm 15.2 --sd-cm2 346 --vas-l 52 \\
        --le-mh 0.9 --qms 5.1 --qes 0.38 --qts 0.35 --xmax-mm 4 \\
        --power-w 300 --f-low-hz 50 --f-high-hz 3000 \\
        --source-url https://example.com/ac10.pdf
"""

import argparse
import datetime
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from horn_drivers.loader import C_AIR, RHO_AIR
from horn_drivers.scraper import RECORD_SCHEMA_VERSION, slugify
from horn_drivers.validator import validate_driver

# (option, record key, factor to SI, help)
PARAMETER_OPTIONS: List[Tuple[str, str, float, str]] = [
    ("fs-hz", "fs_hz", 1.0, "free-air resonance, Hz"),
    ("re-ohm", "re_ohm", 1.0, "DC resistance, ohm"),
    ("bl-tm", "bl_tm", 1.0, "force factor, T*m"),
    ("sd-cm2", "sd_m2", 1e-4, "effective diaphragm area, cm^2"),
    ("mms-g", "mms_kg", 1e-3, "moving mass including air load, g"),
    ("mmd-g", "mmd_kg", 1e-3, "diaphragm moving mass without air load, g"),
    ("le-mh", "le_h", 1e-3, "voice-coil inductance, mH"),
    ("qms", "qms", 1.0, "mechanical Q"),
    ("qes", "qes", 1.0, "electrical Q"),
    ("qts", "qts", 1.0, "total Q"),
    ("rms-kg-s", "rms_kg_per_s", 1.0, "mechanical resistance, kg/s (N*s/m)"),
    ("cms-um-n", "cms_m_per_n", 1e-6, "suspension compliance, um/N"),
    ("xmax-mm", "xmax_m", 1e-3, "linear excursion (one way), mm"),
    ("z-ohm", "nominal_impedance_ohm", 1.0, "nominal impedance, ohm"),
    ("power-w", "power_w", 1.0, "nominal / AES power, W"),
    ("peak-power-w", "peak_power_w", 1.0, "program or peak power, W"),
    ("sensitivity-db", "sensitivity_db_1w1m", 1.0, "sensitivity, dB at 1 W / 1 m"),
    ("front-chamber-cm3", "front_chamber_volume_m3", 1e-6,
     "air volume between diaphragm and throat (adapter or phase plug), cm^3; 0 means none"),
    ("front-chamber-area-cm2", "front_chamber_area_m2", 1e-4,
     "cross-section of the front chamber when it acts as a short duct, cm^2"),
    ("rear-chamber-cm3", "rear_chamber_volume_m3", 1e-6, "sealed rear chamber volume, cm^3; omit for open back"),
]
# Zero is a meaningful value for these: no front chamber at all.
ZERO_ALLOWED = {"front_chamber_volume_m3"}
RECORD_OPTIONS = [
    ("f-low-hz", "usable_f_low_hz", "lower limit of the manufacturer's frequency range, Hz"),
    ("f-high-hz", "usable_f_high_hz", "upper limit of the manufacturer's frequency range, Hz"),
]


def _manufacturer_dir(manufacturer: str) -> str:
    name = re.sub(r"[^A-Za-z0-9]+", "", manufacturer)
    if not name:
        raise ValueError("manufacturer must contain letters or digits")
    return name


def build_record(args: argparse.Namespace) -> Tuple[dict, List[str], List[str]]:
    """Convert CLI values to an SI record. Returns (record, supplied, derived)."""
    params: Dict[str, float] = {}
    supplied: List[str] = []
    derived: List[str] = []
    for option, key, factor, _ in PARAMETER_OPTIONS:
        value = getattr(args, option.replace("-", "_"))
        if value is not None:
            if not math.isfinite(value) or value < 0 or (value == 0 and key not in ZERO_ALLOWED):
                raise ValueError(f"--{option} must be finite and positive")
            params[key] = float(f"{value * factor:.12g}")
            supplied.append(key)

    if args.vas_l is not None:
        if args.mms_g is not None:
            raise ValueError("give either --mms-g or --vas-l, not both")
        if "sd_m2" not in params or "fs_hz" not in params:
            raise ValueError("--vas-l needs --sd-cm2 and --fs-hz to derive Mms")
        if not math.isfinite(args.vas_l) or args.vas_l <= 0:
            raise ValueError("--vas-l must be finite and positive")
        params["vas_l"] = args.vas_l
        # Cms = Vas / (rho c^2 Sd^2); Mms = 1 / ((2 pi fs)^2 Cms)
        cms = args.vas_l * 1e-3 / (RHO_AIR * C_AIR ** 2 * params["sd_m2"] ** 2)
        params["mms_kg"] = float(f"{1.0 / ((2 * math.pi * params['fs_hz']) ** 2 * cms):.6g}")
        supplied.append("vas_l")
        derived.append("mms_kg")

    if args.exit_diameter_mm is not None:
        if args.exit_diameter_mm <= 0:
            raise ValueError("--exit-diameter-mm must be positive")
        params["exit_area_m2"] = float(f"{math.pi * (args.exit_diameter_mm * 1e-3 / 2) ** 2:.12g}")
        derived.append("exit_area_m2")
        supplied.append("exit_diameter_mm")

    model = args.model.strip()
    manufacturer = args.manufacturer.strip()
    record: dict = {
        "driver_id": args.driver_id or slugify(manufacturer, model),
        "manufacturer": manufacturer,
        "model_name": model,
        "driver_type": args.type,
        "parameters": params,
    }
    if args.type != "unknown":
        record["type_source"] = args.type_source or f"datasheet: {args.source_url}"
    if args.nominal_diameter:
        record["nominal_diameter"] = args.nominal_diameter
    for option, key, _ in RECORD_OPTIONS:
        value = getattr(args, option.replace("-", "_"))
        if value is not None:
            record[key] = value
            supplied.append(key)
    record["parameter_source"] = args.source_url
    record["interface_model"] = None
    record["record_schema_version"] = RECORD_SCHEMA_VERSION
    provenance = {
        "retrieved_on": args.retrieved_on,
        "source_kind": args.source_kind,
        "capture_method": "Entered with horn-add-driver from the datasheet values given on the command line.",
        "supplied_fields": sorted(supplied),
        "derived_fields": sorted(derived),
    }
    if args.exit_diameter_mm is not None:
        provenance["exit_diameter_mm"] = args.exit_diameter_mm
    if "usable_f_low_hz" in record:
        provenance["band_definition"] = args.band_definition
    if args.note:
        provenance["notes"] = list(args.note)
    record["provenance"] = provenance
    return record, supplied, derived


def check_record(record: dict) -> Tuple[List[str], List[str]]:
    """Run the loader and validator; return (errors, warnings)."""
    result = validate_driver(record)
    return result.errors, result.warnings


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add a driver from datasheet values. Each option states its unit.")
    parser.add_argument("--db", default="data/drivers-curated",
                        help="driver database directory (default: data/drivers-curated)")
    parser.add_argument("--manufacturer", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--driver-id", help="default: slug of manufacturer and model")
    parser.add_argument("--type", choices=("cone", "compression", "unknown"), default="unknown")
    parser.add_argument("--type-source", help="where the driver type comes from (default: the datasheet)")
    parser.add_argument("--nominal-diameter", help="e.g. 12in")
    for option, _, _, text in PARAMETER_OPTIONS:
        parser.add_argument(f"--{option}", type=float, help=text)
    parser.add_argument("--vas-l", type=float, help="equivalent volume, litres (used to derive Mms when Mms is not given)")
    parser.add_argument("--exit-diameter-mm", type=float, help="compression-driver exit (throat) diameter, mm")
    for option, _, text in RECORD_OPTIONS:
        parser.add_argument(f"--{option}", type=float, help=text)
    parser.add_argument("--band-definition",
                        default="Manufacturer nominal frequency range, not a validated horn-loaded band.",
                        help="what --f-low-hz/--f-high-hz mean in the source")
    parser.add_argument("--source-url", required=True, help="datasheet or product page the values came from")
    parser.add_argument("--source-kind", default="manufacturer_datasheet",
                        help="e.g. manufacturer_datasheet, distributor_copy_of_manufacturer_datasheet, measurement")
    parser.add_argument("--retrieved-on", default=datetime.date.today().isoformat(),
                        help="date the source was read, YYYY-MM-DD (default: today)")
    parser.add_argument("--note", action="append", help="free-text note, repeatable")
    parser.add_argument("--force", action="store_true", help="overwrite an existing record with the same id")
    parser.add_argument("--dry-run", action="store_true", help="print the record instead of writing it")
    args = parser.parse_args(argv)
    try:
        datetime.date.fromisoformat(args.retrieved_on)
    except ValueError:
        parser.error("--retrieved-on must be YYYY-MM-DD")
    if not re.match(r"^https?://", args.source_url):
        parser.error("--source-url must be an http(s) URL")
    return args


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    try:
        record, _, _ = build_record(args)
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    errors, warnings = check_record(record)
    for message in warnings:
        print(f"WARN: {message}", file=sys.stderr)
    if errors:
        for message in errors:
            print(f"ERROR: {message}", file=sys.stderr)
        print("Nothing written.", file=sys.stderr)
        return 1
    text = json.dumps(record, indent=2, ensure_ascii=False) + "\n"
    if args.dry_run:
        print(text, end="")
        return 0
    path = Path(args.db) / _manufacturer_dir(record["manufacturer"]) / f"{record['driver_id']}.json"
    if path.exists() and not args.force:
        print(f"ERROR: {path} exists; pass --force to replace it", file=sys.stderr)
        return 1
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
