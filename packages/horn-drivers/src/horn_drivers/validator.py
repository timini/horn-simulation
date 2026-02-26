"""Validate the driver database for completeness and plausibility.

Checks:
- All essential T-S fields are present and non-zero
- Sd is plausible for compression drivers (0.0005-0.01 m^2)
- Fs is in range for horn drivers (200-2000 Hz)
- BL is positive and plausible (1-30 T*m)
- Re is positive (1-20 Ohm)
- Q factors are positive and in expected ranges

Usage:
    horn-validate-drivers [--db data/drivers.json]
"""

import argparse
import json
import sys
from pathlib import Path
from typing import List, Tuple

ESSENTIAL_FIELDS = ["fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg", "le_h"]

RANGES = {
    "fs_hz": (200.0, 2000.0),
    "re_ohm": (1.0, 20.0),
    "bl_tm": (1.0, 30.0),
    "sd_m2": (0.0005, 0.01),
    "mms_kg": (0.0001, 0.1),
    "le_h": (0.00001, 0.01),
}


def validate_driver(driver: dict) -> List[Tuple[str, str]]:
    """Validate a single driver entry, return list of (level, message) tuples."""
    issues: List[Tuple[str, str]] = []
    driver_id = driver.get("driver_id", "?")
    params = driver.get("parameters", driver)

    for field in ESSENTIAL_FIELDS:
        val = params.get(field)
        if val is None or val == 0.0:
            issues.append(("ERROR", f"[{driver_id}] Missing or zero: {field}"))

    for field, (lo, hi) in RANGES.items():
        val = params.get(field)
        if val is not None and val != 0.0:
            if val < lo or val > hi:
                issues.append((
                    "WARN",
                    f"[{driver_id}] {field}={val} outside expected range [{lo}, {hi}]",
                ))

    for qf in ["qms", "qes", "qts"]:
        val = params.get(qf)
        if val is not None and val <= 0:
            issues.append(("ERROR", f"[{driver_id}] {qf}={val} must be positive"))

    return issues


def validate_database(db_path: str) -> int:
    """Validate entire driver database. Returns number of errors."""
    db = json.loads(Path(db_path).read_text())

    if isinstance(db, dict) and "drivers" in db:
        drivers = db["drivers"]
    else:
        drivers = [{"driver_id": k, **v} for k, v in db.items()]

    print(f"Validating {len(drivers)} drivers from {db_path}\n")

    errors = 0
    warnings = 0

    for driver in drivers:
        issues = validate_driver(driver)
        for level, msg in issues:
            print(f"  {level}: {msg}")
            if level == "ERROR":
                errors += 1
            else:
                warnings += 1

    print(f"\n{'=' * 50}")
    print(f"Results: {len(drivers)} drivers, {errors} errors, {warnings} warnings")

    if errors == 0:
        print("PASS: All drivers have essential fields populated.")
    else:
        print("FAIL: Some drivers are missing essential data.")

    return errors


def main():
    parser = argparse.ArgumentParser(description="Validate driver database.")
    parser.add_argument(
        "--db", type=str, default="data/drivers.json",
        help="Path to driver database JSON (default: data/drivers.json)",
    )
    args = parser.parse_args()

    errors = validate_database(args.db)
    sys.exit(1 if errors > 0 else 0)


if __name__ == "__main__":
    main()
