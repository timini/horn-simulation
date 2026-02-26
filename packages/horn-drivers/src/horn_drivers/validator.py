"""Physics-based validation of loudspeaker driver T-S parameters.

Three tiers of validation:
  Tier 1 — Exact identities (must hold within numerical precision)
  Tier 2 — Cross-parameter consistency (10% tolerance)
  Tier 3 — Size-aware plausibility ranges per nominal diameter

Usage:
    horn-validate-drivers [--db data/drivers.json]
"""

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Physical constants
RHO = 1.204    # kg/m^3  air density at 20 °C
C = 343.0      # m/s     speed of sound at 20 °C

ESSENTIAL_FIELDS = ["fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg", "le_h"]

# ---------------------------------------------------------------------------
# Tier 3 — size-aware plausibility ranges
# ---------------------------------------------------------------------------

RANGES_BY_DIAMETER: Dict[str, Dict[str, Tuple[float, float]]] = {
    "6in": {
        "fs_hz": (60.0, 150.0),
        "sd_m2": (0.012, 0.015),
        "mms_kg": (0.006, 0.025),
        "bl_tm": (4.0, 12.0),
        "le_h": (0.0001, 0.002),
    },
    "8in": {
        "fs_hz": (50.0, 120.0),
        "sd_m2": (0.019, 0.023),
        "mms_kg": (0.012, 0.050),
        "bl_tm": (5.0, 18.0),
        "le_h": (0.0002, 0.003),
    },
    "10in": {
        "fs_hz": (40.0, 80.0),
        "sd_m2": (0.031, 0.038),
        "mms_kg": (0.020, 0.090),
        "bl_tm": (8.0, 20.0),
        "le_h": (0.0003, 0.004),
    },
    "12in": {
        "fs_hz": (35.0, 65.0),
        "sd_m2": (0.048, 0.056),
        "mms_kg": (0.035, 0.200),
        "bl_tm": (10.0, 28.0),
        "le_h": (0.0004, 0.005),
    },
    "15in": {
        "fs_hz": (25.0, 50.0),
        "sd_m2": (0.078, 0.092),
        "mms_kg": (0.060, 0.250),
        "bl_tm": (12.0, 35.0),
        "le_h": (0.0005, 0.006),
    },
    "18in": {
        "fs_hz": (20.0, 45.0),
        "sd_m2": (0.110, 0.135),
        "mms_kg": (0.100, 0.400),
        "bl_tm": (15.0, 40.0),
        "le_h": (0.0008, 0.008),
    },
}

RANGES_COMPRESSION: Dict[str, Tuple[float, float]] = {
    "fs_hz": (200.0, 2000.0),
    "sd_m2": (0.0005, 0.005),
    "mms_kg": (0.001, 0.010),
    "bl_tm": (3.0, 15.0),
    "le_h": (0.00005, 0.002),
}

# Fallback combined ranges (used when driver_type/diameter unknown)
RANGES: Dict[str, Tuple[float, float]] = {
    "fs_hz": (15.0, 2000.0),
    "re_ohm": (1.0, 40.0),
    "bl_tm": (0.5, 60.0),
    "sd_m2": (0.0005, 0.25),
    "mms_kg": (0.0001, 1.0),
    "le_h": (0.00001, 0.02),
}

# Per-type ranges (kept for backward compatibility with tests)
RANGES_BY_TYPE = {
    "compression": {
        "fs_hz": (200.0, 2000.0),
        "re_ohm": (1.0, 20.0),
        "bl_tm": (1.0, 30.0),
        "sd_m2": (0.0005, 0.01),
        "mms_kg": (0.0001, 0.1),
        "le_h": (0.00001, 0.005),
    },
    "cone": {
        "fs_hz": (20.0, 300.0),
        "re_ohm": (1.0, 20.0),
        "bl_tm": (4.0, 40.0),
        "sd_m2": (0.005, 0.15),
        "mms_kg": (0.005, 0.5),
        "le_h": (0.0001, 0.01),
    },
}


@dataclass
class ValidationResult:
    """Result of validating a single driver."""
    driver_id: str
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    computed: Dict[str, float] = field(default_factory=dict)
    score: float = 1.0


def compute_derived_params(params: dict) -> dict:
    """Compute derived physical quantities from T-S parameters.

    Returns dict with: Cms, Vas, eta0, SPL_1W, EBP, Qes_computed, BL_Mms_ratio.
    """
    derived: Dict[str, float] = {}

    fs = params.get("fs_hz", 0.0)
    re = params.get("re_ohm", 0.0)
    bl = params.get("bl_tm", 0.0)
    sd = params.get("sd_m2", 0.0)
    mms = params.get("mms_kg", 0.0)

    if fs <= 0 or mms <= 0:
        return derived

    omega_s = 2.0 * math.pi * fs

    # Cms = 1 / (Mms * omega_s^2)
    cms = 1.0 / (mms * omega_s ** 2)
    derived["Cms"] = cms

    # Vas = rho * c^2 * Cms * Sd^2
    if sd > 0:
        vas = RHO * C ** 2 * cms * sd ** 2
        derived["Vas_m3"] = vas
        derived["Vas_litres"] = vas * 1000.0

    # Qes from T-S: Qes = (omega_s * Mms * Re) / BL^2
    if bl > 0 and re > 0:
        qes_computed = (omega_s * mms * re) / (bl ** 2)
        derived["Qes_computed"] = qes_computed

    # Efficiency (eta0) — method A: from BL, Sd, Mms, Re
    if bl > 0 and sd > 0 and re > 0:
        eta0_a = (RHO * bl ** 2 * sd ** 2) / (2.0 * math.pi * C * mms ** 2 * re)
        derived["eta0_A"] = eta0_a

    # Efficiency (eta0) — method B: from fs, Vas, Qes
    qes = params.get("qes")
    if qes and qes > 0 and "Vas_m3" in derived:
        eta0_b = (4.0 * math.pi ** 2 * fs ** 3 * derived["Vas_m3"]) / (C ** 3 * qes)
        derived["eta0_B"] = eta0_b

    # SPL_1W = 112.2 + 10*log10(eta0)
    if "eta0_A" in derived and derived["eta0_A"] > 0:
        derived["SPL_1W"] = 112.2 + 10.0 * math.log10(derived["eta0_A"])

    # EBP = fs / Qes
    if qes and qes > 0:
        derived["EBP"] = fs / qes

    # BL/Mms ratio
    if bl > 0:
        derived["BL_Mms_ratio"] = bl / mms

    return derived


def validate_driver(driver: dict) -> ValidationResult:
    """Validate a single driver with tiered physics checks.

    Returns a ValidationResult with errors, warnings, computed values,
    and a confidence score.
    """
    driver_id = driver.get("driver_id", "?")
    params = driver.get("parameters", driver)
    driver_type = driver.get("driver_type", "unknown")
    nominal_diameter = driver.get("nominal_diameter")

    result = ValidationResult(driver_id=driver_id)
    penalties = 0.0

    # --- Essential field checks ---
    for fld in ESSENTIAL_FIELDS:
        val = params.get(fld)
        if val is None or val == 0.0:
            result.errors.append(f"Missing or zero: {fld}")
            penalties += 0.2

    # If essential fields missing, can't do physics checks
    if result.errors:
        result.score = max(0.0, 1.0 - penalties)
        return result

    # --- Tier 1: Exact identities ---

    # 1. Q factor consistency: 1/Qts = 1/Qms + 1/Qes
    qms = params.get("qms")
    qes = params.get("qes")
    qts = params.get("qts")
    if qms and qes and qts and qms > 0 and qes > 0 and qts > 0:
        reciprocal_sum = 1.0 / qms + 1.0 / qes
        q_diff = abs(1.0 / qts - reciprocal_sum)
        if q_diff > 0.1:
            result.errors.append(
                f"Q factor identity violated: 1/Qts={1.0/qts:.4f} != "
                f"1/Qms+1/Qes={reciprocal_sum:.4f} (diff={q_diff:.4f})"
            )
            penalties += 0.15
        elif q_diff > 0.01:
            result.warnings.append(
                f"Q factor rounding: 1/Qts={1.0/qts:.4f} vs "
                f"1/Qms+1/Qes={reciprocal_sum:.4f} (diff={q_diff:.4f})"
            )
            penalties += 0.02

    # 2. Cms = 1 / (Mms * omega_s^2) — this is a definition, always holds
    #    (just compute and store, no check needed as it's derived)

    # --- Tier 2: Cross-parameter consistency (10% tolerance) ---

    derived = compute_derived_params(params)
    result.computed = derived

    # 3. Qes cross-check
    if qes and qes > 0 and "Qes_computed" in derived:
        ratio = derived["Qes_computed"] / qes
        if abs(ratio - 1.0) > 0.10:
            result.warnings.append(
                f"Qes cross-check: stated={qes:.3f}, "
                f"computed={derived['Qes_computed']:.3f} "
                f"(ratio={ratio:.2f})"
            )
            penalties += 0.05

    # 5. Dual efficiency cross-check (must agree within 0.5 dB)
    if "eta0_A" in derived and "eta0_B" in derived:
        if derived["eta0_A"] > 0 and derived["eta0_B"] > 0:
            spl_a = 10.0 * math.log10(derived["eta0_A"])
            spl_b = 10.0 * math.log10(derived["eta0_B"])
            if abs(spl_a - spl_b) > 0.5:
                result.warnings.append(
                    f"Efficiency mismatch: method A={spl_a + 112.2:.1f} dB, "
                    f"method B={spl_b + 112.2:.1f} dB "
                    f"(delta={abs(spl_a - spl_b):.2f} dB)"
                )
                penalties += 0.05

    # --- Tier 3: Size-aware plausibility ---

    # Select ranges based on driver type and diameter
    if driver_type == "compression":
        ranges = RANGES_COMPRESSION
    elif nominal_diameter and nominal_diameter in RANGES_BY_DIAMETER:
        ranges = RANGES_BY_DIAMETER[nominal_diameter]
    elif driver_type == "cone":
        ranges = RANGES_BY_TYPE["cone"]
    else:
        ranges = RANGES

    for fld, (lo, hi) in ranges.items():
        val = params.get(fld)
        if val is not None and val != 0.0:
            if val < lo or val > hi:
                result.warnings.append(
                    f"{fld}={val} outside expected range [{lo}, {hi}]"
                )
                penalties += 0.03

    # Q factor positivity
    for qf in ["qms", "qes", "qts"]:
        val = params.get(qf)
        if val is not None and val <= 0:
            result.errors.append(f"{qf}={val} must be positive")
            penalties += 0.1

    # 8. Re vs nominal impedance
    re = params.get("re_ohm", 0.0)
    z_nom = params.get("nominal_impedance_ohm")
    if z_nom and z_nom > 0 and re > 0:
        ratio = re / z_nom
        if ratio < 0.55 or ratio > 0.85:
            result.warnings.append(
                f"Re/Z_nom ratio={ratio:.2f} outside [0.55, 0.85]"
            )
            penalties += 0.02

    # 10. BL/Mms ratio for cone drivers
    # Smaller drivers (6"-10") have higher BL/Mms ratios due to lower Mms
    bl = params.get("bl_tm", 0.0)
    mms = params.get("mms_kg", 0.0)
    if driver_type == "cone" and bl > 0 and mms > 0:
        bl_mms = bl / mms
        if bl_mms < 50.0 or bl_mms > 800.0:
            result.warnings.append(
                f"BL/Mms ratio={bl_mms:.1f} outside [50, 800] T·m/kg"
            )
            penalties += 0.02

    # 11. Air-load mass floor: Mms > 2.67 * rho * (Sd/pi)^1.5
    sd = params.get("sd_m2", 0.0)
    if sd > 0 and mms > 0:
        air_load_floor = 2.67 * RHO * (sd / math.pi) ** 1.5
        if mms < air_load_floor:
            result.warnings.append(
                f"Mms={mms:.6f} kg below air-load floor {air_load_floor:.6f} kg"
            )
            penalties += 0.05

    result.score = max(0.0, 1.0 - penalties)
    return result


def validate_database(db_path: str) -> List[ValidationResult]:
    """Validate entire driver database. Returns list of ValidationResults.

    Accepts either a directory (v3 per-driver files) or a single JSON file.
    """
    from horn_drivers.loader import load_drivers_raw

    drivers = load_drivers_raw(db_path)
    print(f"Validating {len(drivers)} drivers from {db_path}\n")

    results = []
    for driver in drivers:
        result = validate_driver(driver)
        results.append(result)

    errors = sum(len(r.errors) for r in results)
    warnings = sum(len(r.warnings) for r in results)

    for r in results:
        for msg in r.errors:
            print(f"  ERROR [{r.driver_id}]: {msg}")
        for msg in r.warnings:
            print(f"  WARN  [{r.driver_id}]: {msg}")

    print(f"\n{'=' * 60}")
    print(f"Results: {len(drivers)} drivers, {errors} errors, {warnings} warnings")

    scores = [r.score for r in results]
    if scores:
        print(f"Confidence scores: min={min(scores):.2f}, "
              f"mean={sum(scores)/len(scores):.2f}, max={max(scores):.2f}")

    if errors == 0:
        print("PASS: All drivers pass validation.")
    else:
        print("FAIL: Some drivers have validation errors.")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Validate driver database with physics-based checks."
    )
    parser.add_argument(
        "--db", type=str, default="data/drivers",
        help="Path to driver database directory or JSON file (default: data/drivers)",
    )
    args = parser.parse_args()

    results = validate_database(args.db)
    errors = sum(len(r.errors) for r in results)
    sys.exit(1 if errors > 0 else 0)


if __name__ == "__main__":
    main()
