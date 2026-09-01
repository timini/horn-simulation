"""Driver pre-screening for horn optimization.

Filters a driver database to candidates suitable for a target horn
specification before running FEM simulations.
"""

import argparse
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np

from horn_core.parameters import DriverParameters

_SPEED_OF_SOUND = 343.0  # m/s at ~20°C


@dataclass
class PrescreenConfig:
    """Configuration for driver pre-screening."""
    target_f_low_hz: float
    target_f_high_hz: float
    mouth_radius_m: Optional[float] = None
    length_m: Optional[float] = None
    min_ebp: float = 50.0
    ka_max: float = 2 * math.pi  # ~6.28, absolute cap on throat ka at f_high
    min_nominal_diameter_in: Optional[float] = None
    max_nominal_diameter_in: Optional[float] = None


@dataclass
class PrescreenResult:
    """Result of driver pre-screening."""
    drivers: List[DriverParameters]
    throat_radius_m: float
    throat_radii_m: List[float]
    count: int

    def to_dict(self) -> dict:
        return {
            "drivers": [d.driver_id for d in self.drivers],
            "throat_radius_m": self.throat_radius_m,
            "throat_radii_m": self.throat_radii_m,
            "count": self.count,
        }


def _parse_diameter_inches(d: Optional[str]) -> Optional[float]:
    """Parse a nominal diameter string like '4in' or '6.5' to inches."""
    if not d:
        return None
    m = re.match(r"(\d+(?:\.\d+)?)", d)
    return float(m.group(1)) if m else None


def prescreen_drivers(
    drivers: List[DriverParameters],
    config: PrescreenConfig,
) -> PrescreenResult:
    """Filter drivers to candidates suitable for the target horn.

    Throat radius is decoupled from driver size — the throat is a property
    of the horn (constrained by acoustics at f_high), not the driver.

    Filtering criteria:
    1. fs_hz < target_f_low_hz * 1.5 -- driver resonance in or near target band
    2. fs_hz / qes > min_ebp -- horn suitability (Efficiency Bandwidth Product)
    3. Optional nominal diameter filter (user-specified min/max)
    4. Driver must physically fit in the horn mouth

    Throat radii are derived from acoustic constraints (ka ≤ ka_max at f_high),
    not from driver Sd.

    Args:
        drivers: Full list of drivers from the database.
        config: Pre-screening configuration.

    Returns:
        PrescreenResult with filtered drivers and representative throat radius.
    """
    candidates = []

    # Max driver radius: driver must fit inside the horn mouth
    if config.mouth_radius_m is not None and config.mouth_radius_m > 0:
        max_driver_radius = config.mouth_radius_m
    else:
        ideal_mouth = _SPEED_OF_SOUND / (2 * math.pi * config.target_f_low_hz)
        max_driver_radius = ideal_mouth

    for drv in drivers:
        # 1. Resonance frequency check
        if drv.fs_hz >= config.target_f_low_hz * 1.5:
            continue

        # 2. EBP check (horn suitability)
        if drv.qes is not None and drv.qes > 0:
            ebp = drv.fs_hz / drv.qes
            if ebp < config.min_ebp:
                continue
        else:
            continue

        # 3. Optional nominal diameter filter
        if config.min_nominal_diameter_in is not None or config.max_nominal_diameter_in is not None:
            dia_in = _parse_diameter_inches(drv.nominal_diameter)
            if dia_in is not None:
                if config.min_nominal_diameter_in is not None and dia_in < config.min_nominal_diameter_in:
                    continue
                if config.max_nominal_diameter_in is not None and dia_in > config.max_nominal_diameter_in:
                    continue

        # 4. Driver must fit in the horn mouth
        drv_radius = math.sqrt(drv.effective_throat_area / math.pi)
        if drv_radius > max_driver_radius:
            continue

        candidates.append(drv)

    if not candidates:
        return PrescreenResult(drivers=[], throat_radius_m=0.0, throat_radii_m=[], count=0)

    # Derive throat radii from acoustic constraints (ka ≤ ka_max at f_high)
    a_acoustic_max = _SPEED_OF_SOUND * config.ka_max / (
        2 * math.pi * config.target_f_high_hz
    )

    # Acoustic range: fractions of the maximum acoustic throat radius
    acoustic_radii = [a_acoustic_max * f for f in [0.3, 0.65, 1.0]]

    # Also include driver-matched radii for direct-coupling scenarios
    driver_radii = [
        math.sqrt(d.effective_throat_area / math.pi)
        for d in candidates if d.effective_throat_area > 0
    ]
    direct_radii = [r for r in driver_radii if r <= a_acoustic_max]

    all_radii = sorted(set(
        round(r, 6) for r in acoustic_radii + direct_radii
    ))

    # Keep at most 5 to limit combinatorial explosion
    if len(all_radii) > 5:
        indices = np.linspace(0, len(all_radii) - 1, 5, dtype=int)
        all_radii = [all_radii[i] for i in indices]

    representative_radius = all_radii[len(all_radii) // 2]
    throat_radii_m = all_radii

    return PrescreenResult(
        drivers=candidates,
        throat_radius_m=representative_radius,
        throat_radii_m=throat_radii_m,
        count=len(candidates),
    )


def main():
    """CLI for driver pre-screening."""
    parser = argparse.ArgumentParser(
        description="Pre-screen drivers for a target horn specification.",
    )
    parser.add_argument("--drivers-db", required=True, help="Driver database JSON file.")
    parser.add_argument("--target-f-low", type=float, required=True, help="Target low frequency (Hz).")
    parser.add_argument("--target-f-high", type=float, required=True, help="Target high frequency (Hz).")
    parser.add_argument("--mouth-radius", type=float, default=None, help="Horn mouth radius (m). Optional for fullauto mode.")
    parser.add_argument("--length", type=float, default=None, help="Horn length (m). Optional for fullauto mode.")
    parser.add_argument("--min-ebp", type=float, default=50.0, help="Minimum EBP threshold.")
    parser.add_argument("--ka-max", type=float, default=2 * math.pi,
                        help="Absolute cap on throat ka at f_high (default 2π ≈ 6.28).")
    parser.add_argument("--min-diameter", type=float, default=None,
                        help="Minimum driver nominal diameter (inches).")
    parser.add_argument("--max-diameter", type=float, default=None,
                        help="Maximum driver nominal diameter (inches).")
    parser.add_argument("--output", type=str, default="prescreen_result.json", help="Output JSON file.")
    args = parser.parse_args()

    from horn_drivers.loader import load_drivers

    drivers = load_drivers(args.drivers_db)
    print(f"Loaded {len(drivers)} drivers from {args.drivers_db}")

    config = PrescreenConfig(
        target_f_low_hz=args.target_f_low,
        target_f_high_hz=args.target_f_high,
        mouth_radius_m=args.mouth_radius,
        length_m=args.length,
        min_ebp=args.min_ebp,
        ka_max=args.ka_max,
        min_nominal_diameter_in=args.min_diameter,
        max_nominal_diameter_in=args.max_diameter,
    )

    result = prescreen_drivers(drivers, config)

    output = json.dumps(result.to_dict(), indent=2)
    Path(args.output).write_text(output)
    print(f"Pre-screening complete: {result.count} drivers passed")
    print(f"Representative throat radius: {result.throat_radius_m:.4f} m")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
