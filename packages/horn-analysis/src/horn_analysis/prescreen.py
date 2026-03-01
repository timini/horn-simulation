"""Driver pre-screening for horn optimization.

Filters a driver database to candidates suitable for a target horn
specification before running FEM simulations.
"""

import argparse
import json
import math
from dataclasses import dataclass, asdict
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
    horn_load_factor: float = 10.0
    sd_ratio_range: Tuple[float, float] = (0.3, 3.0)


@dataclass
class PrescreenResult:
    """Result of driver pre-screening."""
    drivers: List[DriverParameters]
    throat_radius_m: float
    count: int

    def to_dict(self) -> dict:
        return {
            "drivers": [d.driver_id for d in self.drivers],
            "throat_radius_m": self.throat_radius_m,
            "count": self.count,
        }


def prescreen_drivers(
    drivers: List[DriverParameters],
    config: PrescreenConfig,
) -> PrescreenResult:
    """Filter drivers to candidates suitable for the target horn.

    Filtering criteria:
    1. fs_hz < target_f_low_hz * 1.5 -- driver resonance in or near target band
    2. fs_hz / qes > min_ebp -- horn suitability (Efficiency Bandwidth Product)
    3. Upper freq capability: f_piston * horn_load_factor >= target_f_high
    4. Representative throat radius = median(sqrt(Sd/pi)) of passing drivers
    5. Filter drivers whose effective radius is outside sd_ratio_range of representative

    Args:
        drivers: Full list of drivers from the database.
        config: Pre-screening configuration.

    Returns:
        PrescreenResult with filtered drivers and representative throat radius.
    """
    candidates = []

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
            # Cannot compute EBP, skip
            continue

        # 3. Upper frequency capability from Sd
        #    f_piston = c / (2π·a_eff) where a_eff = √(Sd/π)
        #    Horn loading extends usable range by ~10× above piston breakup
        if drv.sd_m2 > 0:
            a_eff = math.sqrt(drv.sd_m2 / math.pi)
            f_piston = _SPEED_OF_SOUND / (2 * math.pi * a_eff)
            if f_piston * config.horn_load_factor < config.target_f_high_hz:
                continue

        candidates.append(drv)

    if not candidates:
        return PrescreenResult(drivers=[], throat_radius_m=0.0, count=0)

    # 4. Compute representative throat radius
    radii = [math.sqrt(d.sd_m2 / math.pi) for d in candidates if d.sd_m2 > 0]
    if not radii:
        return PrescreenResult(drivers=[], throat_radius_m=0.0, count=0)

    representative_radius = float(np.median(radii))

    # 5. Filter by Sd ratio relative to representative
    lo, hi = config.sd_ratio_range
    filtered = []
    for drv in candidates:
        if drv.sd_m2 <= 0:
            continue
        drv_radius = math.sqrt(drv.sd_m2 / math.pi)
        ratio = drv_radius / representative_radius
        if lo <= ratio <= hi:
            filtered.append(drv)

    # Recompute representative from final set
    if filtered:
        final_radii = [math.sqrt(d.sd_m2 / math.pi) for d in filtered]
        representative_radius = float(np.median(final_radii))

    return PrescreenResult(
        drivers=filtered,
        throat_radius_m=representative_radius,
        count=len(filtered),
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
    parser.add_argument("--horn-load-factor", type=float, default=10.0,
                        help="Multiplier on piston frequency for horn-loaded upper limit estimate.")
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
        horn_load_factor=args.horn_load_factor,
    )

    result = prescreen_drivers(drivers, config)

    output = json.dumps(result.to_dict(), indent=2)
    Path(args.output).write_text(output)
    print(f"Pre-screening complete: {result.count} drivers passed")
    print(f"Representative throat radius: {result.throat_radius_m:.4f} m")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
