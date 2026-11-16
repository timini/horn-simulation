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

# Zeros of J1': the first non-axisymmetric (1,1) duct mode cuts on at ka = 1.841
# and the first radial (0,1) mode at ka = 3.832. Above 3.832 inside the band a
# plane-wave (Webster) prediction of the throat is no longer meaningful, so
# default throat sampling stops there. Above 1.841 the throat is flagged:
# an axisymmetric drive does not excite the (1,1) mode, but any asymmetry does.
THROAT_KA_FIRST_RADIAL_MODE = 3.8317
THROAT_KA_FIRST_ASYMMETRIC_MODE = 1.8412
# Radii within this relative distance of a kept radius are treated as duplicates.
THROAT_RADIUS_DEDUP_RELATIVE = 0.01


@dataclass
class PrescreenConfig:
    """Configuration for driver pre-screening."""
    target_f_low_hz: float
    target_f_high_hz: float
    mouth_radius_m: Optional[float] = None
    length_m: Optional[float] = None
    min_ebp: float = 50.0
    # Cap on throat ka at the top of the target band, where the score is measured.
    ka_max: float = THROAT_KA_FIRST_RADIAL_MODE
    ka_soft: float = THROAT_KA_FIRST_ASYMMETRIC_MODE
    min_nominal_diameter_in: Optional[float] = None
    max_nominal_diameter_in: Optional[float] = None
    throat_fractions: Optional[List[float]] = None  # override default [0.3,0.65,1.0]


@dataclass
class PrescreenResult:
    """Result of driver pre-screening."""
    drivers: List[DriverParameters]
    throat_radius_m: float
    throat_radii_m: List[float]
    count: int
    throat_flags: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "drivers": [d.driver_id for d in self.drivers],
            "throat_radius_m": self.throat_radius_m,
            "throat_radii_m": self.throat_radii_m,
            "count": self.count,
            "throat_flags": self.throat_flags,
        }


def dedupe_radii(radii, relative=THROAT_RADIUS_DEDUP_RELATIVE):
    """Sorted radii with any radius within ``relative`` of the previous kept one dropped."""
    kept = []
    for radius in sorted(radii):
        if not kept or radius > kept[-1] * (1 + relative):
            kept.append(radius)
    return kept


def chamber_matched_throat_radius(driver, f_high_hz):
    """Throat radius whose front-chamber corner c St / (2 pi V) sits at f_high.

    None when the driver has no front chamber (supplied as zero, or unknown
    type). With an estimated or supplied chamber, a smaller throat puts the
    chamber's low-pass corner inside the band.
    """
    from horn_core.driver_interface import resolve_front_chamber
    volume, _, _ = resolve_front_chamber(driver)
    if volume <= 0:
        return None
    area = 2 * math.pi * f_high_hz * volume / _SPEED_OF_SOUND
    return math.sqrt(area / math.pi)


def throat_flags(radii, drivers, f_high_hz, ka_soft=THROAT_KA_FIRST_ASYMMETRIC_MODE):
    """Per sampled radius: ka at f_high, the soft-flag state and drivers with Sd < St."""
    flags = {}
    for radius in radii:
        ka = 2 * math.pi * f_high_hz * radius / _SPEED_OF_SOUND
        area = math.pi * radius**2
        flags[f"{radius:.6f}"] = {
            "ka_at_f_high": ka,
            "above_first_asymmetric_mode": ka > ka_soft,
            "drivers_with_sd_below_throat": sorted(d.driver_id for d in drivers if d.sd_m2 < area),
        }
    return flags


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

    Throat radii are derived from acoustic constraints (ka <= ka_max at
    f_high, default 3.83, the first radial duct mode), plus driver-matched
    radii (St = Sd) and chamber-matched radii (front-chamber corner at
    f_high) that fall under the cap. Radii within 1% are merged. Radii above
    ka_soft (1.84) and drivers with Sd < St are flagged, not removed.

    Args:
        drivers: Full list of drivers from the database.
        config: Pre-screening configuration.

    Returns:
        PrescreenResult with filtered drivers and representative throat radius.
    """
    if config.throat_fractions is not None and (
        not config.throat_fractions or any(not math.isfinite(f) or not 0 < f <= 1 for f in config.throat_fractions)
    ):
        raise ValueError("Throat fractions must be finite and within (0, 1]")
    if not math.isfinite(config.ka_max) or config.ka_max <= 0:
        raise ValueError("Throat ka cap must be positive and finite")
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
    fractions = config.throat_fractions if config.throat_fractions else [0.3, 0.65, 1.0]
    acoustic_radii = [a_acoustic_max * f for f in fractions]

    # Also include driver-matched radii for direct-coupling scenarios
    driver_radii = [
        math.sqrt(d.effective_throat_area / math.pi)
        for d in candidates if d.effective_throat_area > 0
    ]
    direct_radii = [r for r in driver_radii if r <= a_acoustic_max]
    # A front chamber (estimated or supplied) favours throats large enough to
    # put its low-pass corner at or above f_high; capped like every other radius.
    chamber_radii = [r for r in (chamber_matched_throat_radius(d, config.target_f_high_hz) for d in candidates)
                     if r is not None]
    chamber_radii = [min(r, a_acoustic_max) for r in chamber_radii]

    # round() can lift a radius above the cap; floor to the cap after rounding.
    all_radii = dedupe_radii(
        min(round(r, 6), a_acoustic_max) for r in acoustic_radii + direct_radii + chamber_radii
    )

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
        throat_flags=throat_flags(throat_radii_m, candidates, config.target_f_high_hz, config.ka_soft),
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
    parser.add_argument("--ka-max", type=float, default=THROAT_KA_FIRST_RADIAL_MODE,
                        help="Cap on sampled throat ka at the target high frequency (default 3.83, first radial mode).")
    parser.add_argument("--min-diameter", type=float, default=None,
                        help="Minimum driver nominal diameter (inches).")
    parser.add_argument("--max-diameter", type=float, default=None,
                        help="Maximum driver nominal diameter (inches).")
    parser.add_argument("--throat-fractions", type=str, default=None,
                        help="Comma-separated throat radius fractions of a_max (default '0.3,0.65,1.0'). "
                             "Use smaller values for phase-plug-tiny throats, e.g. '0.1,0.25,0.5,0.75,1.0'.")
    parser.add_argument("--output", type=str, default="prescreen_result.json", help="Output JSON file.")
    args = parser.parse_args()

    from horn_drivers.loader import load_drivers

    drivers = load_drivers(args.drivers_db)
    print(f"Loaded {len(drivers)} drivers from {args.drivers_db}")

    throat_fractions = None
    if args.throat_fractions:
        throat_fractions = [float(x) for x in args.throat_fractions.split(",")]

    config = PrescreenConfig(
        target_f_low_hz=args.target_f_low,
        target_f_high_hz=args.target_f_high,
        mouth_radius_m=args.mouth_radius,
        length_m=args.length,
        min_ebp=args.min_ebp,
        ka_max=args.ka_max,
        min_nominal_diameter_in=args.min_diameter,
        max_nominal_diameter_in=args.max_diameter,
        throat_fractions=throat_fractions,
    )

    result = prescreen_drivers(drivers, config)

    output = json.dumps({**result.to_dict(), "ka_max": config.ka_max, "ka_soft": config.ka_soft,
                         "ka_reference": "target_f_high"}, indent=2)
    Path(args.output).write_text(output)
    print(f"Pre-screening complete: {result.count} drivers passed")
    print(f"Representative throat radius: {result.throat_radius_m:.4f} m")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
