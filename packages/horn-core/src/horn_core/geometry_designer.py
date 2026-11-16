"""Analytical geometry derivation for fullauto horn design.

Given only a target frequency band, derives optimal horn geometry
(mouth radius, length) using horn acoustics formulas, then generates
a focused grid of candidates for FEM evaluation.

Key formulas:
  - Mouth radius from low-freq target: r_mouth = c0 / (2*pi*f_low)
    (mouth circumference = wavelength at f_low)
  - Length range: quarter-wave to half-wave at f_low
    L_min = c0 / (4*f_low), L_max = c0 / (2*f_low)
  - Simulation freq range: one octave below and half an octave above the
    target band, so the f3 points that the extension score uses are found
    (see derive_simulation_freq_range)
"""

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np

from horn_core.candidates import CandidateGeometry, write_candidates_csv
from horn_core.profiles import profile_geometry_error

# Speed of sound in air at ~20C
C0 = 343.0

DEFAULT_PROFILES = ["conical", "exponential", "hyperbolic", "tractrix", "os", "scaled_tractrix", "cd"]

# With only a maximum size given, search down to this fraction of it.
SIZE_CAP_SEARCH_FRACTION = 0.4


@dataclass
class DerivedGeometry:
    """Stores analytical geometry derivation results."""

    target_f_low: float
    target_f_high: float
    ideal_mouth_radius: float
    mouth_radius_range: tuple  # (min, max)
    length_range: tuple  # (min, max)
    sim_freq_range: tuple  # (min, max) extended for rolloff
    candidate_count: int
    status: str = "candidates_generated"
    reason: Optional[str] = None
    bound_sources: Optional[dict] = None
    size_compromise: Optional[dict] = None


def _search_range(heuristic, fixed, lower, upper):
    """Searched (lo, hi) and the source of each side for one dimension.

    A user minimum or maximum replaces that side of the heuristic range. With
    only a maximum given, the search runs from 0.4 x max (or the heuristic
    minimum, if smaller) up to the maximum, so a cap below the heuristic
    range still returns designs.
    """
    if fixed is not None:
        return (fixed, fixed), {"low": "user_fixed", "high": "user_fixed"}
    lo, hi = heuristic
    sources = {"low": "heuristic", "high": "heuristic"}
    if lower is not None:
        lo, sources["low"] = lower, "user"
    if upper is not None:
        hi, sources["high"] = upper, "user"
        if lower is None:
            lo = min(heuristic[0], SIZE_CAP_SEARCH_FRACTION * upper)
            sources["low"] = "size_cap_fraction" if lo < heuristic[0] else "heuristic"
    return (lo, hi), sources


def size_compromise(target_f_low, mouth_range, length_range):
    """Low-frequency cost of size caps below the heuristic mouth or length range.

    Uses the same rules that set the heuristic ranges: a mouth circumference
    of one wavelength (ka = 1) and a quarter-wave length. Each capped
    dimension gives the lowest frequency it supports; the larger of those is
    the expected low-frequency limit and the shortfall is in octaves above
    the requested f_low. This is a sizing rule, not a prediction; the
    ranked f3_low values are the predictions. None when no cap binds.
    """
    ideal_mouth = derive_mouth_radius(target_f_low)
    quarter_wave = C0 / (4 * target_f_low)
    limits = {}
    if mouth_range[1] < ideal_mouth * (1 - 1e-9):
        limits["mouth_radius"] = {"max_m": mouth_range[1], "ka1_mouth_radius_m": ideal_mouth,
                                  "supported_f_low_hz": C0 / (2 * math.pi * mouth_range[1])}
    if length_range[1] < quarter_wave * (1 - 1e-9):
        limits["length"] = {"max_m": length_range[1], "quarter_wave_length_m": quarter_wave,
                            "supported_f_low_hz": C0 / (4 * length_range[1])}
    if not limits:
        return None
    f_supported = max(v["supported_f_low_hz"] for v in limits.values())
    return {"limits": limits, "expected_f_low_hz": f_supported, "target_f_low_hz": target_f_low,
            "shortfall_octaves": math.log2(f_supported / target_f_low),
            "basis": "ka = 1 mouth and quarter-wave length sizing rules; see ranked f3_low for predictions"}


def derive_mouth_radius(f_low: float) -> float:
    """Derive ideal mouth radius from low-frequency target.

    The mouth circumference should equal the wavelength at f_low
    for efficient low-frequency radiation.

    r_mouth = c0 / (2 * pi * f_low)
    """
    return C0 / (2 * math.pi * f_low)


def derive_mouth_radius_range(
    f_low: float, spread: float = 0.3
) -> tuple:
    """Derive mouth radius range as +/- spread around ideal.

    Returns (min_radius, max_radius).
    """
    ideal = derive_mouth_radius(f_low)
    return (ideal * (1 - spread), ideal * (1 + spread))


def derive_length_range(f_low: float) -> tuple:
    """Derive horn length range from quarter-wave to half-wave at f_low.

    L_min = c0 / (4 * f_low)  (quarter wavelength)
    L_max = c0 / (2 * f_low)  (half wavelength)
    """
    wavelength = C0 / f_low
    return (wavelength / 4, wavelength / 2)


# Sweep margins beyond the target band, in octaves. The extension score term
# clip(0.5 + E/2, 0, 1) saturates at E = 1 octave, so the low side is swept a
# full octave: an f3 below f_low/2 already earns the whole term on its own and
# reporting it as a bound costs nothing. Half an octave was not enough: in the
# wave-2 runs nearly every feasible row had f3_low at the sweep minimum. The
# low side is cheap because the mesh is sized for the top frequency and the
# point count is fixed by num_intervals. The high side stays at half an
# octave because raising the top frequency refines the mesh.
SIM_LOW_MARGIN_OCTAVES = 1.0
SIM_HIGH_MARGIN_OCTAVES = 0.5


def derive_simulation_freq_range(
    f_low: float, f_high: float
) -> tuple:
    """Solver sweep for a target band: SIM_LOW_MARGIN_OCTAVES below f_low
    and SIM_HIGH_MARGIN_OCTAVES above f_high.

    f3 points beyond these limits are reported as sweep bounds
    (f3_low_is_bound, f3_high_is_bound) and the extension as a lower bound.
    """
    sim_min = f_low / 2 ** SIM_LOW_MARGIN_OCTAVES
    sim_max = f_high * 2 ** SIM_HIGH_MARGIN_OCTAVES
    return (sim_min, sim_max)


def _empty_grid_reason(candidates, mouth_radii, lengths, profile_skipped) -> Optional[str]:
    """Why an auto grid produced no candidates, or None if it produced some."""
    if candidates:
        return None
    if not mouth_radii or not lengths:
        return "size_constraints_exclude_search_range"
    if profile_skipped:
        # Some mouth/length points were larger than the throat but every one
        # was impossible for its profile (a tractrix below its natural length).
        return "profile_geometry_infeasible"
    return "no_mouth_larger_than_throat"


def generate_auto_candidates(
    target_f_low: float,
    target_f_high: float,
    throat_radii: List[float],
    mouth_radius: Optional[float] = None,
    length: Optional[float] = None,
    num_mouth_radii: int = 3,
    num_lengths: int = 3,
    profiles: Optional[List[str]] = None,
    max_length: Optional[float] = None,
    min_length: Optional[float] = None,
    max_mouth_radius: Optional[float] = None,
    min_mouth_radius: Optional[float] = None,
    ka_max: float = 2 * math.pi,
) -> tuple:
    """Generate geometry candidates for unified auto mode.

    When mouth_radius or length is provided, it is fixed (1 value).
    When None, the range is derived from target_f_low.

    Args:
        target_f_low: Target low-frequency cutoff (Hz).
        target_f_high: Target high-frequency cutoff (Hz).
        throat_radii: Throat radii (m) from prescreen (1-3 values).
        mouth_radius: Fixed mouth radius (m), or None to derive range.
        length: Fixed horn length (m), or None to derive range.
        num_mouth_radii: Grid points for mouth radius (used when mouth_radius is None).
        num_lengths: Grid points for length (used when length is None).
        profiles: Horn profile types. Defaults to all seven.

    Returns:
        Tuple of (candidates, derived_geometry).
    """
    from horn_core.acoustics import validate_band
    validate_band(target_f_low, target_f_high)
    if num_mouth_radii < 1 or num_lengths < 1:
        raise ValueError("Geometry grid sizes must be positive")
    if not throat_radii or any(not np.isfinite(r) or r <= 0 for r in throat_radii):
        raise ValueError("No usable throat radii: no feasible drivers or invalid radii")
    if not np.isfinite(ka_max) or ka_max <= 0:
        raise ValueError("Throat ka cap must be finite and positive")
    acoustic_cap = C0 * ka_max / (2 * math.pi * target_f_high)
    if any(r > acoustic_cap * (1 + 1e-12) for r in throat_radii):
        raise ValueError("Throat radius exceeds the acoustic ka cap at target high frequency")
    if profiles is None:
        profiles = DEFAULT_PROFILES

    for name, value, lower, upper in (("mouth radius", mouth_radius, min_mouth_radius, max_mouth_radius), ("length", length, min_length, max_length)):
        if any(bound is not None and (not np.isfinite(bound) or bound <= 0) for bound in (lower, upper)):
            raise ValueError(f"Invalid {name} bounds")
        if lower is not None and upper is not None and lower > upper:
            raise ValueError(f"Invalid {name} bounds: minimum exceeds maximum")
        if value is not None and (not np.isfinite(value) or value <= 0 or (lower is not None and value < lower) or (upper is not None and value > upper)):
            raise ValueError(f"Fixed {name} conflicts with valid dimensions or requested bounds")

    # Mouth radius and length: fixed, or a heuristic range with user limits.
    mouth_range, mouth_sources = _search_range(derive_mouth_radius_range(target_f_low), mouth_radius,
                                               min_mouth_radius, max_mouth_radius)
    length_range, length_sources = _search_range(derive_length_range(target_f_low), length, min_length, max_length)
    if mouth_radius is not None:
        mouth_radii = [mouth_radius]
    else:
        lo, hi = mouth_range
        mouth_radii = np.linspace(lo, hi, num_mouth_radii).tolist() if lo <= hi else []
    if length is not None:
        lengths = [length]
    else:
        lo, hi = length_range
        lengths = np.linspace(lo, hi, num_lengths).tolist() if lo <= hi else []

    sim_range = derive_simulation_freq_range(target_f_low, target_f_high)

    candidates = []
    profile_skipped = False
    idx = 0
    for profile in profiles:
        for r_throat in throat_radii:
            for r_mouth in mouth_radii:
                for horn_length in lengths:
                    if r_mouth <= r_throat:
                        continue
                    if profile_geometry_error(profile, r_throat, r_mouth, horn_length):
                        profile_skipped = True
                        continue  # e.g. a tractrix shorter than its natural length
                    candidate_id = f"auto_{profile[:3]}_{idx:04d}"
                    candidates.append(
                        CandidateGeometry(
                            candidate_id=candidate_id,
                            profile=profile,
                            throat_radius=r_throat,
                            mouth_radius=float(r_mouth),
                            length=float(horn_length),
                        )
                    )
                    idx += 1

    ideal_mouth = derive_mouth_radius(target_f_low)
    derived = DerivedGeometry(
        target_f_low=target_f_low,
        target_f_high=target_f_high,
        ideal_mouth_radius=ideal_mouth,
        mouth_radius_range=mouth_range,
        length_range=length_range,
        sim_freq_range=sim_range,
        candidate_count=len(candidates),
        status="candidates_generated" if candidates else "no_feasible_design",
        reason=_empty_grid_reason(candidates, mouth_radii, lengths, profile_skipped),
        bound_sources={"mouth_radius": mouth_sources, "length": length_sources},
        size_compromise=size_compromise(target_f_low, mouth_range, length_range) if candidates else None,
    )

    return candidates, derived


def generate_fullauto_candidates(
    target_f_low: float,
    target_f_high: float,
    throat_radii: List[float],
    num_mouth_radii: int = 3,
    num_lengths: int = 3,
    profiles: Optional[List[str]] = None,
    ka_max: float = 2 * math.pi,
) -> tuple:
    """Generate geometry candidates from frequency band specification.

    Derives mouth radius and length ranges analytically, then creates
    a focused grid of candidates.

    Args:
        target_f_low: Target low-frequency cutoff (Hz).
        target_f_high: Target high-frequency cutoff (Hz).
        throat_radii: List of throat radii (m) from prescreen.
        num_mouth_radii: Number of mouth radius grid points.
        num_lengths: Number of length grid points.
        profiles: Horn profile types. Defaults to all three.

    Returns:
        Tuple of (candidates, derived_geometry).
    """
    from horn_core.acoustics import validate_band
    validate_band(target_f_low, target_f_high)
    if num_mouth_radii < 1 or num_lengths < 1:
        raise ValueError("Geometry grid sizes must be positive")
    if not throat_radii or any(not np.isfinite(r) or r <= 0 for r in throat_radii):
        raise ValueError("No usable throat radii: no feasible drivers or invalid radii")
    if not np.isfinite(ka_max) or ka_max <= 0:
        raise ValueError("Throat ka cap must be finite and positive")
    acoustic_cap = C0 * ka_max / (2 * math.pi * target_f_high)
    if any(r > acoustic_cap * (1 + 1e-12) for r in throat_radii):
        raise ValueError("Throat radius exceeds the acoustic ka cap at target high frequency")
    if profiles is None:
        profiles = DEFAULT_PROFILES

    mouth_range = derive_mouth_radius_range(target_f_low)
    length_range = derive_length_range(target_f_low)
    sim_range = derive_simulation_freq_range(target_f_low, target_f_high)

    mouth_radii = np.linspace(mouth_range[0], mouth_range[1], num_mouth_radii).tolist()
    lengths = np.linspace(length_range[0], length_range[1], num_lengths).tolist()

    candidates = []
    idx = 0
    for profile in profiles:
        for r_throat in throat_radii:
            for r_mouth in mouth_radii:
                for horn_length in lengths:
                    if r_mouth <= r_throat:
                        continue
                    if profile_geometry_error(profile, r_throat, r_mouth, horn_length):
                        continue  # e.g. a tractrix shorter than its natural length
                    candidate_id = f"fa_{profile[:3]}_{idx:04d}"
                    candidates.append(
                        CandidateGeometry(
                            candidate_id=candidate_id,
                            profile=profile,
                            throat_radius=r_throat,
                            mouth_radius=round(r_mouth, 6),
                            length=round(horn_length, 6),
                        )
                    )
                    idx += 1

    derived = DerivedGeometry(
        target_f_low=target_f_low,
        target_f_high=target_f_high,
        ideal_mouth_radius=derive_mouth_radius(target_f_low),
        mouth_radius_range=mouth_range,
        length_range=length_range,
        sim_freq_range=sim_range,
        candidate_count=len(candidates),
    )

    return candidates, derived


def main():
    """CLI for fullauto geometry derivation."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Derive horn geometry from target frequency band.",
    )
    parser.add_argument(
        "--target-f-low", type=float, required=True, help="Target low frequency (Hz)."
    )
    parser.add_argument(
        "--target-f-high", type=float, required=True, help="Target high frequency (Hz)."
    )
    parser.add_argument(
        "--prescreen-json",
        type=str,
        required=True,
        help="Prescreen result JSON (provides throat_radius_m).",
    )
    parser.add_argument(
        "--mouth-radius", type=float, default=None,
        help="Fixed mouth radius (m). If omitted, derive range from f_low.",
    )
    parser.add_argument(
        "--length", type=float, default=None,
        help="Fixed horn length (m). If omitted, derive range from f_low.",
    )
    parser.add_argument(
        "--num-mouth-radii", type=int, default=3, help="Mouth radius grid points."
    )
    parser.add_argument(
        "--num-lengths", type=int, default=3, help="Length grid points."
    )
    parser.add_argument(
        "--max-length", type=float, default=None,
        help="Maximum length (m). Without --min-length the search runs from 0.4x this value up to it.",
    )
    parser.add_argument(
        "--min-length", type=float, default=None,
        help="Override lower bound of the derived length range (m).",
    )
    parser.add_argument(
        "--max-mouth-radius", type=float, default=None,
        help="Maximum mouth radius (m). Without --min-mouth-radius the search runs from 0.4x this value up to it.",
    )
    parser.add_argument(
        "--min-mouth-radius", type=float, default=None,
        help="Override lower bound of derived mouth radius range (m).",
    )
    parser.add_argument(
        "--output", type=str, default="candidates.csv", help="Output CSV path."
    )
    parser.add_argument(
        "--design-json",
        type=str,
        default="design.json",
        help="Output design summary JSON.",
    )
    parser.add_argument("--throat-radius", type=float, default=None, help="Fix throat radius")
    parser.add_argument("--profiles", type=str, default=None,
                        help="Comma list of profiles to search (default: all seven).")
    args = parser.parse_args()
    profiles = None
    if args.profiles:
        profiles = [name.strip() for name in args.profiles.split(",") if name.strip()]
        unknown = sorted(set(profiles) - set(DEFAULT_PROFILES))
        if not profiles or unknown:
            parser.error(f"Unknown profiles {unknown}; choose from {DEFAULT_PROFILES}")

    prescreen = json.loads(Path(args.prescreen_json).read_text())
    throat_radii = prescreen.get("throat_radii_m") or [prescreen["throat_radius_m"]]
    # The prescreen ka cap limits default throat sampling. A throat the user
    # fixes keeps the wider absolute cap of ka = 2 pi at the target f_high.
    ka_max = prescreen.get("ka_max", 2 * math.pi)
    if args.throat_radius is not None:
        ka_max = max(ka_max, 2 * math.pi)

    candidates, derived = generate_auto_candidates(
        target_f_low=args.target_f_low,
        target_f_high=args.target_f_high,
        throat_radii=[args.throat_radius] if args.throat_radius is not None else throat_radii,
        mouth_radius=args.mouth_radius,
        length=args.length,
        num_mouth_radii=args.num_mouth_radii,
        num_lengths=args.num_lengths,
        max_length=args.max_length,
        min_length=args.min_length,
        max_mouth_radius=args.max_mouth_radius,
        min_mouth_radius=args.min_mouth_radius,
        profiles=profiles,
        ka_max=ka_max,
    )

    write_candidates_csv(candidates, args.output)
    print(f"Generated {len(candidates)} candidates -> {args.output}")

    design = {
        "target_f_low": derived.target_f_low,
        "target_f_high": derived.target_f_high,
        "ideal_mouth_radius": derived.ideal_mouth_radius,
        "throat_radius_range": [min(c.throat_radius for c in candidates), max(c.throat_radius for c in candidates)] if candidates else [],
        "mouth_radius_range": list(derived.mouth_radius_range),
        "length_range": list(derived.length_range),
        "sim_freq_range": list(derived.sim_freq_range),
        "candidate_count": derived.candidate_count,
        "status": derived.status,
        "reason": derived.reason,
        "bound_sources": {**(derived.bound_sources or {}),
                          "throat_radius": ({"low": "user_fixed", "high": "user_fixed"} if args.throat_radius is not None
                                            else {"low": "driver_prescreen", "high": "driver_prescreen"})},
        "throat_ka_max": ka_max,
        "size_compromise": derived.size_compromise,
    }
    Path(args.design_json).write_text(json.dumps(design, indent=2))
    print(f"Design summary -> {args.design_json}")

    print(f"\nDerived geometry for {args.target_f_low}-{args.target_f_high} Hz:")
    print(f"  Ideal mouth radius: {derived.ideal_mouth_radius:.4f} m")
    print(
        f"  Mouth radius range: {derived.mouth_radius_range[0]:.4f} - "
        f"{derived.mouth_radius_range[1]:.4f} m"
    )
    print(
        f"  Length range: {derived.length_range[0]:.4f} - "
        f"{derived.length_range[1]:.4f} m"
    )
    print(
        f"  Simulation freq range: {derived.sim_freq_range[0]:.1f} - "
        f"{derived.sim_freq_range[1]:.1f} Hz"
    )
    if derived.size_compromise:
        c = derived.size_compromise
        print(f"  Size caps are below the sizing rules for {args.target_f_low:g} Hz: horn loading is expected "
              f"to weaken below about {c['expected_f_low_hz']:.0f} Hz ({c['shortfall_octaves']:.2f} octaves; "
              "rule of thumb, see ranked f3_low for predictions)")


if __name__ == "__main__":
    main()
