"""Analytical geometry derivation for fullauto horn design.

Given only a target frequency band, derives optimal horn geometry
(mouth radius, length) using horn acoustics formulas, then generates
a focused grid of candidates for FEM evaluation.

Key formulas:
  - Mouth radius from low-freq target: r_mouth = c0 / (2*pi*f_low)
    (mouth circumference = wavelength at f_low)
  - Length range: quarter-wave to half-wave at f_low
    L_min = c0 / (4*f_low), L_max = c0 / (2*f_low)
  - Simulation freq range: extend +/-0.5 octave beyond target band
    to capture rolloff for accurate f3 detection
"""

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np

from horn_core.candidates import CandidateGeometry, write_candidates_csv

# Speed of sound in air at ~20C
C0 = 343.0

DEFAULT_PROFILES = ["conical", "exponential", "hyperbolic"]


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


def derive_simulation_freq_range(
    f_low: float, f_high: float
) -> tuple:
    """Extend target band by +/-0.5 octave for rolloff capture.

    This ensures the f3 points (where SPL drops 3dB) can be
    accurately detected even if they fall outside the target band.
    """
    sim_min = f_low / (2 ** 0.5)  # -0.5 octave
    sim_max = f_high * (2 ** 0.5)  # +0.5 octave
    return (sim_min, sim_max)


def generate_fullauto_candidates(
    target_f_low: float,
    target_f_high: float,
    throat_radii: List[float],
    num_mouth_radii: int = 3,
    num_lengths: int = 3,
    profiles: Optional[List[str]] = None,
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
        "--num-mouth-radii", type=int, default=3, help="Mouth radius grid points."
    )
    parser.add_argument(
        "--num-lengths", type=int, default=3, help="Length grid points."
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
    args = parser.parse_args()

    prescreen = json.loads(Path(args.prescreen_json).read_text())
    throat_radius = prescreen["throat_radius_m"]

    candidates, derived = generate_fullauto_candidates(
        target_f_low=args.target_f_low,
        target_f_high=args.target_f_high,
        throat_radii=[throat_radius],
        num_mouth_radii=args.num_mouth_radii,
        num_lengths=args.num_lengths,
    )

    write_candidates_csv(candidates, args.output)
    print(f"Generated {len(candidates)} candidates -> {args.output}")

    design = {
        "target_f_low": derived.target_f_low,
        "target_f_high": derived.target_f_high,
        "ideal_mouth_radius": derived.ideal_mouth_radius,
        "mouth_radius_range": list(derived.mouth_radius_range),
        "length_range": list(derived.length_range),
        "sim_freq_range": list(derived.sim_freq_range),
        "candidate_count": derived.candidate_count,
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


if __name__ == "__main__":
    main()
