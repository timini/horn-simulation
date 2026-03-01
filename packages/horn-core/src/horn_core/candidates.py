"""Geometry candidate generation for auto-select pipeline.

Generates a grid of candidate horn geometries (profiles x throat radii x
mouth radii x lengths) constrained by the user's target frequency range
and geometry limits.
"""

import csv
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import List, Optional

import numpy as np


# Standard 1" compression driver throat radii
DEFAULT_THROAT_RADII = [0.0125, 0.0175, 0.025]  # metres

DEFAULT_PROFILES = ["conical", "exponential", "hyperbolic", "tractrix", "os", "lecleach", "cd"]


@dataclass
class CandidateGeometry:
    """A single candidate horn geometry for screening."""
    candidate_id: str
    profile: str
    throat_radius: float  # m
    mouth_radius: float   # m
    length: float          # m


def generate_candidates(
    target_f_low: float,
    target_f_high: float,
    max_length: float,
    max_mouth_radius: float,
    profiles: Optional[List[str]] = None,
    throat_radii: Optional[List[float]] = None,
    num_lengths: int = 4,
    num_mouth_radii: int = 4,
) -> List[CandidateGeometry]:
    """Generate a grid of candidate horn geometries.

    Args:
        target_f_low: Target low-frequency cutoff (Hz), used for context.
        target_f_high: Target high-frequency cutoff (Hz), used for context.
        max_length: Maximum horn length (m).
        max_mouth_radius: Maximum mouth radius (m).
        profiles: List of profile names. Defaults to all three.
        throat_radii: List of throat radii (m). Defaults to standard 1" driver sizes.
        num_lengths: Number of length grid points.
        num_mouth_radii: Number of mouth radius grid points.

    Returns:
        List of CandidateGeometry objects.
    """
    if profiles is None:
        profiles = DEFAULT_PROFILES
    if throat_radii is None:
        throat_radii = DEFAULT_THROAT_RADII

    # Generate evenly-spaced grids within constraints
    # Mouth radius: from a reasonable minimum to max
    min_mouth = max(max_mouth_radius * 0.25, 0.03)
    mouth_radii = np.linspace(min_mouth, max_mouth_radius, num_mouth_radii).tolist()

    # Length: from a reasonable minimum to max
    min_length = max(max_length * 0.3, 0.1)
    lengths = np.linspace(min_length, max_length, num_lengths).tolist()

    candidates = []
    idx = 0
    for profile, r_throat, r_mouth, horn_length in product(
        profiles, throat_radii, mouth_radii, lengths
    ):
        # Skip invalid geometries: mouth must be larger than throat
        if r_mouth <= r_throat:
            continue

        candidate_id = f"{profile[:3]}_{idx:04d}"
        candidates.append(CandidateGeometry(
            candidate_id=candidate_id,
            profile=profile,
            throat_radius=r_throat,
            mouth_radius=r_mouth,
            length=horn_length,
        ))
        idx += 1

    return candidates


def write_candidates_csv(candidates: List[CandidateGeometry], output_path: str) -> Path:
    """Write candidates to a CSV file for Nextflow consumption."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "candidate_id", "profile", "throat_radius", "mouth_radius", "length",
        ])
        for c in candidates:
            writer.writerow([
                c.candidate_id, c.profile,
                c.throat_radius, c.mouth_radius, c.length,
            ])
    return path


def main():
    """CLI for candidate generation."""
    import argparse

    parser = argparse.ArgumentParser(description="Generate candidate horn geometries.")
    parser.add_argument("--target-f-low", type=float, default=500, help="Target low freq (Hz).")
    parser.add_argument("--target-f-high", type=float, default=4000, help="Target high freq (Hz).")
    parser.add_argument("--max-length", type=float, default=0.5, help="Max horn length (m).")
    parser.add_argument("--max-mouth-radius", type=float, default=0.2, help="Max mouth radius (m).")
    parser.add_argument("--profiles", type=str, nargs="*", default=None, help="Profile names.")
    parser.add_argument("--output", type=str, default="candidates.csv", help="Output CSV path.")
    args = parser.parse_args()

    candidates = generate_candidates(
        target_f_low=args.target_f_low,
        target_f_high=args.target_f_high,
        max_length=args.max_length,
        max_mouth_radius=args.max_mouth_radius,
        profiles=args.profiles,
    )

    path = write_candidates_csv(candidates, args.output)
    print(f"Generated {len(candidates)} candidates -> {path}")


if __name__ == "__main__":
    main()
