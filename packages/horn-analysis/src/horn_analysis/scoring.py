"""Scoring and ranking for driver-horn selection.

Scores each driver-horn combination on bandwidth coverage, passband
ripple, and average sensitivity, then ranks candidates.
"""

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Optional

import numpy as np

from horn_analysis.kpi import HornKPI


@dataclass
class TargetSpec:
    """User-specified target frequency range and geometry constraints."""
    f_low_hz: float
    f_high_hz: float
    max_length_m: Optional[float] = None
    max_mouth_radius_m: Optional[float] = None


@dataclass
class SelectionScore:
    """Composite score for a single driver-horn combination."""
    driver_id: str
    horn_label: str
    bandwidth_coverage: float       # 0-1: fraction of target range covered
    passband_ripple_db: float       # ripple in dB (lower is better)
    avg_sensitivity_db: float       # mean SPL in passband
    composite_score: float          # weighted combination, 0-1

    def to_dict(self) -> dict:
        return asdict(self)


def compute_selection_score(
    kpi: HornKPI,
    target: TargetSpec,
    driver_id: str = "",
    horn_label: str = "",
    weights: Optional[dict] = None,
) -> SelectionScore:
    """Score a single driver-horn combination against a target specification.

    Args:
        kpi: Extracted KPIs for the coupled response.
        target: User's target frequency range.
        driver_id: Identifier for the driver.
        horn_label: Identifier for the horn geometry.
        weights: Optional dict with keys ``bandwidth``, ``ripple``,
                 ``sensitivity``.  Defaults to 0.50/0.25/0.25.

    Returns:
        SelectionScore with composite score in [0, 1].
    """
    if weights is None:
        weights = {"bandwidth": 0.50, "ripple": 0.25, "sensitivity": 0.25}

    target_width = target.f_high_hz - target.f_low_hz

    # --- Bandwidth coverage ---
    if kpi.f3_low_hz is not None and kpi.f3_high_hz is not None:
        overlap_low = max(kpi.f3_low_hz, target.f_low_hz)
        overlap_high = min(kpi.f3_high_hz, target.f_high_hz)
        overlap = max(0.0, overlap_high - overlap_low)
        bandwidth_coverage = overlap / target_width if target_width > 0 else 0.0
    else:
        bandwidth_coverage = 0.0

    bandwidth_coverage = min(bandwidth_coverage, 1.0)

    # --- Ripple score ---
    ripple_db = kpi.passband_ripple_db if kpi.passband_ripple_db is not None else 6.0
    ripple_score = max(0.0, 1.0 - ripple_db / 6.0)

    # --- Sensitivity score ---
    avg_sens = kpi.average_sensitivity_db if kpi.average_sensitivity_db is not None else 80.0
    sensitivity_score = np.clip((avg_sens - 80.0) / 40.0, 0.0, 1.0)

    # --- Composite ---
    composite = (
        weights["bandwidth"] * bandwidth_coverage
        + weights["ripple"] * ripple_score
        + weights["sensitivity"] * float(sensitivity_score)
    )

    return SelectionScore(
        driver_id=driver_id,
        horn_label=horn_label,
        bandwidth_coverage=bandwidth_coverage,
        passband_ripple_db=ripple_db,
        avg_sensitivity_db=avg_sens,
        composite_score=composite,
    )


def rank_candidates(
    scores: List[SelectionScore],
    top_n: int = 5,
) -> List[SelectionScore]:
    """Rank candidates by composite score and return the top N."""
    return sorted(scores, key=lambda s: s.composite_score, reverse=True)[:top_n]


def main():
    """CLI for scoring a single KPI JSON against a target spec."""
    parser = argparse.ArgumentParser(description="Score a horn-driver combination.")
    parser.add_argument("--kpi-json", required=True, help="KPI JSON file.")
    parser.add_argument("--target-f-low", type=float, required=True, help="Target low freq (Hz).")
    parser.add_argument("--target-f-high", type=float, required=True, help="Target high freq (Hz).")
    parser.add_argument("--driver-id", type=str, default="", help="Driver identifier.")
    parser.add_argument("--horn-label", type=str, default="", help="Horn geometry label.")
    parser.add_argument("--output", type=str, default=None, help="Output JSON file.")
    args = parser.parse_args()

    kpi_data = json.loads(Path(args.kpi_json).read_text())
    kpi = HornKPI(**kpi_data)
    target = TargetSpec(f_low_hz=args.target_f_low, f_high_hz=args.target_f_high)

    score = compute_selection_score(
        kpi, target,
        driver_id=args.driver_id,
        horn_label=args.horn_label,
    )

    result = json.dumps(score.to_dict(), indent=2)
    if args.output:
        Path(args.output).write_text(result)
        print(f"Score written to {args.output}")
    else:
        print(result)


if __name__ == "__main__":
    main()
