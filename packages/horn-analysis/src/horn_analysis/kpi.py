"""KPI extraction from horn frequency response data.

Computes key performance indicators from a CSV with (frequency, spl) columns:
- f3_low / f3_high: -3 dB cutoff frequencies
- bandwidth_hz / bandwidth_octaves: usable bandwidth
- passband_ripple_db: max - min SPL within the -3 dB band
- average_sensitivity_db: mean SPL in the passband
- peak_spl_db / peak_frequency_hz: maximum SPL point
"""

import json
import argparse
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import brentq


@dataclass
class HornKPI:
    """Key performance indicators for a horn frequency response."""
    peak_spl_db: float
    peak_frequency_hz: float
    f3_low_hz: Optional[float]
    f3_high_hz: Optional[float]
    bandwidth_hz: Optional[float]
    bandwidth_octaves: Optional[float]
    passband_ripple_db: Optional[float]
    average_sensitivity_db: Optional[float]

    def to_dict(self) -> dict:
        return asdict(self)


def extract_kpis_from_arrays(freq: np.ndarray, spl: np.ndarray) -> HornKPI:
    """Extract KPIs from frequency/SPL arrays directly.

    Args:
        freq: Array of frequency values in Hz.
        spl: Array of SPL values in dB.

    Returns:
        HornKPI dataclass with computed metrics.
    """
    # Peak
    peak_idx = np.argmax(spl)
    peak_spl = float(spl[peak_idx])
    peak_freq = float(freq[peak_idx])

    # -3 dB threshold
    threshold = peak_spl - 3.0

    # Interpolate for root-finding
    spl_interp = interp1d(freq, spl, kind="linear", fill_value="extrapolate")

    def spl_minus_threshold(f):
        return float(spl_interp(f)) - threshold

    # Find f3_low: search from lowest freq up to peak
    f3_low = _find_crossing(spl_minus_threshold, freq[0], freq[peak_idx], direction="rising")

    # Find f3_high: search from peak to highest freq
    f3_high = _find_crossing(spl_minus_threshold, freq[peak_idx], freq[-1], direction="falling")

    # Derived KPIs
    bandwidth_hz = None
    bandwidth_octaves = None
    passband_ripple = None
    avg_sensitivity = None

    if f3_low is not None and f3_high is not None:
        bandwidth_hz = f3_high - f3_low
        bandwidth_octaves = np.log2(f3_high / f3_low) if f3_low > 0 else None

        # Passband: frequencies within [f3_low, f3_high]
        mask = (freq >= f3_low) & (freq <= f3_high)
        if np.any(mask):
            passband_spl = spl[mask]
            passband_ripple = float(np.max(passband_spl) - np.min(passband_spl))
            avg_sensitivity = float(np.mean(passband_spl))

    return HornKPI(
        peak_spl_db=peak_spl,
        peak_frequency_hz=peak_freq,
        f3_low_hz=f3_low,
        f3_high_hz=f3_high,
        bandwidth_hz=bandwidth_hz,
        bandwidth_octaves=float(bandwidth_octaves) if bandwidth_octaves is not None else None,
        passband_ripple_db=passband_ripple,
        average_sensitivity_db=avg_sensitivity,
    )


def extract_kpis(csv_path: str) -> HornKPI:
    """Extract KPIs from a frequency response CSV file.

    Args:
        csv_path: Path to CSV with 'frequency' and 'spl' columns.

    Returns:
        HornKPI dataclass with computed metrics.
    """
    df = pd.read_csv(csv_path)
    return extract_kpis_from_arrays(df["frequency"].values, df["spl"].values)


def _find_crossing(func, f_start, f_end, direction="rising"):
    """Find where func crosses zero between f_start and f_end.

    For 'rising', we look for a transition from negative to positive.
    For 'falling', we look for a transition from positive to negative.

    Returns the crossing frequency, or None if not found.
    """
    # Sample the function at multiple points to find a bracket
    num_samples = 200
    f_samples = np.linspace(f_start, f_end, num_samples)
    vals = np.array([func(f) for f in f_samples])

    if direction == "rising":
        # Find first interval where val goes from < 0 to >= 0
        for i in range(len(vals) - 1):
            if vals[i] < 0 and vals[i + 1] >= 0:
                try:
                    return float(brentq(func, f_samples[i], f_samples[i + 1]))
                except ValueError:
                    continue
    else:  # falling
        # Find last interval where val goes from >= 0 to < 0
        for i in range(len(vals) - 2, -1, -1):
            if vals[i] >= 0 and vals[i + 1] < 0:
                try:
                    return float(brentq(func, f_samples[i], f_samples[i + 1]))
                except ValueError:
                    continue

    return None


def main():
    """CLI for KPI extraction."""
    parser = argparse.ArgumentParser(description="Extract KPIs from horn frequency response CSV.")
    parser.add_argument("csv_file", type=str, help="Input CSV with frequency,spl columns.")
    parser.add_argument("--output", type=str, default=None, help="Output JSON file (default: stdout).")
    args = parser.parse_args()

    kpis = extract_kpis(args.csv_file)
    result = json.dumps(kpis.to_dict(), indent=2)

    if args.output:
        Path(args.output).write_text(result)
        print(f"KPIs written to {args.output}")
    else:
        print(result)


if __name__ == "__main__":
    main()
