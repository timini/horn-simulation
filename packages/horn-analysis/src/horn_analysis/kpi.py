"""KPI extraction from horn frequency response data.

The ``spl`` column is the area-averaged RMS pressure on the horn mouth
plane, in dB re 20 uPa. It is a mouth-plane level, not a 1 m sensitivity;
see docs/model-assumptions.md. Every KPI below inherits that definition.

Computes key performance indicators from a CSV with (frequency, spl) columns:
- f3_low / f3_high: -3 dB cutoff frequencies
- bandwidth_hz / bandwidth_octaves: usable bandwidth
- passband_ripple_db: max - min SPL within the -3 dB band
- average_level_db: mean mouth-plane level in the passband
- peak_spl_db / peak_frequency_hz: maximum mouth-plane level
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
from scipy.signal import savgol_filter


@dataclass
class HornKPI:
    """Key performance indicators for a horn frequency response.

    All dB values are mouth-plane levels (dB re 20 uPa), not 1 m
    sensitivities. See docs/model-assumptions.md.
    """
    peak_spl_db: float
    peak_frequency_hz: float
    f3_low_hz: Optional[float]
    f3_high_hz: Optional[float]
    bandwidth_hz: Optional[float]
    bandwidth_octaves: Optional[float]
    passband_ripple_db: Optional[float]
    average_level_db: Optional[float]

    def to_dict(self) -> dict:
        return asdict(self)


def extract_kpis_from_arrays(freq: np.ndarray, spl: np.ndarray) -> HornKPI:
    """Extract KPIs from frequency/SPL arrays directly.

    Args:
        freq: Array of frequency values in Hz.
        spl: Array of mouth-plane SPL values (dB re 20 uPa).

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
    if f3_low is None and spl_minus_threshold(freq[0]) >= 0:
        # Response is above -3 dB at the low measurement boundary
        f3_low = float(freq[0])

    # Find f3_high: search from peak to highest freq
    f3_high = _find_crossing(spl_minus_threshold, freq[peak_idx], freq[-1], direction="falling")
    if f3_high is None and spl_minus_threshold(freq[-1]) >= 0:
        # Response is above -3 dB at the high measurement boundary
        f3_high = float(freq[-1])

    # Derived KPIs
    bandwidth_hz = None
    bandwidth_octaves = None
    passband_ripple = None
    avg_level = None

    if f3_low is not None and f3_high is not None:
        bandwidth_hz = f3_high - f3_low
        bandwidth_octaves = np.log2(f3_high / f3_low) if f3_low > 0 else None

        # Passband: resample onto uniform log grid and smooth to remove
        # band-stitching artifacts before computing ripple/level
        n_passband = max(200, len(freq) * 2)
        freq_uniform = np.geomspace(f3_low, f3_high, n_passband)
        spl_passband = np.array([float(spl_interp(f)) for f in freq_uniform])

        # Light Savitzky-Golay smoothing to suppress band-boundary glitches
        # Window must be odd and < n_passband; 11 points is ~5% of 200
        sg_window = min(11, n_passband if n_passband % 2 == 1 else n_passband - 1)
        if sg_window >= 5:
            spl_passband = savgol_filter(spl_passband, sg_window, polyorder=3)

        passband_ripple = float(np.max(spl_passband) - np.min(spl_passband))
        avg_level = float(np.mean(spl_passband))

    return HornKPI(
        peak_spl_db=peak_spl,
        peak_frequency_hz=peak_freq,
        f3_low_hz=f3_low,
        f3_high_hz=f3_high,
        bandwidth_hz=bandwidth_hz,
        bandwidth_octaves=float(bandwidth_octaves) if bandwidth_octaves is not None else None,
        passband_ripple_db=passband_ripple,
        average_level_db=avg_level,
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
