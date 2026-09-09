"""KPI extraction from horn frequency response data.

Computes key performance indicators from a CSV with (frequency, spl) columns:
- f3_low / f3_high: -3 dB crossings or explicitly flagged sweep bounds
- bandwidth_hz / bandwidth_octaves: contiguous peak lobe, possibly a lower bound
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
from horn_core.acoustics import validate_response


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
    f3_low_is_bound: bool = False
    f3_high_is_bound: bool = False
    bandwidth_is_lower_bound: bool = False

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
    freq = validate_response(freq, spl)
    spl = np.asarray(spl, dtype=float)
    # Peak
    peak_idx = np.argmax(spl)
    peak_spl = float(spl[peak_idx])
    peak_freq = float(freq[peak_idx])

    # -3 dB threshold
    threshold = peak_spl - 3.0

    # Interpolate for root-finding
    spl_interp = interp1d(freq, spl, kind="linear", fill_value="extrapolate")

    # Walk the original knots outward from the peak. Resampling before
    # crossing detection can hide a narrow notch or join disconnected lobes.
    values = spl-threshold
    f3_low, low_bound = _peak_lobe_edge(freq, values, peak_idx, -1)
    f3_high, high_bound = _peak_lobe_edge(freq, values, peak_idx, 1)

    # Derived KPIs
    bandwidth_hz = None
    bandwidth_octaves = None
    passband_ripple = None
    avg_sensitivity = None

    if f3_low is not None and f3_high is not None:
        bandwidth_hz = f3_high - f3_low
        bandwidth_octaves = np.log2(f3_high / f3_low) if f3_low > 0 else None

        # Keep original knots and interpolated edges for exact linear extrema.
        # Uniform resampling alone can miss narrow resonances or notches.
        knots = np.r_[f3_low, freq[(freq > f3_low) & (freq < f3_high)], f3_high]
        passband_ripple = float(np.ptp(spl_interp(knots)))
        freq_uniform = np.geomspace(f3_low, f3_high, max(200, len(freq) * 2))
        spl_passband = spl_interp(freq_uniform)
        avg_sensitivity = float(np.mean(spl_passband))

    return HornKPI(
        peak_spl_db=peak_spl,
        peak_frequency_hz=peak_freq,
        f3_low_hz=f3_low,
        f3_high_hz=f3_high,
        bandwidth_hz=bandwidth_hz,
        bandwidth_octaves=float(bandwidth_octaves) if bandwidth_octaves is not None else None,
        passband_ripple_db=passband_ripple,
        average_sensitivity_db=avg_sensitivity,
        f3_low_is_bound=low_bound,
        f3_high_is_bound=high_bound,
        bandwidth_is_lower_bound=low_bound or high_bound,
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


def _peak_lobe_edge(freq, values, peak, direction):
    """Find the first threshold crossing outward from the sampled peak."""
    index = peak
    while 0 <= index+direction < len(freq):
        neighbour = index+direction
        if values[neighbour] <= 0:
            fraction = values[index]/(values[index]-values[neighbour])
            return float(freq[index]+fraction*(freq[neighbour]-freq[index])), False
        index = neighbour
    return float(freq[index]), bool(values[index] > 0)


def format_kpi(kpi, key, precision='.0f', missing='N/A'):
    """Render sweep limits explicitly in every report using KPI metadata."""
    if isinstance(kpi, HornKPI):
        kpi = kpi.to_dict()
    value = kpi.get(key)
    if value is None:
        return missing
    prefix = ''
    if key == 'f3_low_hz' and kpi.get('f3_low_is_bound'):
        prefix = '≤'
    elif key == 'f3_high_hz' and kpi.get('f3_high_is_bound'):
        prefix = '≥'
    elif key.startswith('bandwidth_') and kpi.get('bandwidth_is_lower_bound'):
        prefix = '≥'
    return prefix+format(value, precision)


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
