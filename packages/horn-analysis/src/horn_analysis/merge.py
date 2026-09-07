"""Fail-closed merge of the expected frequency bands of one candidate."""
import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd
from horn_core.acoustics import validate_band, validate_response, inlet_area_from_frame


def merge_bands(paths, *, num_bands, min_freq, max_freq, points_per_band, output):
    validate_band(min_freq, max_freq)
    if num_bands < 1 or points_per_band < 2:
        raise ValueError("Need at least one band and two points per band")
    indexed = {}
    for path in paths:
        match = re.search(r"_(\d+)\.csv$", str(path))
        if not match or int(match[1]) in indexed:
            raise ValueError("Missing or duplicate band identity")
        indexed[int(match[1])] = path
    if set(indexed) != set(range(num_bands)):
        raise ValueError(f"Incomplete simulation: expected bands 0..{num_bands-1}, got {sorted(indexed)}")
    frames = []
    for index, path in sorted(indexed.items()):
        df = pd.read_csv(path)
        for col in ("frequency", "spl", "z_real", "z_imag"):
            if col not in df:
                raise ValueError(f"Missing solver column {col}: {path}")
        validate_response(df.frequency, *[df[c] for c in df.select_dtypes(include="number")])
        width = (max_freq - min_freq) / num_bands
        expected = np.geomspace(min_freq+index*width, min_freq+(index+1)*width, points_per_band)
        if len(df) != len(expected) or not np.allclose(df.frequency, expected, rtol=1e-9, atol=1e-8):
            raise ValueError(f"Incomplete or incorrect frequency grid in band {index}")
        if "schema_version" not in df or not (df.schema_version == 2).all():
            raise ValueError("Old solver results cannot be merged into a schema-2 run")
        for col in ("inlet_area_m2", "mouth_area_m2", "mouth_u_real", "mouth_u_imag", "mouth_p_real", "mouth_p_imag", "radiation_model", "phasor_convention", "bc_mode"):
            if col not in df or df[col].isna().any():
                raise ValueError(f"Missing acoustic contract column {col}")
        frames.append(df)
    joined = pd.concat(frames, ignore_index=True).sort_values("frequency")
    inlet_area_from_frame(joined)
    for col in ("mouth_area_m2",):
        if (joined[col] <= 0).any() or not np.allclose(joined[col], joined[col].iloc[0], rtol=1e-6):
            raise ValueError(f"Inconsistent {col} across bands")
    for col in ("radiation_model", "phasor_convention", "bc_mode"):
        if joined[col].nunique() != 1:
            raise ValueError(f"Inconsistent {col} across bands")
    for col in ("loss_model", "flange_width_m", "air_c_m_s", "air_rho_kg_m3",
                "air_gamma", "air_viscosity_pa_s", "air_conductivity_w_m_k",
                "air_heat_capacity_j_kg_k", "element_degree"):
        if col in joined and (joined[col].isna().any() or joined[col].nunique() != 1):
            raise ValueError(f"Inconsistent {col} across bands")
    # An SPL match alone does not constrain the complex transfers used by the
    # motor. Check impedance and aperture pressure/volume velocity as well.
    from horn_core.duct import DEFAULT_AIR
    characteristic = (float(joined.air_c_m_s.iloc[0])*float(joined.air_rho_kg_m3.iloc[0])
                      if "air_c_m_s" in joined and "air_rho_kg_m3" in joined
                      else DEFAULT_AIR.c*DEFAULT_AIR.rho)
    scales = {"z": characteristic, "mouth_p": 1.,
              "mouth_u": float(joined.mouth_area_m2.iloc[0])/characteristic}
    for _, rows in joined.groupby("frequency"):
        if len(rows) < 2:
            continue
        if rows.spl.max()-rows.spl.min() > 0.5:
            raise ValueError("Band boundary pressure mismatch exceeds 0.5 dB; refine mesh")
        for prefix, scale in scales.items():
            values = rows[prefix+"_real"].to_numpy()+1j*rows[prefix+"_imag"].to_numpy()
            # 5% relative complex difference plus 1% of a fixed unit-pressure
            # reference scale near zeros; phase jumps cannot hide behind SPL.
            tolerance = .05*np.maximum(np.abs(values),abs(values[0])) + .01*scale
            if np.any(np.abs(values-values[0]) > tolerance):
                raise ValueError(f"Band boundary complex {prefix} mismatch; refine mesh")
    joined = joined.drop_duplicates("frequency", keep="first")
    validate_response(joined.frequency, joined.spl)
    joined.to_csv(output, index=False)
    return joined


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--files", nargs="+", required=True)
    p.add_argument("--num-bands", type=int, required=True)
    p.add_argument("--min-freq", type=float, required=True)
    p.add_argument("--max-freq", type=float, required=True)
    p.add_argument("--points-per-band", type=int, required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    merge_bands(a.files, num_bands=a.num_bands, min_freq=a.min_freq, max_freq=a.max_freq,
                points_per_band=a.points_per_band, output=a.output)


if __name__ == "__main__":
    main()
