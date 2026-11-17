"""Fail-closed merge of the expected frequency bands of one candidate."""
import argparse
from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd
from horn_core.acoustics import validate_band, validate_response, inlet_area_from_frame
from horn_core.frequency_grid import DEFAULT_MIN_POINTS_PER_OCTAVE, band_edges, points_per_band as grid_points_per_band


SPL_JUMP_LIMIT_DB = 0.5
COMPLEX_RELATIVE_TOLERANCE = .05
COMPLEX_SCALE_TOLERANCE = .01


class BandBoundaryMismatch(ValueError):
    """Neighbouring bands disagree at their shared frequency.

    The bands of one candidate are solved on separate meshes, so a jump at a
    shared frequency means at least one mesh is not resolved well enough.
    ``quantity`` is ``spl`` (``jump`` in dB) or a complex column prefix
    (``jump`` relative to the larger magnitude). ``excess_ratio`` is the jump
    divided by the allowed difference, so anything above 1 failed.
    """
    reason = "band_boundary_mismatch"

    def __init__(self, worst, failing_checks):
        self.quantity, self.frequency_hz, self.jump, self.excess_ratio = worst
        self.failing_checks = failing_checks
        if self.quantity == "spl":
            what = (f"Band boundary pressure mismatch exceeds {SPL_JUMP_LIMIT_DB} dB: "
                    f"{self.jump:.2f} dB jump")
        else:
            what = (f"Band boundary complex {self.quantity} mismatch: {self.jump:.1%} relative jump, "
                    f"{self.excess_ratio:.2f} times the allowed difference")
        super().__init__(f"{what} at {self.frequency_hz:.1f} Hz "
                         f"({failing_checks} failing boundary checks); refine mesh")

    def summary(self):
        return {"quantity": self.quantity, "frequency_hz": self.frequency_hz,
                "jump": self.jump, "jump_unit": "dB" if self.quantity == "spl" else "relative",
                "excess_ratio": self.excess_ratio, "failing_checks": self.failing_checks}


def boundary_rejection(error, *, candidate_id, output=None):
    """Rejection record for a candidate whose bands do not join."""
    import json
    record = {"candidate_id": candidate_id, "status": "rejected",
              "rejection_reasons": [error.reason], "rejection_detail": str(error),
              "band_boundary_mismatch": error.summary(), "bands": []}
    if output is not None:
        Path(output).write_text(json.dumps(record, indent=2))
    return record


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
    edges = band_edges(min_freq, max_freq, num_bands)
    frames = []
    for index, path in sorted(indexed.items()):
        df = pd.read_csv(path)
        for col in ("frequency", "spl", "z_real", "z_imag"):
            if col not in df:
                raise ValueError(f"Missing solver column {col}: {path}")
        validate_response(df.frequency, *[df[c] for c in df.select_dtypes(include="number")])
        expected = np.geomspace(edges[index], edges[index+1], points_per_band)
        if len(df) != len(expected) or not np.allclose(df.frequency, expected, rtol=1e-9, atol=1e-8):
            raise ValueError(f"Incomplete or incorrect frequency grid in band {index}")
        if "schema_version" not in df or not (df.schema_version == 2).all():
            raise ValueError("Old solver results cannot be merged into a schema-2 run")
        for col in ("inlet_area_m2", "mouth_area_m2", "mouth_u_real", "mouth_u_imag", "mouth_p_real", "mouth_p_imag", "radiation_model", "phasor_convention", "bc_mode"):
            if col not in df or df[col].isna().any():
                raise ValueError(f"Missing acoustic contract column {col}")
        if (df.radiation_model=='modal_baffled').any():
            if not (df.radiation_model=='modal_baffled').all() or 'modal_mode_count' not in df or not (df.modal_mode_count==16).all():
                raise ValueError('Inconsistent modal aperture mode count')
            for col in [f'modal_v_{n}_{part}' for n in range(16) for part in ('real','imag')]:
                if col not in df or df[col].isna().any():
                    raise ValueError('Incomplete modal aperture coefficients')
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
    if joined.radiation_model.iloc[0]=='modal_baffled':
        scales.update({f'modal_v_{n}':1/characteristic for n in range(16)})
    failures = []
    for frequency, rows in joined.groupby("frequency"):
        if len(rows) < 2:
            continue
        spl_jump = float(rows.spl.max()-rows.spl.min())
        if spl_jump > SPL_JUMP_LIMIT_DB:
            failures.append(("spl", float(frequency), spl_jump, spl_jump/SPL_JUMP_LIMIT_DB))
        for prefix, scale in scales.items():
            values = rows[prefix+"_real"].to_numpy()+1j*rows[prefix+"_imag"].to_numpy()
            # 5% relative complex difference plus 1% of a fixed unit-pressure
            # reference scale near zeros; phase jumps cannot hide behind SPL.
            difference = np.abs(values-values[0])
            tolerance = COMPLEX_RELATIVE_TOLERANCE*np.maximum(np.abs(values),abs(values[0])) + COMPLEX_SCALE_TOLERANCE*scale
            if np.any(difference > tolerance):
                worst = int(np.argmax(difference/tolerance))
                failures.append((prefix, float(frequency),
                                 float(difference[worst]/max(np.abs(values).max(), 1e-300)),
                                 float(difference[worst]/tolerance[worst])))
    if failures:
        raise BandBoundaryMismatch(max(failures, key=lambda f: f[3]), len(failures))
    joined = joined.drop_duplicates("frequency", keep="first")
    validate_response(joined.frequency, joined.spl)
    joined.to_csv(output, index=False)
    return joined


def merge_rejections(paths, *, candidate_id, output):
    """Combine per-band solver rejections into one record for the candidate.

    One rejected band means the candidate has no complete response, so the
    whole candidate is rejected; the reasons and band details are kept.
    """
    import json
    records = [json.loads(Path(path).read_text()) for path in paths]
    if not records or any(r.get("status") != "rejected" or not r.get("rejection_reasons") for r in records):
        raise ValueError("Rejection files must each name at least one reason")
    reasons = sorted({reason for r in records for reason in r["rejection_reasons"]})
    record = {"candidate_id": candidate_id, "status": "rejected", "rejection_reasons": reasons,
              "rejection_detail": "; ".join(r.get("rejection_detail", "") for r in records),
              "bands": records}
    Path(output).write_text(json.dumps(record, indent=2))
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--files", nargs="+")
    p.add_argument("--rejections", nargs="+", help="Per-band solver rejection JSON files")
    p.add_argument("--candidate-id", help="Candidate label for --rejections or --rejection-file")
    p.add_argument("--rejection-file", help="With --files: write a band_boundary_mismatch rejection here instead of failing")
    p.add_argument("--num-bands", type=int)
    p.add_argument("--min-freq", type=float)
    p.add_argument("--max-freq", type=float)
    p.add_argument("--num-intervals", type=int)
    p.add_argument("--min-points-per-octave", type=float, default=DEFAULT_MIN_POINTS_PER_OCTAVE)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    if a.rejections:
        if not a.candidate_id:
            p.error("--rejections needs --candidate-id")
        merge_rejections(a.rejections, candidate_id=a.candidate_id, output=a.output)
        return
    if not a.files or None in (a.num_bands, a.min_freq, a.max_freq, a.num_intervals):
        p.error("--files needs --num-bands, --min-freq, --max-freq and --num-intervals")
    if a.rejection_file and not a.candidate_id:
        p.error("--rejection-file needs --candidate-id")
    points = grid_points_per_band(a.min_freq, a.max_freq, a.num_bands, a.num_intervals, a.min_points_per_octave)
    try:
        merge_bands(a.files, num_bands=a.num_bands, min_freq=a.min_freq, max_freq=a.max_freq,
                    points_per_band=points, output=a.output)
    except BandBoundaryMismatch as error:
        # In auto mode one unresolved candidate is rejected; the run continues.
        if not a.rejection_file:
            raise
        boundary_rejection(error, candidate_id=a.candidate_id, output=a.rejection_file)
        print(f"Rejected {a.candidate_id}: {error}", file=sys.stderr)


if __name__ == "__main__":
    main()
