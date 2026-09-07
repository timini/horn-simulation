"""Compare independent CSV responses without fitting an arbitrary level offset."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from horn_core.acoustics import validate_band, validate_response


def compare_curves(reference, prediction, low, high, *, relative=False):
    validate_band(low, high)
    for frame in (reference, prediction):
        validate_response(frame.frequency, frame.spl)
        if frame.frequency.iloc[0] > low or frame.frequency.iloc[-1] < high:
            raise ValueError("Both datasets must bracket the entire comparison band")
    # A fixed grid weights equal log-frequency bandwidth equally, independently
    # of how densely either source sampled a particular part of the band.
    grid = np.geomspace(low, high, 10001)
    a = np.interp(np.log(grid), np.log(reference.frequency), reference.spl)
    b = np.interp(np.log(grid), np.log(prediction.frequency), prediction.spl)
    if relative:
        # Explicit shape-only comparison, never satisfies the absolute gate.
        a, b = a-a[0], b-b[0]
    error = np.abs(a-b)
    # Preserve exact extrema at every original knot, without letting those
    # knots change the bandwidth-weighted percentile gates.
    knots = np.unique(np.r_[low, high, reference.frequency, prediction.frequency])
    knots = knots[(knots >= low) & (knots <= high)]
    knot_error = np.interp(np.log(knots), np.log(reference.frequency), reference.spl) - np.interp(np.log(knots), np.log(prediction.frequency), prediction.spl)
    if relative:
        knot_error -= knot_error[0]
    median, p95 = float(np.median(error)), float(np.percentile(error,95))
    return {"mode": "relative_shape" if relative else "absolute_level",
            "median_absolute_error_db":median, "p95_absolute_error_db":p95,
            "max_absolute_error_db":float(np.abs(knot_error).max()), "samples":len(grid),
            "comparison_grid": "10001 points uniform in log frequency",
            "comparison_band_hz":[low,high], "level_gate_passed": not relative and median<=2 and p95<=4,
            "shape_gate_passed": (median<=2 and p95<=4) if relative else None,
            "physical_validation_passed": False,
            "note":"Curve agreement alone does not validate geometry, driver, calibration or ranking."}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',required=True)
    p.add_argument('--prediction',required=True)
    p.add_argument('--metadata',required=True,help='Reference provenance and measurement conditions JSON')
    p.add_argument('--output',required=True)
    p.add_argument('--low',type=float,required=True)
    p.add_argument('--high',type=float,required=True)
    p.add_argument('--relative',action='store_true')
    a=p.parse_args();meta=json.loads(Path(a.metadata).read_text())
    for field in ('source_url','reference_kind','quantity','conditions','limitations'):
        if field not in meta: p.error(f'Missing reference metadata: {field}')
    if not a.relative and not meta.get('absolute_level_calibrated',False):
        p.error('Reference does not establish absolute calibration; use --relative for shape only')
    result=compare_curves(pd.read_csv(a.reference),pd.read_csv(a.prediction),a.low,a.high,relative=a.relative)
    result['reference_metadata']=meta
    result['sha256']={key:hashlib.sha256(Path(path).read_bytes()).hexdigest() for key,path in [('reference',a.reference),('prediction',a.prediction)]}
    Path(a.output).write_text(json.dumps(result,indent=2))
    # In --relative mode success means shape agreement only.
    gate = 'shape_gate_passed' if a.relative else 'level_gate_passed'
    return 0 if result[gate] else 1


if __name__=='__main__':
    raise SystemExit(main())
