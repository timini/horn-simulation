#!/usr/bin/env python3
"""Compare an exported Hornresp results file with this project's output.

Usage, after entering records/<id>.txt in Hornresp and exporting its chart data:

    .venv/bin/python validation/hornresp/compare_hornresp.py --horn conical --export path/to/export.txt
    .venv/bin/python validation/hornresp/compare_hornresp.py --horn example_800_1600 --export export.txt --spl

The export is read as a table whose header line contains a frequency column
("Freq ..."). Columns are found by name: Ra and Xa (throat acoustic impedance
normalised to rho c / S1) and SPL. Tab, comma, semicolon and runs of two or
more spaces are accepted as separators. If your Hornresp version names the
columns differently, pass --freq-col, --ra-col, --xa-col or --spl-col.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# Frozen with the records, before any Hornresp file was available.
# Impedance: both codes are plane-wave 1D models with a baffled piston mouth.
# SPL: Hornresp (plane-wave mouth) against our modal_baffled FEM, which is
# expected to differ by the nonuniform mouth velocity.
# spl_halfspace was added in the review fix pass, still before any Hornresp
# data: our own 1D model built the way Hornresp builds SPL (TMM horn, ideal
# transformer, mouth volume velocity into half space). Our 1D half-space
# curve itself differs from the FEM on-axis curve by median 1.06 dB and max
# 1.41 dB over 800 to 1600 Hz, so the "spl" check above is expected to fail
# even if Hornresp matches our 1D model exactly. It is kept as frozen.
TOLERANCES = {
    "impedance": {"db_median": 0.2, "db_p95": 1.0, "db_max": 3.0, "phase_p95_deg": 5.0, "peak_freq_rel": 0.01,
                  "band_hz": [20.0, 20000.0]},
    "spl": {"db_median": 0.5, "db_max": 2.0, "band_hz": [800.0, 1600.0]},
    "spl_halfspace": {"db_median": 0.5, "db_max": 2.0, "band_hz": [800.0, 1600.0]},
}


def _split(line, sep):
    if sep == "spaces":
        return [c.strip() for c in re.split(r"\s{2,}", line.strip()) if c.strip()]
    return [c.strip() for c in line.strip().split(sep)]


def read_export(path, overrides=None):
    """Return {column name: float array} from a Hornresp text export."""
    lines = Path(path).read_text(errors="replace").splitlines()
    for i, line in enumerate(lines):
        if re.search(r"\bfreq", line, re.IGNORECASE):
            header_index = i
            break
    else:
        raise ValueError("No header line with a frequency column found")
    header = lines[header_index]
    sep = next((s for s in ("\t", ",", ";") if s in header), "spaces")
    names = _split(header, sep)
    rows = []
    for line in lines[header_index + 1:]:
        if not line.strip():
            continue
        cells = _split(line, sep) if sep != "spaces" else line.split()
        try:
            values = [float(c.replace("D", "E")) for c in cells[:len(names)]]
        except ValueError:
            continue
        if len(values) == len(names):
            rows.append(values)
    if not rows:
        raise ValueError("No numeric rows after the header")
    data = np.array(rows)
    table = {name: data[:, j] for j, name in enumerate(names)}
    return table, names


def pick(names, wanted, override=None):
    if override is not None:
        return names[int(override)] if str(override).isdigit() else override
    for name in names:
        token = name.split()[0].strip("()").lower() if name.split() else ""
        if token == wanted.lower():
            return name
    for name in names:
        if name.lower().startswith(wanted.lower()):
            return name
    raise ValueError(f"No '{wanted}' column in {names}; pass the column explicitly")


def stats(f, ours, ref, tol):
    lo, hi = tol["band_hz"]
    m = (f >= lo) & (f <= hi)
    db = np.abs(20 * np.log10(np.abs(ours[m]) / np.abs(ref[m])))
    out = {"band_hz": [float(f[m][0]), float(f[m][-1])], "points": int(m.sum()),
           "db_median": float(np.median(db)), "db_p95": float(np.percentile(db, 95)), "db_max": float(db.max())}
    if np.iscomplexobj(ours):
        ph = np.abs(np.degrees(np.angle(ours[m] / ref[m])))
        out.update(phase_median_deg=float(np.median(ph)), phase_p95_deg=float(np.percentile(ph, 95)),
                   phase_max_deg=float(ph.max()))
    out["checks"] = {k: out[k] <= v for k, v in tol.items() if k in out and k != "band_hz"}
    return out


def peak_errors(f, ours, ref, tol):
    from scipy.signal import find_peaks
    lo, hi = tol["band_hz"]
    m = (f >= lo) & (f <= hi)
    pr, _ = find_peaks(20 * np.log10(np.abs(ref[m])), prominence=1.0)
    po, _ = find_peaks(20 * np.log10(np.abs(ours[m])), prominence=1.0)
    fr, fo = f[m][pr], f[m][po]
    # One to one, closest pairs first, so one of our peaks is never used twice.
    pairs = sorted((abs(q - x), i, j) for i, x in enumerate(fr) for j, q in enumerate(fo))
    matched, used = {}, set()
    for _, i, j in pairs:
        if i not in matched and j not in used:
            matched[i] = j
            used.add(j)
    errors = [float((fo[matched[i]] - x) / x) if i in matched else None for i, x in enumerate(fr)]
    return {"reference_peaks_hz": fr.tolist(), "our_peaks_hz": fo.tolist(), "peak_rel_errors": errors,
            "peak_count_match": len(fr) == len(fo),
            "peak_freq_rel_max": max((abs(e) for e in errors if e is not None), default=0.0)}


def compare_impedance(horn_id, table, names, args, air):
    sys.path.insert(0, str(HERE))
    from hornresp_records import our_output
    sheets = {h["id"]: h for h in json.loads((HERE / "records/horns.json").read_text())["horns"]}
    f = table[pick(names, "freq", args.freq_col)]
    ref = table[pick(names, "Ra", args.ra_col)] + 1j * table[pick(names, "Xa", args.xa_col)]
    keep = f > 0
    f, ref = f[keep], ref[keep]
    _, ours = our_output(sheets[horn_id], air=air, freqs=f)
    tol = TOLERANCES["impedance"]
    result = stats(f, ours, ref, tol)
    result.update(peak_errors(f, ours, ref, tol))
    result["checks"]["peak_freq_rel"] = result["peak_freq_rel_max"] <= tol["peak_freq_rel"]
    result["checks"]["peak_count"] = result["peak_count_match"]
    return result


def _spl_against(f_ours, spl_ours, table, names, args, tol):
    f = table[pick(names, "freq", args.freq_col)]
    spl = table[pick(names, "SPL", args.spl_col)]
    lo, hi = max(f.min(), f_ours.min()), min(f.max(), f_ours.max())
    m = (f_ours >= lo) & (f_ours <= hi)
    ref = np.interp(np.log(f_ours[m]), np.log(f), spl)
    p_ours, p_ref = 10 ** (spl_ours[m] / 20), 10 ** (ref / 20)
    return stats(f_ours[m], p_ours, p_ref, tol)


def compare_spl(table, names, args):
    """Hornresp SPL against our modal FEM on-axis response (frozen check)."""
    import csv
    rows = list(csv.DictReader(open(ROOT / "examples/800-1600-modal/response.csv")))
    f_ours = np.array([float(r["frequency"]) for r in rows])
    spl_ours = np.array([float(r["spl"]) for r in rows])
    return _spl_against(f_ours, spl_ours, table, names, args, TOLERANCES["spl"])


def compare_spl_halfspace(table, names, args):
    """Hornresp SPL against our 1D model built the same way (half-space point source)."""
    data = np.loadtxt(HERE / "reference/example_800_1600_spl_halfspace.csv", delimiter=",", skiprows=1)
    return _spl_against(data[:, 0], data[:, 1], table, names, args, TOLERANCES["spl_halfspace"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--horn", required=True, help="Record id, e.g. conical or example_800_1600")
    parser.add_argument("--export", required=True, type=Path, help="Hornresp exported text file")
    parser.add_argument("--spl", action="store_true",
                        help="Also compare SPL with examples/800-1600-modal/response.csv and with "
                             "reference/example_800_1600_spl_halfspace.csv")
    parser.add_argument("--c", type=float, help="Speed of sound used in Hornresp (default from records/horns.json)")
    parser.add_argument("--rho", type=float, help="Air density used in Hornresp (default from records/horns.json)")
    parser.add_argument("--output", type=Path, help="Write the comparison JSON here")
    for col in ("freq", "ra", "xa", "spl"):
        parser.add_argument(f"--{col}-col")
    args = parser.parse_args(argv)
    from horn_core.duct import AirProperties
    air_json = json.loads((HERE / "records/horns.json").read_text())["air"]
    air = AirProperties(c=args.c or air_json["c"], rho=args.rho or air_json["rho"])
    table, names = read_export(args.export)
    result = {"horn": args.horn, "export": str(args.export), "air": {"c": air.c, "rho": air.rho},
              "tolerances": TOLERANCES, "impedance": compare_impedance(args.horn, table, names, args, air)}
    if args.spl:
        result["spl"] = compare_spl(table, names, args)
        result["spl_halfspace"] = compare_spl_halfspace(table, names, args)
    result["pass"] = all(all(v["checks"].values()) for k, v in result.items() if isinstance(v, dict) and "checks" in v)
    text = json.dumps(result, indent=1)
    if args.output:
        args.output.write_text(text + "\n")
    print(text)
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
