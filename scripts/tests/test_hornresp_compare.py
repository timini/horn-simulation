"""The Hornresp comparison must read an export and catch a wrong curve."""
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2] / "validation/hornresp"


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


records = load("hornresp_records")
compare = load("compare_hornresp")


def write_export(path, f, z, spl=None, sep="\t"):
    header = ["Freq (hertz)", "Ra (norm)", "Xa (norm)", "Za (norm)"] + (["SPL (dB)"] if spl is not None else [])
    lines = ["Hornresp export (synthetic test file)", "", sep.join(header)]
    for i, (fi, zi) in enumerate(zip(f, z)):
        cells = [f"{fi:.4f}", f"{zi.real:.6E}", f"{zi.imag:.6E}", f"{abs(zi):.6E}"]
        if spl is not None:
            cells.append(f"{spl[i]:.4f}")
        lines.append(sep.join(cells))
    path.write_text("\n".join(lines) + "\n")


def test_own_curve_passes_and_a_shifted_curve_fails(tmp_path):
    sheet = json.loads((HERE / "records/horns.json").read_text())["horns"][0]
    assert sheet["id"] == "conical"
    f = np.geomspace(20, 20000, 400)
    _, z = records.our_output(sheet, freqs=f)
    good = tmp_path / "good.txt"
    write_export(good, f, z)
    assert compare.main(["--horn", "conical", "--export", str(good)]) == 0
    bad = tmp_path / "bad.txt"
    write_export(bad, f * 1.03, z, sep="  ")       # resonances 3 % off
    assert compare.main(["--horn", "conical", "--export", str(bad)]) == 1


def test_columns_are_found_by_name():
    names = ["Freq (hertz)", "Ra (norm)", "Xa (norm)", "SPL (dB)"]
    assert compare.pick(names, "Ra") == "Ra (norm)"
    assert compare.pick(names, "SPL") == "SPL (dB)"
    assert compare.pick(names, "freq") == "Freq (hertz)"


def test_segmented_os_is_close_to_the_profile():
    sheets = {h["id"]: h for h in json.loads((HERE / "records/horns.json").read_text())["horns"]}
    assert sheets["os"]["max_radius_error"] < 0.02
    assert len(sheets["os"]["segments"]) == 4


def test_spl_against_our_halfspace_curve_passes_and_the_fem_check_is_expected_to_fail(tmp_path):
    data = np.loadtxt(HERE / "reference/example_800_1600_spl_halfspace.csv", delimiter=",", skiprows=1)
    sheets = {h["id"]: h for h in json.loads((HERE / "records/horns.json").read_text())["horns"]}
    f = data[:, 0]
    _, z = records.our_output(sheets["example_800_1600"], freqs=f)
    export = tmp_path / "ex.txt"
    write_export(export, f, z, spl=data[:, 1])
    out = tmp_path / "ex.json"
    assert compare.main(["--horn", "example_800_1600", "--export", str(export), "--spl", "--output", str(out)]) == 1
    result = json.loads(out.read_text())
    assert all(result["spl_halfspace"]["checks"].values())
    assert all(result["impedance"]["checks"].values())
    # Our own 1D curve differs from the FEM on-axis response by about 1 dB median.
    assert not result["spl"]["checks"]["db_median"]


def test_peaks_are_matched_one_to_one():
    f = np.linspace(100, 1000, 901)
    ref = 1 + np.exp(-((f - 300) / 10) ** 2) * 10 + np.exp(-((f - 330) / 10) ** 2) * 10
    ours = 1 + np.exp(-((f - 310) / 10) ** 2) * 10
    out = compare.peak_errors(f, ours, ref, {"band_hz": [100, 1000]})
    assert len(out["reference_peaks_hz"]) == 2 and len(out["our_peaks_hz"]) == 1
    assert out["peak_rel_errors"].count(None) == 1
    assert not out["peak_count_match"]


def test_band_limits_are_not_checked_as_a_tolerance():
    f = np.geomspace(25, 18000, 300)
    z = np.ones_like(f) * (1 + 0.5j)
    out = compare.stats(f, z, z, {"db_median": 0.2, "band_hz": [20.0, 20000.0]})
    assert out["checks"] == {"db_median": True}
