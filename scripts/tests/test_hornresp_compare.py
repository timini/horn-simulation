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
