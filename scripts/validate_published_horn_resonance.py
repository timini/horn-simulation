#!/usr/bin/env python3
"""Limited independent check of Sa & Park's explicitly reported sixth resonance.

The paper supplies one numerical measured peak in its discussion of Fig. 13.
No figure tracing, measured amplitude, motor parameters or fitted correction
filters are used. Unknown air and flange conditions prevent certification.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from horn_core.duct import DEFAULT_AIR
from horn_core.webster import compute_horn_transfer_tmm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    out = parser.parse_args().output_dir
    out.mkdir(parents=True, exist_ok=False)
    frequency = np.linspace(110, 2000, 18901)
    result = {
        "source_url": "https://doi.org/10.5050/KSNVE.2014.24.7.537",
        "source_location": "Section 2.2, discussion of Fig. 13, page 545",
        "reference_kind": "author_reported_measured_scalar",
        "measured_sixth_resonance_hz": 1712.,
        "geometry": {"length_m": .535, "throat_diameter_m": .018, "mouth_diameter_m": .080},
        "assumed_air": asdict(DEFAULT_AIR), "assumed_flange_width_m": 0.,
        "loss_model": "boundary_layer", "radiation_model": "finite_flange",
        "frequency_step_hz": .1, "fitted_parameters": [],
        "limitations": ["Temperature, humidity and exact lip geometry are not established.",
            "The source reports 44.1 kHz / 4096-sample blocks; peak estimation uncertainty is not supplied.",
            "One resonance frequency does not validate impedance magnitude, bandwidth, output SPL or driver coupling."],
        "acceptance_gate": None, "physical_assembly_validated": False,
    }
    (out/"protocol.json").write_text(json.dumps(result, indent=2)+"\n")
    result["predictions"] = {}
    for count in (400, 800):
        transfer = compute_horn_transfer_tmm(frequency, lambda z: .009+.031*z/.535,
            .535, .009, .04, n_segments=count, loss_model="boundary_layer",
            radiation_model="finite_flange", flange_width=0., air=DEFAULT_AIR)
        impedance = transfer["z_real"]+1j*transfer["z_imag"]
        peaks = find_peaks(np.abs(impedance))[0]
        if len(peaks) != 6:
            raise ValueError("Expected six resonances in the fixed analysis band")
        sixth = float(frequency[peaks[5]])
        csv = out/f"prediction-{count}.csv"
        pd.DataFrame({"frequency": frequency, "specific_z_real": impedance.real,
                      "specific_z_imag": impedance.imag}).to_csv(csv, index=False)
        result["predictions"][str(count)] = {
            "resonances_hz": frequency[peaks].tolist(),
            "sixth_resonance_error_hz": sixth-1712.,
            "sixth_resonance_error_percent": 100*(sixth-1712.)/1712.,
            "prediction_sha256": hashlib.sha256(csv.read_bytes()).hexdigest()}
    result["source_code_sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (Path(__file__), Path("packages/horn-core/src/horn_core/webster.py"),
         Path("packages/horn-core/src/horn_core/duct.py"))}
    (out/"comparison.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result["predictions"], indent=2))


if __name__ == "__main__":
    main()
