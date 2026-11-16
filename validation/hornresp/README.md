# Hornresp comparison kit

Hornresp is a Windows program. It is not installed here and Wine is not available on the machine that produced these files, so no Hornresp result has been compared yet. This folder lets you run Hornresp yourself and check its output against ours.

## Contents

- `records/<id>.txt`: the values to type into the Hornresp input window for each horn, in Hornresp units (cm2, cm, Hz). One sheet per profile at throat radius 25 mm, mouth radius 150 mm, length 300 mm (the tractrix uses its natural length, 223.8 mm), plus `example_800_1600.txt`, the 800 to 1600 Hz worked example horn with the 18Sound 6NMB420 driver.
- `records/horns.json`: the same values in machine-readable form, with the air properties used.
- `reference/<id>_tmm.csv`: our throat impedance (Ra, Xa normalised to rho c / S1) for the geometry exactly as entered in Hornresp, baffled piston mouth, 2000-segment TMM.
- `hornresp_records.py`: regenerates the sheets and reference curves.
- `compare_hornresp.py`: compares an exported Hornresp file with our output.

The sheets are not Hornresp record files. The record file layout differs between Hornresp versions and could not be checked without the program, so the values are listed for manual entry instead.

## Geometry

| Horn | Hornresp flare | Match |
| --- | --- | --- |
| conical | Con | exact |
| exponential | Exp | exact |
| hyperbolic (hypex T = 0) | Hyp, T = 0 | exact |
| tractrix | Tra | exact if the length Hornresp derives from S1, S2 and F12 equals 22.38 cm |
| os | 4 segments (Exp, Exp, Exp, Con) | 1.30 % largest radius error |
| scaled_tractrix | 4 segments | 1.36 % largest radius error |
| cd | Exp then Con | exact |

For OS and the scaled tractrix our reference curve uses the segmented geometry, so the comparison checks the solvers on the same shape. The segment approximation is a separate, stated error.

## Steps

1. Set Hornresp's air to c = 344 m/s and rho = 1.205 kg/m3 (or note what your copy uses and pass `--c` and `--rho`).
2. Enter one sheet. Use Ang = 2 x Pi so the mouth is in an infinite baffle, which matches our `flanged_piston` load.
3. For the impedance sheets, read the throat acoustic impedance chart and export the chart data as text.
4. For `example_800_1600`, also export the SPL chart (2.83 V, 1 m).
5. Run, from the repository root:

```sh
.venv/bin/python validation/hornresp/compare_hornresp.py --horn conical --export ~/Hornresp/Export/conical.txt
.venv/bin/python validation/hornresp/compare_hornresp.py --horn example_800_1600 --export ex.txt --spl
```

The script finds the header line that contains a frequency column and picks the Ra, Xa and SPL columns by name. If your export names them differently, pass `--freq-col`, `--ra-col`, `--xa-col` or `--spl-col` with the column name or index. It exits 1 if any frozen tolerance fails.

## Tolerances (frozen before any Hornresp data)

| Check | Limit |
| --- | --- |
| Throat impedance magnitude, 20 Hz to 20 kHz | median 0.2 dB, p95 1.0 dB, max 3.0 dB |
| Throat impedance phase | p95 5 degrees |
| Impedance peak frequencies | 1 % |
| SPL against `examples/800-1600-modal/response.csv`, 800 to 1600 Hz | median 0.5 dB, max 2.0 dB |

Both codes are plane-wave models for the impedance, so near agreement is expected. The SPL reference is our modal FEM with a nonuniform mouth velocity, while Hornresp assumes a plane-wave piston mouth; part of any SPL difference comes from that.

## Driver assumptions for the example

Our coupling uses the published Mms as the whole moving mass, no rear chamber (open back) and an ideal area transformer from Sd to the throat. In Hornresp enter Mms as Mmd, Vrc = 0, Vtc = 0. Cms and Rms are derived from fs, Mms and Qms exactly as the project does. Hornresp may add its own air load to Mmd, which would shift the low-frequency response slightly.
