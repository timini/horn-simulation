# Existing horn validation data

Latest expanded run: [8 September reference validation report](REFERENCE_VALIDATION_REPORT.md).

Research and archive inspection: 7 September 2026. Machine-readable source catalog: [`data/validation/references.json`](../data/validation/references.json).

There is substantial existing DIY measurement data. It can support several distinct checks; a successful download is not itself validation of our predictions.

## Downloaded numeric data

| Source | Acquired data | Useful comparison | Conditions still to establish |
| --- | --- | --- | --- |
| [Bill Waslo — SEOS15](https://libinst.com/SEOS/SEOS15/SEOS15.zip) | 26 FRD curves: horizontal and vertical, 0–90° at 7.5° increments | Relative directivity, frequency-response shape | Exact driver, baffle, drive voltage, microphone distance and processing for this particular archive. Phase columns contain zeros; do not treat these as measured phase. |
| [Waslo — SEOS12 / Delta12LFA](https://libinst.com/SEOS/SEOSD12LFAAug2012.zip) | Three FRD and two ZMA files; crossover workbooks also present | Electrical impedance and system/driver response, once the crossover status is established | Raw versus filtered response, voltage, distance, exact baffle/interface. ZMA is **electrical**, not throat acoustic impedance. |
| [Waslo — SEOS12 / Designer12](https://libinst.com/SEOS/SEOSDES12Aug2012.zip) | Three FRD and two ZMA files | Comparison system and woofer data | Its tweeter FRD and DNA360 ZMA are byte-identical to those in the previous archive. They are not independent horn measurements. |

The importer successfully processed **383 curves, 378 unique files from six archives**: 44 DIY response/electrical-impedance curves (39 unique files), plus 299 measured and 40 independently simulated acoustic-impedance curves from the pipe benchmark below. Repeated measurements are not independent assemblies. It preserves original levels, phase values and frequencies. Original archive and member hashes are recorded. It ignores the workbooks/macros and images. These non-axisymmetric SEOS geometries cannot be replaced by a circular horn with a similar mouth area and then presented as a matched physical validation.

Reproduce the acquisition and import:

```sh
python -m horn_analysis.reference_data \
  --catalog data/validation/references.json \
  --output-dir results/references
```

The catalog pins SHA-256 hashes; changed downloads fail instead of silently replacing the reference. Previously fetched archives can be supplied with `--cache-dir`. The local preparation output is `results/validation-81/imported/inventory.json`. Source archives are intentionally not vendored: public download availability does not establish an open redistribution license.

## Matched-dimension pipe benchmark

[Ernoult and colleagues, Zenodo v2](https://zenodo.org/records/20024938), licensed CC BY 4.0, supplies **299 measured impedance files and 40 independent simulation files**. The downloaded 126.5 MB archive matches the publisher's MD5 and our pinned SHA-256. There are five experimenters, multiple specimens and repeat measurements; these are retained separately to expose reproducibility and material differences. Some numeric files have no extension and are intentionally supported by the importer.

- Cylinder: 180 mm long, 14 mm inner diameter. Closed ends, or finite circular flanges of 2 mm (brass) and 7 mm (wood/ABS).
- Cone: 180 mm long, diameter 10–22.6 mm. Closed or 2.7 mm flange measurements.
- Infinitely thin unflanged terminations appear in **simulations only**.
- Measured impedance is dimensionless, normalized by entrance `rho*c/S`. Simulated impedance is acoustic `Pa*s/m³`. Neither is electrical ohms. The importer preserves these distinct units in column names and metadata.
- The authors already transpose measurements to dry air at 25°C. Their inspected analysis source specifies `c = 346.28592 m/s`, `rho = 1.184490 kg/m³` and the positive-time harmonic convention. Source URL and checksum are pinned; downloaded scripts were inspected as text, not executed.

This dataset can test propagation, termination and repeatability independently of an uncertain loudspeaker motor. It does not validate driver choice or far-field horn output. The production solver now supports opt-in thermoviscous side-wall losses, finite flanges and rigid closed outlets. Roughness, porosity, material compliance and cap losses remain unmodeled. The earlier lossless comparison below is retained as a historical diagnostic.

The two additional Waslo archives are [SEOS12 / Delta12A](https://libinst.com/SEOS/SEOS12+Eminence%20pro%20woofers/SEOS12%20Delta%2012A.zip) and [SEOS12 / DeltaLite2512](https://libinst.com/SEOS/SEOS12+Eminence%20pro%20woofers/SEOS12%20DeltaLite2512.zip), each with two FRD and two ZMA files. Their tweeter data include duplicates; filenames suggest a 30-inch response distance, but complete setup and filtering must still be established.

## First measured comparison: a diagnostic failure worth keeping

`scripts/benchmark_measured_impedance.py` now runs the production FEM on the specified 180 mm × 14 mm cylinder at two mesh sizes and compares the shared lossless/unflanged model with **36 brass-pipe measurements from four operators**, over a fixed 110–3900 Hz interval. It fits no parameters. Frequency is mapped using the published sound speed to preserve `kL`, and specific solver impedance is divided by `rho*c` to match the measured normalization.

The preliminary magnitude criterion (median error ≤2 dB and 95th-percentile error ≤4 dB) passes for only **16/36 curves**. Across measurements, the median error is typically 0.82 dB, but the typical 95th-percentile error is 4.24 dB (range 2.79–6.59 dB). The first resonance peak is overpredicted by 25.1–32.2 dB across these measurements; an average-error pass would hide that defect. The model does not reproduce wall losses or the measured 2 mm flange, so even those 16 magnitude matches are **not physical-validation passes**.

FEM versus the same analytical model improves with mesh refinement: the 95th-percentile scaled complex error `|Z_FEM - Z_TMM|/(1+|Z_TMM|)`, using normalized impedance, falls from 4.98% at 4 mm mesh to 1.31% at 2 mm. Its maximum is still 10.54% at 2 mm near sensitive frequencies. This separates numerical consistency from the missing physical effects; it does not attribute all measurement differences to a single cause.

Reproduce inside the solver container with the repository mounted and its source packages available:

```sh
python3 scripts/benchmark_measured_impedance.py \
  --reference-dir results/validation-81/imported/ernoult-pipe-impedance-v2 \
  --output-dir results/validation-81/measured-pipe-comparison
```

The output includes the acoustic CAD, both raw FEM runs, frequency-transposed normalized exports, dense analytical prediction, per-measurement errors, provenance and a comparison plot. A compact frozen result is checked in as [`pipe_diagnostic_baseline.json`](../data/validation/pipe_diagnostic_baseline.json), including source hashes. Keep this failure as a baseline when adding loss/termination physics; do not tune thresholds until it passes. Reserve other materials and the conical measurements for held-out checks.

## Additional primary sources inspected

| Source | Available evidence | Intended use / limitation |
| --- | --- | --- |
| [MTG Great Waveguide Shootout](https://www.mtg-designs.com/tips-tricks-tests/waveguide-shootout) | Large collection of named commercial and DIY horns, with off-axis measurements and setup notes | Strong comparative source. Individual pages must be checked; early distortion data were withdrawn by the author as inaccurate. |
| [MTG ST260](https://www.mtg-designs.com/tips-tricks-tests/waveguide-shootout/1in-bolt-horns/ath-st260) | Celestion CDX1-1731, 1 m, 4.5 ms gate, 1/48-octave plots | Most promising first circular DIY case when combined with exact ST260 geometry. Numeric response files not acquired. |
| [MTG Tritonia-M](https://www.mtg-designs.com/tips-tricks-tests/waveguide-shootout/1in-bolt-horns/tritonia-m) | Same horn with CDX1-1731 and B&C DH450-8; 1 m, 4.5 ms, 1/48 octave | Useful driver substitution comparison. Non-axisymmetric; measured plots, not downloaded FRD. |
| [ATH ST260](https://at-horns.eu/ST260.html) | CAD, ABEC project and five linked DIY measurement reports; CC BY-NC-SA 4.0 | Downloaded ABEC archive contains model inputs, not computed response curves. Includes free-standing and baffled variants; do not mix them. |
| [ATH Gen2 measurements](https://at-horns.eu/gen2m.html) | Many named driver, throat-adapter and horn-body combinations | Excellent shape/adapter comparisons. Author explicitly states uncalibrated SPL and varying drive voltage; not an absolute-output or sensitivity reference. |
| [ATH A460D](https://at-horns.eu/A460D.html) | Prototype measurements with five named drivers and different adapters | Same calibration caveats as Gen2. Exact adapter is part of the geometry. |
| [JW Sound Paraflex 2x12 CRAM](https://www.jwsound.live/paraflex-2x12-measurements) | Downloadable REW `.mdat`, named drivers, outdoor 10 m ground-plane setup, environmental and calibration notes | Rich response/compression data. No input-normalized sensitivity supplied. Folded multi-path enclosure is outside the current front-loaded horn model. Archive link identified, not imported by the FRD/ZMA importer. |
| [Ampslab/JBL 2414HC–Dayton H6512 archive](https://audiokarma.org/forums/threads/excellent-crossover-design-simulation-for-jbl-2414hc-on-jbl-6x12-clone-waveguide.1019895/) | Published FRD/ZMA attachment supplied by the measurer | Archive retrieval returned HTTP 403. The discussion explicitly distinguishes H6512 from JBL clones; geometry substitution would invalidate comparison. |
| [Rasetshwane & Neely, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4617734/) | Calibrated reflectance measurements of known conical, exponential and parabolic horns | Promising input-impedance reference. Sealed mouth cap and calibrated source must be reproduced; not open-mouth sound output. Raw numeric samples not acquired. |
| [Kolbrek exponential midrange](https://kolbrek.hornspeakersystems.info/index.php/horns/exponential-midrange-1) | Author-built rectangular exponential horn, Altec 288B, outdoor response and impedance | Useful later rectangular-horn reference; recover full dimensions/setup and digitize curves with uncertainty. |

## Turning these sources into acceptance tests

1. First reproduce a specified source assembly, including the driver exit, adapter, baffle/rollback, and measurement reference plane. Keep unsupported geometries in the research catalog until the model can represent them.
2. Import raw FRD/ZMA/REW exports where supplied. For plots, retain the figure identity, axis calibration and digitization uncertainty; never label traced points as raw measurements.
3. Compare only the frequency interval resolved by the measurement gate. Never treat the mere existence of a low-frequency point in a gated plot as valid low-frequency evidence.
4. Compare matching quantities. Electrical ZMA requires a motor/electrical prediction; normalized polars require the same angle normalization; absolute SPL requires known voltage and calibration.
5. Freeze setup, predictions and tolerances before comparison. Reserve another assembly/adapter as a held-out check if any parameters are fitted.
6. Use `horn_analysis.validation` for a response comparison. It rejects uncalibrated absolute comparisons, includes prediction samples to expose notches, and never converts a curve-agreement result into a full physical-validation claim.

The existing repository files named `conical_horn_hornresp.csv` and `exponential_horn_ijert.csv` contain internally generated Webster predictions for published geometries. Their headers identify this. They are **not measured curves or an independently executed Hornresp reference** and must not count toward the physical gate.


### Additional fully coupled assembly lead

[Sa and Park (2014)](https://doi.org/10.5050/KSNVE.2014.24.7.537) describes two conical assemblies with driver/chamber parameters and measured impedance/sensitivity figures. The source catalog records the known dimensions and author-copy link. Text inspection confirms a richer chamber/rear-volume and frequency-dependent motor model than the current direct-coupling approximation. Curves have not been extracted and absolute normalization still needs verification. This is a documented validation candidate, not an acquired numeric dataset or a passing physical comparison.
