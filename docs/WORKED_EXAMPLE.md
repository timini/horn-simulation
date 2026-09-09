# Worked example: an 800–1600 Hz horn candidate

The pipeline can now generate and refine a candidate from an operating band, rank the driver/horn combinations, and produce response data, dimensions, acoustic STEP geometry and a report. This example demonstrates that software workflow. It is **an experimental prediction**, not a physically qualified or manufacturing-ready assembly.

## Reproduce the search

From a clean checkout, install the README prerequisites and run:

```sh
just build
just run-auto --target_f_low 800 --target_f_high 1600 \
  --drivers_db data/drivers-curated --lem_top_n 3 --refinement_budget 2 \
  --num_bands 2 --num_intervals 101 --mesh_size 0.01
just latest-run
```

The last command returns the completed run directory. Open `outputs/auto/report/auto_report.html` inside it. The report is self-contained. The default drive is 2.83 V RMS and the observation point is 1 m on axis, using the lossless horn and uniform circular piston in an infinite baffle. The full resolved specification is stored in the run.

The bounded example deliberately evaluates three screened geometries and two refinement proposals. A larger search is possible, but this example does not prove a global optimum or that physical differences between close scores are meaningful.

## Generated candidate

The completed run and a fresh run after the report fixes selected the following leading numerical candidate:

| Item | Value |
| --- | --- |
| Driver | 18Sound 6NMB420, curated manufacturer-parameter record |
| Profile | Hyperbolic |
| Throat diameter | 64.996 mm |
| Mouth diameter | 177.418 mm |
| Acoustic length | 117.906 mm |
| Target-band predicted ripple | Approximately 0.99 dB |
| Predicted mean target-band output | Approximately 104.94 dB SPL at 1 m, 2.83 V RMS |
| Evidence status | Insufficient evidence for a physical recommendation |

The STEP file is `outputs/auto/refinement/refine_hyperbolic_0001.step`. Its matching raw FEM data are `refine_hyperbolic_0001_results.csv` in the same directory. Driver-coupled SPL data are in `outputs/auto/report/coupled_01_18sound-6nmb420_refine_hyperbolic_0001.csv`. Full rankings retain near-tie information, rejection reasons, model assumptions and evidence gaps.

The throat/mouth dimensions describe the air passage. They do not specify horn wall thickness, driver mounting, front chamber, adapter/phase-plug geometry or a separately characterized rear load. The curated record lacks a verified interface and separate moving-mass/rear-load decomposition. An ideal piston coupled at the throat cannot establish how this real cone driver behaves through a physical compression chamber. Do not fabricate this as a supposedly validated assembly from the air-volume STEP alone.

## Resolution and reporting

The [resolution study](CANDIDATE_RESOLUTION.md) freezes the geometry, driver and band, then separately increases mesh, loft and frequency resolution. Numerical convergence does not improve the missing physical evidence or establish the order of nearly tied designs.

Reports distinguish observed −3 dB crossings from sweep bounds. For this example the upper crossing is outside the simulated range: `≥2263 Hz` is a bound, not an exact cutoff. [Response metric definitions](REPORT_METRICS.md) explain the difference between target-band acceptance and secondary peak-lobe KPIs.

## Resume and preserve evidence

Resume with the run path returned by the launcher:

```sh
python3 scripts/run_pipeline.py --run-dir results/<run-id> -resume
```

The source, inputs and container identities must still match. Successful completion and a cached resume have been exercised for the original example. Preserve the manifest and all outputs together; do not edit old evidence to match a newer checkout. See [run management](RUN_MANAGEMENT.md) for lookup and archival cleanup.

The remaining step to a dependable physical recommendation is the [reference/interface/radiation sequence](SINGLE_HORN_ROADMAP.md), tracked in [#81](https://github.com/timini/horn-simulation/issues/81). Neither this successful workflow nor its numerical comparisons closes that gate.

## Archived workflow evidence

The [2.6 MB worked-example archive](../data/validation/worked_example_800_1600.tar.gz) includes the corrected self-contained HTML report, acoustic STEP files, raw and coupled responses, search audit, full rankings, resolved specification, source/input/container manifest and launch/resume logs. Its [identity manifest](../data/validation/worked_example_800_1600_manifest.json) records the exact source revision. The fresh run completed all 18 processes; its subsequent resume reused every process from cache. Earlier invocation failures are explicitly excluded from this pass.

The separate [resolution archive](CANDIDATE_RESOLUTION.md#executed-8001600-hz-example) passes all six fixed comparisons. Neither archive contains a measured physical horn/driver assembly or a manufacturing specification.
