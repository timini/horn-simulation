# 6NMB420: 320–5000 Hz, maximum 300 mm horn length

**Selected experimental design:** OS flare, 300 mm acoustic length, 76 mm throat diameter and 244.80268 mm mouth diameter, using the **8 Ω** 18Sound 6NMB420. The 239 mm mouth candidate is effectively tied; this is a finite search, not proof of a global optimum.

Open [report.html](report.html), download [horn.step](horn.step) for CAD, or inspect [design-overview.png](design-overview.png). The 300 mm length excludes the driver, mounting adapter and rear enclosure. The STEP is the acoustic air volume; walls, mounting, adapter and phase plug are not supplied.

The final prediction is **4.42 dB peak-to-peak variation** over 320–5000 Hz, averaging **105.52 dB SPL at 2.83 V RMS / 1 m on axis**, without EQ. [Response and electrical impedance](response.csv), [nominal profile coordinates in mm](profile.csv), [full assessment](design.json), [numerical checks](resolution.json) and [CAD dimensions](cad-check.json) are included.

Mesh and frequency checks passed: maximum in-band SPL change 0.244 dB, impedance magnitude 0.397 dB and phase 1.609°, with 0.088 dB interpolation difference after doubling the frequency samples. This compares 8 and 6 mm meshes using the same 40-section CAD; it does not establish every form of discretization error.

## Limits

The ideal model uses a lossless horn, infinite rigid baffle, modal mouth radiation, and [manufacturer driver parameters](https://www.eighteensound.it/es/products/lf-driver/6-5/8/6NMB420). The real adapter, rear arrangement, separated moving/rear air masses and cone breakup are uncharacterized. No off-axis coverage, distortion, maximum-SPL or physical assembly qualification is claimed. Source measurements would be needed before treating this as a dependable production assembly.

## Reproduce

From the repository root, build the images with `just build`. The search that found the selected family is:

```bash
just run-auto --target_f_low 320 --target_f_high 5000 \
  --drivers_db examples/6nmb420-320-5000/driver-database.json \
  --throat_radius 0.038 --min_length 0.1 --max_length 0.3 \
  --num_lengths 5 --num_mouth_radii 5 --lem_top_n 2 --refinement_budget 6 \
  --num_bands 2 --num_intervals 161 --mesh_size 0.008 \
  --num_sections 40 --radiation_model modal_baffled
```

A broad 700-geometry analytical search preceded the focused 175-geometry search. Its eight FEM finalists failed the ripple limit; the narrow-throat search found feasible designs. Six bounded refinement trials produced five valid FEM sweeps and one recorded CAD rejection. In total, 15 distinct geometries received complete FEM sweeps. The automatic broad shortlist did not find this family by itself.

Reproduce the selected geometry at the finer resolution:

```bash
just run --step_file examples/6nmb420-320-5000/solver-input.step \
  --profile os --length 0.3 --throat_radius 0.038 --mouth_radius 0.12240134 \
  --drivers_db examples/6nmb420-320-5000/driver-database.json \
  --driver_id 18sound-6nmb420 --target_f_low 320 --target_f_high 5000 \
  --min_freq 226.2741699796952 --max_freq 7071.067811865476 \
  --num_bands 2 --num_intervals 321 --mesh_size 0.006 \
  --num_sections 40 --radiation_model modal_baffled
```

The sweep extends beyond the target band to show roll-off. Assess the requested 320–5000 Hz band, rather than the padded sweep or a peak-relative contiguous −3 dB region. New runs use new output directories and do not inherit physical qualification.

## CAD units and provenance

**Use `horn.step` for CAD.** The pipeline's preserved `solver-input.step` has numeric metre coordinates but a legacy millimetre declaration. The physical export changes `SI_UNIT(.MILLI.,.METRE.)` to `SI_UNIT($,.METRE.)`, leaving the geometry coordinates unchanged. Re-importing that physical export in millimetres verifies the dimensions in `cad-check.json`. Use the preserved solver input for the current pipeline, not the physical export. `profile.csv` describes the nominal mathematical curve; the STEP uses its 40-section loft approximation.

`search-manifest.json`, `fine-manifest.json` and `wide-search-manifest.json` identify the actual source, inputs, container images, commands and completion-time output hashes. The executed refined/fine package source is commit `fea77f8`; earlier broad search source is `9febae1`. `example-manifest.json` maps the selected raw results to their sealed originating outputs and hashes the delivered files. Source archives and interrupted exploratory attempts remain in the original local run directories; this directory is a compact result package, not resumable Nextflow state.
