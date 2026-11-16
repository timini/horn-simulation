# Usage

This page lists every `just` command, the launcher options and every pipeline parameter in `main.nf`. Install the prerequisites in the [README](../README.md#quick-start) first.

## Commands

| Command | What it does |
| --- | --- |
| `just` | List the commands. |
| `just build [tag]` | Build the three Docker images (`horn-solver`, `horn-geometry`, `horn-analysis`) with tag `latest` or the tag given. Each image is labelled with a hash of its sources. |
| `just run [params]` | Run the pipeline in single mode through the launcher. Any parameter below can be added. |
| `just run-auto [params]` | Run auto mode (band search and driver ranking) through the launcher. |
| `just run-fullauto [params]` | Same as `just run-auto`. Kept as an alias for older commands. |
| `just latest-run [--status S] [--root DIR]` | Print the folder of the most recently completed run. |
| `just runs [--root DIR]` | Print a JSON inventory of every run folder. |
| `just test` | Build the test images and run every Docker package test suite. |
| `just test-package <pkg>` | Build and test one Docker package, for example `just test-package horn-solver`. |
| `just test-local` | Run the package tests that need no Docker (horn-drivers, horn-core, horn-analysis, horn-geometry and the gmsh-only solver meshing tests) from `.venv`. |
| `just test-nextflow [args]` | Run the nf-test workflow tests through `scripts/nf_test.sh`. |
| `just validate-drivers` | Check `data/drivers` and `data/drivers-curated` with `horn-validate-drivers`. CI requires zero errors. |
| `just validate-cross-software` | Run the openwind, bempp-cl and Hornresp comparisons. Slow; see [CROSS_SOFTWARE_VALIDATION](CROSS_SOFTWARE_VALIDATION.md). |
| `just validate-cross-software-fem <runs>` | Compare two archived FEM runs with openwind. |
| `just clean` | Remove the Docker images and the Nextflow files in the repository root. Run folders under `results/` are not touched. |

`just test-local` and `just validate-drivers` use `.venv`. Create it with `uv sync --all-packages --extra test`. A plain `uv sync` installs no workspace packages.

## Modes

**Single mode** (the default) simulates one horn. Give the geometry with `throat_radius`, `mouth_radius`, `length` and `profile`, or import a STEP file with `step_file` and `length`. Add `driver_id` to couple one driver to the horn.

```sh
just run --profile exponential --throat_radius 0.025 --mouth_radius 0.15 --length 0.3
```

**Auto mode** takes a target band and returns ranked driver and horn pairs. Any dimension you fix is kept; the others are derived from the band and searched. See [SEARCH_AND_SCORING](SEARCH_AND_SCORING.md) for how candidates are chosen and ranked.

```sh
just run-auto --target_f_low 800 --target_f_high 1600 --drivers_db data/drivers-curated
```

`--mode fullauto` is accepted and runs auto mode. It is kept so older commands still work.

## Launcher

`just run` and `just run-auto` call `scripts/run_pipeline.py`. The launcher:

- creates a new folder `results/<UTC time>-<id>/` for each run, with the outputs under `outputs/`;
- records a `manifest.json` with the source revision, a source archive and patch, input file hashes, Docker image ids and their source labels, the Nextflow version and the exact command;
- picks a Java 17 to 22 runtime for Nextflow without changing your shell settings. It tries `NXF_JAVA_HOME`, `JAVA_HOME`, `JAVA_CMD`, the `java` on `PATH`, then on macOS Homebrew `openjdk@21` and `openjdk@17`;
- warns when a Docker image was built from other sources than the checkout. Rebuild with `just build`.

Launcher options (put them before the pipeline parameters):

| Option | Effect |
| --- | --- |
| `--run-dir PATH` | Use this run folder instead of a generated one. A new run refuses an existing folder. |
| `--image-tag TAG` | Run `horn-*:TAG` images built with `just build TAG`. The environment variable `HORN_IMAGE_TAG` does the same. |
| `--require-current-images` | Fail instead of warning when an image does not match the checkout. CI uses this. |
| `-resume` | Resume the run in `--run-dir`. |
| `-params-file FILE` | Read parameters from a JSON or YAML file. Its hash is recorded. A value given on the command line wins. |

The launcher owns `--outdir`, `-work-dir`, `-log` and `-c`; passing them is an error.

Use a named folder when you want one:

```sh
python3 scripts/run_pipeline.py --run-dir results/my-500hz-horn --mode auto \
  --target_f_low 500 --target_f_high 4000 --drivers_db data/drivers-curated
```

### Resume

```sh
python3 scripts/run_pipeline.py --run-dir results/my-500hz-horn -resume
```

Resume checks that the source, driver and input files, Docker images and Nextflow engine are unchanged, and reuses the original Nextflow session. A change needs a new run. Keep the run's `work/` folder while you may want to resume.

### Find runs

```sh
just latest-run                     # last completed run
just latest-run --status running    # or failed, or any
just runs                           # JSON inventory
just latest-run --root /other/parent
```

`latest-run` orders completed runs by their finish time. Running, failed, malformed and unmanaged folders are excluded unless `--status` asks for them. `runs` lists every immediate child folder, including folders without a manifest (`unmanaged`) and unreadable manifests (`invalid_manifest`). Both commands only read. They also run without `just` as `python3 scripts/runs.py latest` and `python3 scripts/runs.py list`.

"Completed" means the workflow exited cleanly. A completed auto run may still report no feasible design.

### Clean up

Nothing under `results/` is deleted automatically. Before removing a run, check with `just runs` that it is not the only copy of a validation result, archive it outside the repository if you need it, then delete that one folder. Do not clear `results/` or `work/` wholesale.

### Direct Nextflow

You can call Nextflow yourself, but you must choose a fresh output folder and you lose the launcher's manifest and checks:

```sh
nextflow run main.nf -profile docker --outdir results/manual-example/outputs
```

## Pipeline parameters

Every parameter is set with `--name value`. The defaults below are the values in `main.nf`. "unset" means the parameter is `null` and the behaviour described applies. A test in `scripts/tests/test_docs.py` fails when a parameter in `main.nf` is missing from these tables or its default differs.

### Mode and output

| Parameter | Default | Meaning |
| --- | --- | --- |
| `mode` | `single` | `single` or `auto`. `fullauto` is an alias of `auto`. |
| `outdir` | unset | Output folder. Required for direct Nextflow runs; the launcher sets it. |

### Horn geometry

| Parameter | Default | Meaning |
| --- | --- | --- |
| `throat_radius` | unset | Throat radius in metres. Single mode uses 0.05 when unset. In auto mode a value fixes the throat; unset means throats are sampled from the drivers and the band. |
| `mouth_radius` | unset | Mouth radius in metres. Single mode uses 0.2 when unset. In auto mode a value fixes the mouth. |
| `length` | unset | Axial length in metres. Single mode uses 0.5 when unset, and needs it with `step_file`. In auto mode a value fixes the length. |
| `profile` | `conical` | Single-mode profile: `conical`, `exponential`, `hyperbolic`, `tractrix`, `os`, `scaled_tractrix` or `cd`. `lecleach` is a deprecated alias of `scaled_tractrix`. |
| `num_sections` | `20` | Minimum number of profile points in the revolved STEP wall. The generator adds points when a profile needs them to stay within 1 %. |
| `step_file` | unset | Single mode: import this STEP instead of generating geometry. Inlet at z = 0, mouth at z = `length`. Units declared in the file are honoured. |
| `horn_3d_png` | unset | Single mode: use this image in the report instead of the parametric render. |

### Auto search

| Parameter | Default | Meaning |
| --- | --- | --- |
| `target_f_low` | `500` | Lower edge of the target band, Hz. |
| `target_f_high` | `4000` | Upper edge of the target band, Hz. |
| `drivers_db` | `data/drivers` | Driver database folder or JSON file. `data/drivers-curated` is the small hand-checked set. |
| `profiles` | unset | Comma list that restricts the profiles searched, for example `os,conical`. Unset searches all seven. |
| `min_length` | unset | Lower bound on the searched length, m. |
| `max_length` | unset | Upper bound on the searched length, m. With only a maximum, the search runs from 0.4 times it up to it. |
| `min_mouth_radius` | unset | Lower bound on the searched mouth radius, m. |
| `max_mouth_radius` | unset | Upper bound on the searched mouth radius, m. Same 0.4 rule as `max_length`. |
| `num_mouth_radii` | `3` | Mouth radius grid points when the mouth is not fixed. |
| `num_lengths` | `3` | Length grid points when the length is not fixed. |
| `throat_fractions` | unset | Comma list of throat radii as fractions of the largest allowed throat, for example `0.1,0.3,0.65,1.0`. Unset uses `0.3,0.65,1.0` plus driver-matched radii. |
| `min_diameter` | unset | Keep only drivers with nominal diameter at least this, inches. |
| `max_diameter` | unset | Keep only drivers with nominal diameter at most this, inches. |
| `lem_top_n` | `10` | Number of screened geometries sent to FEM. |
| `refinement_budget` | `6` | Extra FEM evaluations for refinement after the shortlist. `0` turns refinement off. |
| `top_n` | `10` | Number of ranked results kept. The report details the top five. |

### Drive, observer and limits

| Parameter | Default | Meaning |
| --- | --- | --- |
| `voltage_rms` | `2.83` | Drive voltage, V RMS, used for levels, excursion and power checks. |
| `observation_distance` | `1.0` | On-axis distance from the mouth plane to the observer, m. |
| `max_ripple_db` | `6.0` | Largest allowed peak-to-peak level variation over the target band, dB. |
| `max_compression_ratio` | `10.0` | Largest allowed ratio of diaphragm area to throat area. |
| `min_output_db` | unset | Optional floor on the lowest level anywhere in the target band, dB. |
| `score_weights` | unset | Preference weights, for example `sensitivity=0.625,ripple=0.25,extension=0.125,size=0`. Unset uses those values. Weights are normalised to sum to 1. |
| `test_power_w` | unset | Optional stress power, W into the nominal impedance. Checks excursion over the whole sweep and reports the Xmax margin. |
| `driver_id` | unset | Single mode: couple this driver from `drivers_db` to the horn and add the coupled response to the report. |

### Solver

| Parameter | Default | Meaning |
| --- | --- | --- |
| `min_freq` | `500` | Single-mode sweep start, Hz. Auto mode derives its sweep from the target band. |
| `max_freq` | `8000` | Single-mode sweep end, Hz. |
| `num_intervals` | `100` | Number of frequency points in the sweep. |
| `min_points_per_octave` | `24` | Raise the point count until every octave has at least this many points. `0` disables. |
| `num_bands` | `8` | Number of frequency bands solved in parallel per geometry. |
| `mesh_size` | `0.01` | Largest element size, m. The solver uses the finer of this and the wavelength rule. |
| `max_mesh_nodes` | `100000` | Node limit per band. A P1 mesh is coarsened towards lambda/6 to fit; if it still does not fit, the candidate is rejected. |
| `radiation_model` | unset | Mouth boundary. Unset uses `modal_baffled` for generated geometry with lossless P1 elements and `flanged_piston` otherwise. Choices: `modal_baffled`, `flanged_piston`, `unflanged_piston`, `finite_flange`, `plane_wave`, `closed` (single mode only). |
| `loss_model` | `lossless` | `lossless` or `boundary_layer` (thermoviscous wall losses, opt-in). |
| `minimum_wall_scale` | unset | Required with `boundary_layer`: the smallest wall curvature or gap scale, m. |
| `element_degree` | `1` | Finite element degree, 1 or 2. `modal_baffled` needs 1. |
| `flange_width` | `0.0` | Flange width for `finite_flange`, m. |

See [PHYSICS](PHYSICS.md) for what the solver and radiation models compute.

## Imported STEP geometry

```sh
just run --step_file path/to/horn.step --length 0.3
```

The inlet must be at z = 0 and the mouth at z = `length`. The file's declared length unit is honoured, so millimetre CAD imports at true size. The report uses the actual CAD inlet and mouth areas. An imported STEP uses `flanged_piston` unless you choose another model, and the report prints that model's known throat impedance error.

## Other command-line tools

The packages install console scripts into `.venv/bin`. The pipeline calls most of them as Python modules inside the containers; they are also useful on their own.

| Script | Package | Purpose |
| --- | --- | --- |
| `horn-geometry` | horn-geometry | Write a horn STEP from profile and dimensions. Needs gmsh. |
| `horn-solver` | horn-solver | Solve one STEP over a frequency range. Needs dolfinx, so run it in the `horn-solver` image. |
| `horn-candidates` | horn-core | Write a candidate geometry grid as CSV. |
| `horn-add-driver` | horn-drivers | Add a driver from datasheet values. See [DRIVERS](DRIVERS.md). |
| `horn-validate-drivers` | horn-drivers | Check a driver database. |
| `horn-scrape-drivers` | horn-drivers | Resumable import from the loudspeakerdatabase.com pages. Needs the `scrape` extra. |
| `horn-prescreen` | horn-analysis | Filter drivers and choose throat radii for a band. |
| `horn-lem-prescreen` | horn-analysis | Score a candidate grid with the Webster model. |
| `horn-screen-drivers` | horn-analysis | Screen drivers against one solved horn. |
| `horn-rank` | horn-analysis | Rank drivers against one solved horn. |
| `horn-auto-report` | horn-analysis | Build the auto report from ranked results. |
| `horn-single-report` | horn-analysis | Build the single-mode HTML report. |
| `horn-plotter` | horn-analysis | Plot SPL against frequency from a results CSV. |
| `horn-impedance-plot` | horn-analysis | Plot throat impedance. |
| `horn-phase-plot` | horn-analysis | Plot phase. |
| `horn-dashboard` | horn-analysis | Combined SPL, impedance and phase plot. |
| `horn-kpi` | horn-analysis | Extract peak-lobe KPIs to JSON. |
| `horn-compare` | horn-analysis | Overlay several results CSVs. |
| `horn-render` | horn-analysis | Render the horn in 3D as PNG. |
| `horn-directivity-plot` | horn-analysis | Directivity plots from a directivity CSV. |

Each takes `--help`.

## Scripts

| Script | Purpose |
| --- | --- |
| `scripts/run_pipeline.py` | The launcher described above. |
| `scripts/runs.py` | `list` and `latest` run lookup. |
| `scripts/check_run.py` | Check a run folder's status, expected outputs and provenance. Used by CI. |
| `scripts/image_source.py` | Hash the sources of a Docker image. Used by `just build` and the launcher. |
| `scripts/nf_test.sh` | Download nf-test 0.9.5 to `.tools/`, check its sha256 and run it. `NF_TEST` overrides the binary. |
| `scripts/validate_drivers.py`, `scripts/scrape_drivers.py` | Thin wrappers around `horn-validate-drivers` and `horn-scrape-drivers`. |

The validation scripts are listed in [VALIDATION](VALIDATION.md).
