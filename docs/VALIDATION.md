# Validation

This page indexes every validation record, script and CI job, with its current result. Status is as of 2026-11-16.

The short version: the numerical solver agrees with independent solvers and with one measured horn's input impedance in a bounded range. No complete horn and driver assembly has been compared with measurement. Predicted levels at a listener, and the order in which drivers are ranked, are therefore not validated. Failures are kept and listed; no tolerance has been loosened to obtain a pass.

## What is and is not validated

| Claim | Status | Evidence |
| --- | --- | --- |
| Interior FEM matches an independent FEM and motor model on shared meshes | Pass | [Independent solver](#numerical-solver) |
| Production volume operator gives correct cavity modes | Pass | [Independent solver](#numerical-solver) |
| `modal_baffled` mouth matches an independent mode-matching code (one horn, 800 to 1600 Hz) | Pass | [Radiation](#radiation) |
| `flanged_piston` mouth matches the same reference | Fail on throat impedance (1.05 dB, 5.7 degrees); on-axis pressure within 0.40 dB | [Radiation](#radiation) |
| Measured input impedance of one exponential horn, ka 1 to 5 | Pass | [Measured horns](#measured-horns) |
| Measured reflectance of one conical horn standard (Hajicek and Neely) | Fail on frozen gates | [Measured horns](#measured-horns) |
| Measured pipe impedance with wall losses | Partial: 237 of 299 curves | [Measured pipes and losses](#measured-pipes-and-losses) |
| Webster screen against openwind, same physics | Pass, 35 of 35 | [Cross-software](#cross-software) |
| Screen keeps the best FEM designs on a fixed grid | Pass, four bands, both radiation models | [Search and workflow](#search-and-workflow) |
| Worked candidates are converged in mesh, loft and frequency grid | Pass for the two worked candidates | [Search and workflow](#search-and-workflow) |
| Absolute SPL of a driver on a horn | Not validated | No matched measured assembly yet; see [STATUS](STATUS.md) |
| Driver ranking order | Not validated | Same |
| Directivity | Not validated | Finite-baffle BEM disagrees, dominated by baffle diffraction |

## Numerical solver

| Check | Result | Record | Script | Write-up |
| --- | --- | --- | --- | --- |
| Interior FEM and ideal motor against Boundary Lab (Julia), tube and cone, three shared meshes, 13 frequencies | Pass, max relative complex difference below 9e-13 | `data/validation/boundary_lab_reference.json`, `_manifest.json`, `_protocol.json`, `_reference_artifacts.tar.gz` | `scripts/validate_boundary_lab.py`, `scripts/boundary_lab_fixture.py` | [INDEPENDENT_SOLVER_VALIDATION](INDEPENDENT_SOLVER_VALIDATION.md) |
| Rectangular cavity, 12 modes, three P2 meshes | Pass, finest max error 0.033 % | test `test_cavity_modes.py` | `just test-package horn-solver` | [INDEPENDENT_SOLVER_VALIDATION](INDEPENDENT_SOLVER_VALIDATION.md) |
| Uniform baffled piston against ten Rayleigh integrals | Pass | `packages/horn-core/tests/test_rayleigh_integral.py` | `just test-local` | [ACOUSTIC_CONTRACT](ACOUSTIC_CONTRACT.md) |
| Production acoustic cases (contract, velocity inlet, wall losses, Webster horns, modal production, mesh geometry) | 28 cases required without skips | horn-solver tests | CI job `production_reference_cases` | [ACOUSTIC_CONTRACT](ACOUSTIC_CONTRACT.md) |

## Radiation

| Check | Result | Record | Script | Write-up |
| --- | --- | --- | --- | --- |
| `flanged_piston` horn against MMM toolbox (pinned, Octave), infinite baffle, 33 frequencies | Throat impedance fails (1.05 dB, 5.73 degrees); on-axis pressure within 0.40 dB | `data/validation/mmm_radiation_*` | `scripts/validate_mmm_radiation.py`, `scripts/mmm_radiation_reference.m`, `scripts/mmm_reference.Dockerfile` | [EXTERIOR_RADIATION_VALIDATION](EXTERIOR_RADIATION_VALIDATION.md) |
| `modal_baffled` against the same reference | Pass: 0.017 dB impedance, 0.013 dB pressure, below 0.19 degrees | `data/validation/modal_mouth_*` | `scripts/validate_modal_mouth.py` | [MODAL_APERTURE](MODAL_APERTURE.md), [history](history/MODAL_APERTURE_VALIDATION.md) |
| Finite-baffle pilots | Diagnostic only: several dB from the infinite baffle | `data/validation/exterior_baffle_pilot_*` | `scripts/audit_exterior_pilot.py` | [EXTERIOR_RADIATION_VALIDATION](EXTERIOR_RADIATION_VALIDATION.md) |
| `modal_baffled` domain rejection (2 m mouth) | Completes with an empty report and no FEM tasks | `data/validation/modal_empty_workflow.tar.gz` | | [MODAL_APERTURE](MODAL_APERTURE.md) |

## Measured horns

| Check | Result | Record | Script | Write-up |
| --- | --- | --- | --- | --- |
| Post 1994 exponential horn input impedance, ka 1 to 5, three meshes | Pass: finest max 0.497 dB, 1.873 degrees | `data/validation/post_hixson_*` | `scripts/validate_post_hixson.py` | [POST_HIXSON_VALIDATION](POST_HIXSON_VALIDATION.md) |
| Hajicek and Neely conical horn reflectance (Zenodo 19710026), frozen and unfitted | Fail: analytic cone check passes 4 of 4; windowed reflectance 0 of 4, impedance 0 of 2, FEM against TMM fails | `data/validation/hajicek_horn_reflectance.json` | `scripts/validate_horn_reflectance.py`, `scripts/summarise_horn_reflectance.py` | [VALIDATION_DATA](VALIDATION_DATA.md#conical-horn-reflectance-standard-hajicek-and-neely) |
| Sa and Park 2014 sixth resonance of a conical horn | 1712.6 Hz predicted against 1712 Hz reported; narrow check | `data/validation/reference_audit_2026_09_08.json` | `scripts/validate_published_horn_resonance.py` | [history/REFERENCE_VALIDATION_REPORT](history/REFERENCE_VALIDATION_REPORT.md) |
| Physical assembly references (Pi Speakers, Makarski, Cinanni and others) | None yet supplies a complete driver, interface and calibrated measurement | `data/validation/reference_audit_evidence.json` | | [PHYSICAL_REFERENCE_AUDIT](PHYSICAL_REFERENCE_AUDIT.md) |

## Measured pipes and losses

| Check | Result | Record | Script | Write-up |
| --- | --- | --- | --- | --- |
| Lossless brass pipe baseline | Fail, kept as baseline: first peaks 25 to 32 dB too high | `data/validation/pipe_diagnostic_baseline.json` | `scripts/benchmark_measured_impedance.py` | [LOSS_PHYSICS_VALIDATION](LOSS_PHYSICS_VALIDATION.md) |
| Wall losses against published numerical impedances | Pass, 4 of 4 | `data/validation/loss_physics_validation.json` | `scripts/validate_duct_physics.py` | [LOSS_PHYSICS_VALIDATION](LOSS_PHYSICS_VALIDATION.md) |
| Wall losses against measured pipe impedances | Partial: 237 of 299 (201 of 263 held out); failures kept | same | same | [LOSS_PHYSICS_VALIDATION](LOSS_PHYSICS_VALIDATION.md) |
| Full reference audit, 8 September 2026 | TMM 237 of 299; production FEM 240 of 299; independent simulations 40 of 40 | `data/validation/reference_audit_2026_09_08.json`, `reference_search_definition.json` | `scripts/validate_reference_audit.py`, `scripts/render_reference_audit.py` | [history/REFERENCE_VALIDATION_REPORT](history/REFERENCE_VALIDATION_REPORT.md) |
| Measurement catalogue and importer | 383 curves from six archives, plus the Hajicek record | `data/validation/references.json` | `python -m horn_analysis.reference_data` | [VALIDATION_DATA](VALIDATION_DATA.md) |

## Cross-software

All from [CROSS_SOFTWARE_VALIDATION](CROSS_SOFTWARE_VALIDATION.md). `just validate-cross-software` runs `validation/cross_software/openwind_compare.py`, which writes `data/validation/cross_software_openwind.json`, `validation/cross_software/bempp_compare.py`, which writes `data/validation/cross_software_bempp.json`, and `validation/hornresp/hornresp_records.py`. `just validate-cross-software-fem` runs `validation/cross_software/fem_compare.py` on archived runs and writes `data/validation/cross_software_fem.json`.

| Comparison | Result |
| --- | --- |
| Webster TMM against openwind, same physics, 7 profiles, 5 cases | Pass, 35 of 35 |
| Webster TMM piston load against openwind flanged-pipe fit | Fails 6 of 7 profiles; expected model difference |
| 3D FEM against openwind and TMM, conical and exponential | Fails the frozen low-ka limit (0.60 and 0.40 dB median). The difference between ka 0.5 and 1.5, up to about 1.5 dB, is not explained |
| bempp-cl pulsating sphere | Pass |
| Modal far field against direct Rayleigh quadrature | Pass, below 1e-12 dB |
| bempp-cl piston in a finite baffle | Impedance passes at ka 2 and 4; directivity fails at ka 2 (finite baffle) |
| bempp-cl horn in a finite baffle against `modal_baffled` FEM | Throat impedance passes 5 of 5; directivity and on-axis level fail, dominated by baffle diffraction |
| BEM mesh refinement | 7 of 9 pairs pass |
| Hornresp | Not run (Windows only). Input sheets and a comparison script are in `validation/hornresp/`; see [its README](../validation/hornresp/README.md) |

## Search and workflow

| Check | Result | Record | Script | Write-up |
| --- | --- | --- | --- | --- |
| Screen against exhaustive FEM, 20 horns, 3 motors, 4 bands, `flanged_piston` | Pass: 100 % top-ten recall, winner kept, zero regret | `data/validation/search_benchmark_summary.json` | `scripts/benchmark_search.py` | [SEARCH_AND_SCORING](SEARCH_AND_SCORING.md#search-quality-evidence) |
| Same with `modal_baffled` | Pass | `data/validation/modal_search_*` | `scripts/benchmark_search.py` | [history/MODAL_SEARCH_VALIDATION](history/MODAL_SEARCH_VALIDATION.md) |
| Resolution of the 800 to 1600 Hz candidate, `flanged_piston` | Pass, 6 of 6 | `data/validation/candidate_resolution_*` | `scripts/validate_candidate_resolution.py` | [CANDIDATE_RESOLUTION](CANDIDATE_RESOLUTION.md) |
| Resolution of the same band, `modal_baffled` | Pass, 6 of 6 | `data/validation/modal_candidate_resolution_*` | same | [MODAL_CANDIDATE_RESOLUTION](MODAL_CANDIDATE_RESOLUTION.md) |
| Clean-checkout band-only run | Completed: 235 drivers, 231 geometries, 106 tasks | `data/validation/band_only_release_summary.json` | launcher | [history/IMPLEMENTATION_STATUS](history/IMPLEMENTATION_STATUS.md) |
| Interrupt and resume | Resumed results identical to an uninterrupted run | `data/validation/resume_validation_summary.json` | launcher | [history/IMPLEMENTATION_STATUS](history/IMPLEMENTATION_STATUS.md) |
| Archived worked examples | Complete with cached resume | `data/validation/worked_example_800_1600*`, `worked_example_modal_800_1600*` | launcher | [WORKED_EXAMPLE](WORKED_EXAMPLE.md), [examples](../examples/README.md) |
| LaVoce driver classification audit | Corrected (issue 73) | `data/validation/lavoce_classification_audit.json` | | |

Most records were made with earlier source revisions. Each record or manifest names its revision. Code changes since then (for example the front chamber model) do not affect the solver checks, but the archived worked examples rank differently under the current coupling; see [WORKED_EXAMPLE](WORKED_EXAMPLE.md).

## Tests

| Suite | Command | Last local result (2026-11-16) |
| --- | --- | --- |
| Local package tests | `just test-local` | 745 passed, 1 skipped |
| Script tests, including the documentation checks | `.venv/bin/python -m pytest scripts/tests` | 190 passed |
| Docker package tests | `just test` | 458 passed, 5 skipped |
| Nextflow workflow tests | `just test-nextflow` | 22 tests passed |
| Driver databases | `just validate-drivers` | 0 errors |

## CI jobs

`.github/workflows/ci.yml` runs on pushes to `main` and on pull requests:

| Job | What it checks |
| --- | --- |
| `test_nextflow` | Builds the images and runs nf-test 0.9.5 on Nextflow 24.10.5. |
| `launcher_end_to_end` | Runs a tiny single-mode case and a tiny auto case through the launcher with `--require-current-images`, then checks the outputs with `scripts/check_run.py`. |
| `check_lockfile` | `uv lock --check`: the images install from `uv.lock`, so it must match the `pyproject.toml` files. |
| `test_pure_python` | horn-core and horn-drivers tests on Python 3.11. |
| `validate_driver_databases` | Zero validator errors on `data/drivers` and `data/drivers-curated`. |
| `test_reference_audit` | `pytest scripts/tests`, and the run locator tests on Python 3.10. |
| `test_packages` | horn-geometry, horn-solver and horn-analysis tests in their Docker test images (slow `validation` tests skipped on pull requests). On `main` it also pushes the production images. |

`.github/workflows/acoustic-validation.yml` runs `production_reference_cases`: the 28 supported acoustic cases in the solver image, with no skips allowed. `.github/workflows/pages.yml` publishes `site/` to GitHub Pages.

The `launcher_end_to_end` and `check_lockfile` jobs were added on 2026-11-16 and have been run locally but not yet on GitHub.
