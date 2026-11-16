# Status and roadmap

Updated 2026-11-16. This page replaces the earlier implementation status (7 September 2026) and single-horn roadmap (9 September 2026), which are kept in [history](history/).

## Summary

The software works end to end. Give it a band and it returns ranked driver and horn pairs, dimensions, acoustic STEP files, response, impedance and excursion predictions, rejection reasons and a report, with full run provenance. The predictions are **experimental**: the numerical solver has been checked against independent codes and one measured horn, but no complete horn and driver assembly has been compared with measurement. Every result is labelled accordingly, and a run may legitimately return no feasible design. [Issue #81](https://github.com/timini/horn-simulation/issues/81) tracks the remaining physical qualification.

The scope is conventional single-driver, rigid, axisymmetric, front-loaded horns. Multiple-entry horns belong in [MEH Design Studio](https://github.com/timini/meh-design-studio). Items from the original brief that were not built are listed in [DESCOPED](DESCOPED.md).

## What works

- Single mode: one horn from parameters or an imported STEP (units honoured), optional coupling of one driver. See [USAGE](USAGE.md).
- Auto mode: driver screen, throat sampling, geometry grid over seven profiles, Webster screen, FEM shortlist with profile diversity, bounded refinement with bound widening, ranking with hard limits and a documented score. See [SEARCH_AND_SCORING](SEARCH_AND_SCORING.md).
- Size caps, minimum output, fixed dimensions, profile restriction and custom score weights.
- Driver coupling with front chamber (supplied or estimated for cones), optional chamber duct and sealed rear chamber. Linear excursion and power checks at the drive voltage, and an optional stress power over the whole sweep.
- Per-candidate electrical impedance, excursion, input power and group delay, builder profile CSVs, true-size STEP files, rejection summary with margins, and experimental directivity for `modal_baffled` horns. See [OUTPUTS](OUTPUTS.md).
- `modal_baffled` mouth radiation as the default for generated horns, with `flanged_piston` and other local models as options. Optional thermoviscous wall losses.
- Isolated runs with manifest, source snapshot, input and image hashes; checked resume; run lookup; stale-image warnings.
- Driver databases with provenance, a validator that CI enforces, `horn-add-driver` for datasheet entry, and a resumable scraper. See [DRIVERS](DRIVERS.md).
- Locked image dependencies, pinned Nextflow 24.10.x and nf-test 0.9.5, launcher runs in CI.
- Cross-software comparisons with openwind and bempp-cl, a Hornresp kit, and a frozen comparison with the Hajicek and Neely horn standard. See [VALIDATION](VALIDATION.md).

## What is experimental or not validated

- Absolute level at a listener for a real driver on a real horn.
- The order in which drivers and horns are ranked.
- Driver interfaces: no shipped driver has a verified adapter or phase plug model, or a separated diaphragm mass and rear load. The cone front chamber is a geometric estimate.
- Directivity: infinite-baffle far field only, not compared with measurement; a finite-baffle BEM disagrees.
- Profiles and bands beyond those covered by the resolution studies. One worked candidate does not qualify all seven profiles.
- Wall losses outside the measured pipe regime.
- STEP files are the air volume only, not a buildable part.

## Latest checks

On 2026-11-16, at the merge of the cross-software, Hajicek and CI work:

| Check | Result |
| --- | --- |
| `just test-local` | 745 passed, 1 skipped |
| `pytest scripts/tests` | 139 passed before the documentation tests were added |
| `just test` (Docker) | 458 passed, 5 skipped |
| `just test-nextflow` | 22 tests passed |
| `just validate-drivers` | 0 errors on both databases |
| Image source labels | match the checkout for all three images |

The full openwind and bempp-cl comparison was not rerun at that merge (about 90 and 40 minutes); its committed results are from earlier the same day.

## Remaining work

In order:

1. **Choose a physical reference and a comparison.** Two horn and driver assemblies with exact air and interface geometry, driver variant, rear load, drive voltage, microphone position, environment and reuse terms. [PHYSICAL_REFERENCE_AUDIT](PHYSICAL_REFERENCE_AUDIT.md) lists the candidates; none is complete yet. If the data cannot be found, record what is missing and the measurement that would supply it.
2. **Qualify the driver interface.** Model the chosen adapter and rear load at the fidelity the reference needs, check it against an independent circuit or solver, and qualify a small driver set before growing the catalogue.
3. **Settle radiated output accuracy for the release domain.** `modal_baffled` passes the infinite-baffle reference; finite baffles and free-standing horns do not. Run the resolution harness for each profile and band admitted to the release.
4. **Compare complete assemblies.** Freeze predictions first, then compare absolute response with no fitted level shift, electrical impedance, cutoff and ordering. Differences smaller than the combined uncertainty become ties.
5. **Release a repeatable workflow** with a qualified example, its drivers, dimensions, interface drawing, raw evidence and reproduction commands.

Until step 4 passes, results stay labelled experimental. Adding more profiles, a larger catalogue or more complex geometry is lower priority than this physical evidence.

### Open issues

| Issue | Purpose | Role |
| --- | --- | --- |
| [#81](https://github.com/timini/horn-simulation/issues/81) | Assembly evidence, driver interface, radiation, convergence and final acceptance | Main completion gate |
| [#78](https://github.com/timini/horn-simulation/issues/78) | Independent solver evaluation; interior, cavity and infinite-baffle comparisons done; finite-baffle pilot diagnostic | Numerical part delivered |
| [#74](https://github.com/timini/horn-simulation/issues/74) | LaVoce refresh and catalogue provenance audit | Catalogue maintenance |
| [#45](https://github.com/timini/horn-simulation/issues/45) | Interior field export | Optional, for diagnosis |
| [#49](https://github.com/timini/horn-simulation/issues/49) | Non-axial boundary identification | Deferred |
| [#51](https://github.com/timini/horn-simulation/issues/51) | Folded, back-loaded or phase-plug geometry | Deferred; MEH excluded |

Other known gaps:

- 376 Faital, B&C and 18Sound records with copied pages were removed from `data/drivers` and not yet re-fetched; the source site returned HTTP 429 on 2026-11-15.
- The new CI jobs (`launcher_end_to_end`, `check_lockfile`) have not yet run on GitHub.
- The 3D FEM with `flanged_piston` differs from the 1D codes by up to about 1.5 dB between mouth ka 0.5 and 1.5, which is not explained ([CROSS_SOFTWARE_VALIDATION](CROSS_SOFTWARE_VALIDATION.md)).
- The archived examples were made before the front chamber estimate; under the current coupling their ranking and ripple differ. See [WORKED_EXAMPLE](WORKED_EXAMPLE.md) and [examples](../examples/README.md).

## Decisions

- Driver data are per-driver JSON files in git, not a database server. A single-user local pipeline gains nothing from a server, and git history gives the audit trail.
- No PDF datasheet parser. `horn-add-driver` takes datasheet values with explicit units, checks them and records the source; the person reading the datasheet is the review.
- Local Docker and Nextflow, not a cloud deployment.
- Geometry from gmsh/OpenCASCADE, not FreeCAD; see [history/GEOMETRY_KERNEL_EVALUATION](history/GEOMETRY_KERNEL_EVALUATION.md).
- The FEM-BEM coupling was removed because it coupled the whole boundary rather than the mouth; `modal_baffled` replaced it.

See [DESCOPED](DESCOPED.md) for the reasons in full.
