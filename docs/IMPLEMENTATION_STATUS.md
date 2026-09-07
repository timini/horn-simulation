# Implementation and validation status — issue 81

Updated 7 September 2026. Band-to-driver-and-horn search, bounded refinement, CAD export and reports run automatically. **Physical recommendation validation is not complete**, and issue 81 remains open. Reports distinguish model feasibility from missing evidence and do not certify purchases or builds.

## Implemented

- Isolated runs with source archive/patch, input hashes, actual container IDs, command, run status and exact Nextflow session ID. Resume rejects changed source, data, images or arguments.
- Consistent RMS pressure, velocity, impedance and voltage conventions across analytical screening, FEM and motor coupling. Actual CAD inlet/mouth areas include annular openings. The listening-distance approximation is explicitly a uniform baffled piston.
- Corrected FEM radiation sign, boundary reaction for inlet impedance, quadratic geometry/field elements, fail-closed linear-solve residuals, and exported power accounting. Optional thermoviscous wall losses and finite-flange pipe radiation have matched-reference tests.
- Complete-band merge validation: missing bands, missing samples, inconsistent metadata, nonfinite responses and excessive overlap discrepancies terminate the run.
- Requested-band ripple/output constraints, known driver-band and compression rejection, missing-provenance/interface/mass flags, and linear excursion/input-power rejection at the requested voltage.
- Fixed-budget analytical shortlist with score and geometry diversity, bounded refinement, explicit boundary/budget status, and near-tied score reporting. A result is the best evaluated candidate, never a proven global optimum.
- Explicit empty and rejected outcomes; coupled response files and acoustic air-volume CAD; default assumptions and experimental labels in reports.
- Three manufacturer-sourced motor variants in `data/drivers-curated`, with supplied/derived/missing fields distinguished. Their throat interfaces and separate air loads are unresolved; they are not promoted to physically validated drivers.
- Published reference importer with archive hashes, units, phase availability and duplicate detection. Six acquired archives contain 383 curves, 378 unique files: 299 research measurements, 40 independent simulations and 44 DIY response/impedance files.
- Disabled the invalid legacy FEM–BEM horn coupling at public entry points. It mapped the whole boundary instead of a verified mouth-only exterior problem. Independent backend/operator tests remain available.
- Native container builds, repaired Nextflow tests, pinned Nextflow/nf-test CI versions, and a dedicated production acoustic lane requiring all 13 declared cases without skips.

## Demonstrated evidence

| Check | Result | Evidence |
| --- | --- | --- |
| Clean-checkout local packages before review | 309 passed | `/tmp/horn-clean-local.xml` |
| Local packages after review corrections | 398 passed | `results/validation-81/review9-local.xml` |
| Complete Nextflow suite after review fixes | All eighteen tests passed in one run, including loss/refinement, finite-flange/unflanged domain rejection, custom render paths, imported STEP coupling, default single-mode provenance and missing-band rejection | `results/validation-81/review9-nextflow.xml` |
| Expanded finite-grid search | Four bands, 20 geometries × three manufacturer motors; 32/40/38/28 feasible pairs after the per-mouth fit correction; 100% top-ten recall and zero regret in all cases | `data/validation/search_benchmark_summary.json` |
| Full solver suite after wall-loss/residual changes | 48 passed, one optional straight-tube regime skipped | `results/validation-81/loss-physics/all-solver-tests.xml` |
| Three-mesh production wall-loss case | Complex impedance, finest-mesh output convergence, passivity and RMS energy balance pass | `packages/horn-solver/tests/validation/test_wall_losses.py` |
| Independent published numerical impedances | All four comparisons pass; p95 magnitude error <0.010 dB | `data/validation/loss_physics_validation.json` |
| Published measured impedances | 237/299 pass magnitude limits, including 201/263 held-out curves; failures retained | [Full comparison and limits](LOSS_PHYSICS_VALIDATION.md) |
| Direct production FEM versus measurements | 90/100 comparisons pass on the three open geometries | `results/validation-81/loss-physics/held-out-final/validation.json` |
| Fresh clean-checkout band-only run | 235 screened drivers, 231 geometries; 106 tasks completed in 6m28s. All 80 original simulation bands also pass the stricter complex-join gate | `data/validation/band_only_release_summary.json` |
| Earlier band-only database run | 239 drivers, 231 geometries / 55,209 analytical pairs; 106 tasks completed in 6m23s | `results/validation-81/band-only-final/manifest.json` |
| Clean-checkout interruption/resume | Interrupted after three tasks; resumed the same session; numerical/ranking files identical to an uninterrupted run. Other differences limited to timestamps/output directory | `data/validation/resume_validation_summary.json` |
| Earlier exact resume | All 14 tasks cached; output hashes unchanged | `results/validation-81/reproducible-auto/manifest.json` |
| Earlier finite-grid audit | 20 geometries × three synthetic motors: all six feasible pairs retained, winner regret zero | `results/validation-81/search-benchmark/search_benchmark.json` |

The earlier band-only run produced a Beyma 3FR30Nd/hyperbolic experimental candidate with 19.868 mm throat diameter, 141.934 mm mouth diameter and 171.5 mm length. It had insufficient driver/interface evidence. This is historical execution evidence, not a recommendation or evidence for the latest source snapshot. Likewise the earlier resume and six-feasible-pair benchmark do not replace final release checks.

Review corrections also enforce each independently known driver-band endpoint, constrain custom throat fractions, propagate single-driver voltage, reject inconsistent complex band joins, and report valid-but-empty size searches. The custom run directory is excluded from source hashing so it cannot invalidate its own resume.

## Outstanding gates

1. Reproduce a fully characterized horn-and-driver assembly and a held-out comparison assembly, including chamber/adapter, source voltage and observation environment. Pipe impedance agreement cannot establish absolute loudspeaker output or driver ranking.
2. Establish verified driver interfaces and separate diaphragm/rear-load parameters. Manufacturer data improve traceability but do not fill those missing quantities.
3. Investigate failed material/closed-end impedance comparisons. The new loss model is substantially better but does not pass every measured case; do not broaden its physical domain by adjusting thresholds.
4. Validate the listening-distance/exterior approximation beyond ideal piston references. Directivity and arbitrary free-standing horn radiation are outside the supported release.
5. Retain complete final search and clean-checkout integration evidence with the release. Acoustic CAD still needs mechanical design before fabrication.

## Reproduce

```sh
just build
just test-local
just test-nextflow
just run-auto --target_f_low 1000 --target_f_high 1200
```

Use `--drivers_db data/drivers-curated` for traceable manufacturer inputs. Losses are opt-in with `--loss_model boundary_layer --minimum_wall_scale 0.01 --element_degree 2` only where that scale is conservative for the actual geometry. See [the acoustic contract](ACOUSTIC_CONTRACT.md) and [loss validation](LOSS_PHYSICS_VALIDATION.md).

The launcher selects an available compatible Java without changing global settings. The checked-in validation scripts emit protocols, complete comparisons and hashes. Large third-party archives are downloaded locally from pinned sources rather than vendored without reuse permission.

The second review adds a fixed-radius acoustic cap, per-mouth driver fit in screening and ranking, unsmoothed ripple extrema, explicit imported STEP length, profile-preserving refinement IDs and shared moving-mass selection for Neumann/pressure coupling. All four search bands were reassessed against the frozen FEM responses with the updated eligibility rules; the final summary retains response/code hashes and separates original simulation runtime from reassessment.

The third review corrects imported-geometry and variable-throat reporting, separates shape-only comparison from the absolute-level gate, and records/pins the Nextflow engine and launcher identity for resume. Reports regenerated from saved numerical results retain their original simulation source identity.

The fourth review removes the last fallback-radius label from imported driver-coupling plots/metadata and records effective default dimensions for generated single-mode runs. Imported runs preserve unknown radii in the resolved specification and use actual inlet area for coupling.

The fifth review makes radiation-domain eligibility a per-geometry check in screening and refinement. Mixed finite-flange grids keep supported designs; all-invalid grids produce an explicit report with retained rejection records. New workflow cases verify both outcomes and distinguish rejected trials from FEM evaluations.

Refinement totals count actual FEM evaluations and coupled driver pairs. Analytically rejected trials remain in the audit but do not inflate those totals; all-invalid reports expose analytical attempt counts separately from zero FEM totals.

The final review corrections preserve the exact Le Cléac’h throat in the shared screening/rendering/CAD profile, stage custom render images under a safe fixed name, reject unsupported unflanged candidates individually, and show unavailable dimensions in early-exit reports. General curve comparison uses a fixed 10,001-point log-frequency grid for percentile gates, with original knots retained for exact maximum-error reporting; uneven source sampling cannot change the gate weights. The published pipe benchmark already uses its separately frozen comparison protocol and is unchanged. Earlier full seven-profile search runs retain their historical source identity; their Le Cléac’h geometries predate the inlet correction. The finite-grid search benchmark uses conical and exponential geometries and is unaffected by this profile correction.

The subsequent packaging review updates the frozen workspace lock for the shared geometry dependency and stages imported STEP inputs under a safe fixed name in each solver process. A clean `uv sync --frozen --package horn-geometry` environment successfully exports the 1 mm Le Cléac’h inlet, and the imported assembly workflow exercises a STEP basename containing spaces and a quote.

The next review validates every supplied numeric driver field (positive physical limits and ordered usable-band bounds; zero mechanical loss/rear air mass remain explicit limiting cases). Invalid or impossible derived parameters reject that database record with a warning; requesting an invalid driver directly raises its diagnostic. Refinement now reapplies the requested `top_n` result limit while preserving complete evaluation totals and the best score. The complete refinement fixture requests one result and verifies only one is returned after multiple evaluations.
