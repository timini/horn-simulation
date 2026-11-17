# Horn pipeline completion and validation plan

Prepared 6 September 2026 from inspection of the local working tree, saved runs, and local tests.

## Objective and definition of done

Given a target frequency band, produce ranked, physically credible driver-and-horn recommendations with dimensions, predicted performance, limitations, and reproducible evidence. Preserve the current geometry, solver, database, orchestration, and reporting infrastructure; rebuild selection around validated acoustic outputs.

The band-only interface must remain available. A versioned default specification supplies missing choices such as drive voltage, observation distance, radiation environment, response tolerance, and size bounds. Every report shows these assumptions. Optional user constraints override them. The program must be able to return **no feasible design** or **insufficient evidence**, rather than always declaring a winner.

Two completion gates:

- **Validated software prototype:** analytical, numerical, independent-reference, and search-quality checks pass within an explicitly supported modelling domain. Suitable for comparing candidates, with prediction limits disclosed.
- **Validated design workflow:** at least one recommended physical assembly and one comparison assembly are measured under controlled conditions, with agreement inside the agreed tolerances. Claims remain limited to the tested driver/interface/geometry domain.

A generated acoustic air-volume STEP file is not a manufacturing-ready horn. Report it accordingly. Walls, mounting, adapter/chamber construction, and fabrication drawings require separate mechanical work before a prototype can be built.

## Initial scope

Begin with rigid, axisymmetric, front-loaded horns and a small curated driver set with known parameters and a reproducible throat interface. Validate one simple reference assembly before expanding driver classes. A cone driver requiring a chamber or phase plug cannot be treated as directly connected to an arbitrary throat without modelling or validating that interface.

Validate simple conical and exponential cases first. Retain all seven profile implementations, but only promote a profile into the validated search when it passes the relevant geometry, acoustic, and convergence checks.

Automatic phase-plug design, coaxial tweeter integration, folded horns, manufacturing automation, and nonlinear distortion prediction are later extensions. The existing F1-mid assembly remains an experimental example. Its throat-area and report defects are corrected before interpreting its coupled results.

## Observed baseline

- Band-driven candidate generation, analytical screening, FEM evaluation, coupling, ranking, and reports already exist in `main.nf`.
- The local database loads 1,876 drivers from 21 manufacturers; only nine have a loaded continuous-power field. Data availability is not equivalent to model suitability.
- Local core, driver, and analysis tests: 238 passed, eight directivity tests failed at an import-time NumPy compatibility error.
- Geometry tests require missing local gmsh. Docker was not running, so current container tests and a fresh complete pipeline were not verified.
- The working tree contains substantial uncommitted work. Existing saved reports do not certify the present source and dependencies.
- Fast screening scores throat pressure; final ranking scores mouth pressure. Neither ranking path establishes listening-distance sensitivity.
- Automatic FEM does not forward the radiation-model option and ignores simulation failures; merging does not require all frequency bands.
- Custom single-driver coupling infers a circular throat area even when the imported geometry has an annular inlet.

## Milestone 1 — Reproducible baseline and truthful reporting

**Work**

- Inventory tracked and untracked work; preserve it in a reviewed development snapshot without discarding existing changes.
- Record source revision plus working-tree patch, dependency versions, container image digests, and driver-database digest.
- Fix NumPy compatibility, restore reproducible geometry/solver environments, and align build recipes with configured container names.
- Repair stale Nextflow process tests, including geometry inputs, and add a tiny complete auto-workflow fixture.
- Correct mouth/throat pressure labels immediately; remove unsupported 1 W / 1 m equivalence claims.
- Require valid band bounds, positive dimensions, usable drivers, and nonempty candidate sets. Diagnose malformed data without silently accepting it.
- Fail or explicitly invalidate any candidate missing a simulation band. Validate expected frequencies, duplicates, finite values, and target-band coverage before ranking.
- Make requested radiation settings reach every relevant solver process.

**Exit gate**

All required unit and container integration tests pass. Required validation tests cannot silently skip dependencies in the validation job. A small automatic run completes from a clean checkout and exports a provenance manifest. An injected band failure prevents a valid recommendation for that candidate. No source changes are required to reproduce the run.

**Primary files:** `main.nf`, `nextflow.config`, `justfile`, package Dockerfiles, `tests/`, `.github/workflows/ci.yml`, plotting/report modules.

## Milestone 2 — One consistent acoustic contract

**Work**

- Define units and conventions for complex pressure, velocity, volume velocity, impedance, electrical input, RMS/peak amplitudes, and acoustic power.
- Export actual inlet and mouth areas and boundary identities from geometry/solver results. Driver coupling consumes the actual interface area, including annular openings.
- Extend the analytical transfer-matrix model to return the transfer quantities needed for the same output metric as FEM, rather than scoring throat pressure alone.
- Verify the driver coupling equations and interface assumptions against independently derived limiting cases.
- Select and validate a radiation treatment for the supported domain. Predict absolute pressure at a stated observation point; do not obtain a nominal one-metre value by relabelling mouth pressure.
- Keep internal pressure, radiated output, directivity, voltage sensitivity, and electrical power as distinct outputs. Document any far-field assumptions and their applicable distances.

**Exit gate**

Both models consume the same specification and emit the same physical quantities with explicit conventions. Circular and annular fixtures use verified areas. Tests establish correct voltage scaling, area transformations, impedance signs, and energy accounting under their applicable assumptions. An independently implemented reference agrees with a simple driver/interface case.

**Primary files:** `horn_core/webster.py`, `horn_analysis/transfer_function.py`, `horn_analysis/couple_single.py`, `horn_solver/solver.py`, radiation/BEM modules, shared result schema.

## Milestone 3 — Verify acoustic predictions before optimization

**Work**

- Run existing tube and horn analytical tests through the production solver entry point. Tests that reimplement a miniature solver are supplementary, not substitutes.
- Compare matched boundary conditions and quantities. Separate exact analytical comparisons from approximate one-dimensional models and physical measurements.
- Run at least three mesh resolutions and two frequency-grid densities on the reference cases. Include geometry section-count convergence for lofted profiles.
- Verify frequency-band splitting and merging against the same unsplit sweep.
- Validate radiation against independent piston/open-duct references in their applicable regimes; validate absolute output and directivity separately.
- Acquire an independent numerical or measured reference with matching geometry, driver, interface, excitation, environment, and microphone position. Record source and reuse permission; do not label internally generated Webster data as measured evidence.

**Proposed initial gates**

- Exact straight-tube pressure response: maximum error at most 0.5 dB, retaining the existing target.
- Between the two finest meshes in the declared usable band: output change at most 0.5 dB and cutoff change at most 5%; investigate impedance/phase changes using separately recorded tolerances that remain meaningful near zeros.
- Split versus unsplit sweep: difference at most 0.5 dB on a common frequency grid.
- One-dimensional versus three-dimensional approximate cases: retain the existing 3 dB comparison target only within a documented regime where that approximation applies.

These are proposed engineering acceptance thresholds, not achieved results or universal acoustic standards. Freeze case-specific tolerances before evaluating release candidates. Changes need a written technical reason, rather than being adjusted merely to make a failure pass.

**Exit gate**

A versioned validation report includes expected versus computed curves, error statistics, mesh/frequency convergence, and limitations. Release fixtures satisfy their declared thresholds. Unresolved cases are excluded from the validated domain.

## Milestone 4 — Driver suitability and target-based acceptance

**Work**

- Curate a small initial database with parameter provenance, driver type, units, measured/recommended operating band, and uncertainty or missing-data flags. Review suspicious metadata edits.
- Separate known, derived, estimated, and missing parameters. Never interpret unknown power or excursion capability as unlimited capability.
- Specify the driver-to-throat interface, including compression ratio and any chamber/adapter model. Reject combinations outside that interface model's validated domain.
- Support band, permitted response deviation, drive conditions, maximum dimensions, and optional output/coverage requirements in the target specification.
- Evaluate ripple and output over the requested band. Detect internal response gaps and distinguish measured/simulated sweep bounds from identified cutoffs.
- Apply physical and data-quality eligibility rules before preference scoring. Report rejected candidates and reasons.
- Add excursion/thermal feasibility only where data and models support it. Otherwise explicitly withhold maximum-output claims.

**Exit gate**

Fixtures cover a suitable driver, an unsuitable upper-frequency range, an excessive interface compression ratio, missing essential data, an internal response notch, an impossible size/band combination, and no feasible driver. Every case has the expected eligibility result. A score cannot override a hard requirement.

**Primary files:** `horn_core/parameters.py`, driver loader/validator, prescreen, KPI, scoring, and rank pipeline modules.

## Milestone 5 — Demonstrate that the search finds good designs

**Work**

- Use validated analytical output for the coarse screen; retain profile/dimension diversity and uncertainty-driven candidates rather than relying on an untested top-three cutoff.
- Construct a tractable benchmark grid of roughly 20–50 geometries and 3–5 vetted drivers. Evaluate every geometry with FEM and couple every eligible driver to establish the best result within that finite grid.
- Separate calibration cases from held-out cases used to test screening quality. Compare feasible-design recall and final recommendation quality, not merely rank correlation.
- Add bounded refinement around promising mouth/throat/length values, with an explicit simulation budget, minimum improvement, and boundary checks.
- If a winner lies on a search bound, expand the search where allowed or report that the bound limits the result.
- Report several feasible tradeoffs and ties within numerical uncertainty. Describe the result as the best evaluated candidate, never a proven global optimum.
- Profile screening, coupling, and simulation costs. Set practical run budgets from measured performance and disclose budget exhaustion.

**Proposed initial gate**

On held-out finite-grid benchmarks, accelerated search retains the exhaustive FEM winner or returns a feasible design whose normalized preference score is within 0.02 of it. Score definition must be frozen first. At least 90% of the exhaustive top-ten feasible combinations remain representable by the shortlisted geometries. No infeasible combination is promoted. If these targets fail, increase screening breadth or revise the surrogate before shipping.

Refinement must improve or preserve the best feasible result, obey bounds, terminate within budget, and reproduce with a fixed seed. Increasing resolution must not conceal materially different recommendations.

**Primary files:** candidate/geometry designer, LEM prescreen, rank pipeline, `main.nf`, new search benchmark fixtures.

## Milestone 6 — Physical validation and prediction limits

**Work**

- Choose a readily characterized reference assembly in the supported domain, plus a comparison assembly. Obtain measured data or arrange a build and measurement session; purchasing or fabrication is a separate decision.
- Prepare mechanical drawings with known interface/chamber dimensions, not just an acoustic STEP volume.
- Freeze predicted response and ranking before measurements. Record actual drive voltage, impedance, microphone calibration/position, environment, geometry, and processing/gating limits.
- Measure electrical impedance and on-axis response; measure off-axis response where coverage is claimed. Repeat the setup to estimate measurement variability.
- Compare at a consistent distance and input level. Exclude frequency regions the measurement setup cannot resolve, and disclose the resulting validated band.
- Investigate systematic disagreement before calibrating. Any fitted correction must be evaluated on a separate held-out assembly or configuration.

**Proposed physical gate**

Across the declared measurement-valid band: median absolute level error at most 2 dB, 95th-percentile absolute error at most 4 dB, and comparable cutoff locations within 10%. No arbitrary level shift is permitted for an absolute-output claim. Major impedance/response features must be explained. Recommendation ordering must agree where predicted differences exceed combined uncertainty; otherwise report a tie.

These targets are provisional until the reference setup and uncertainty budget are specified. One matching prototype validates that case, not every driver or profile in the database.

**Exit gate**

Versioned raw measurements, setup description, frozen predictions, comparison results, and uncertainty notes establish the supported physical domain. No physical-validation claim is made if hardware/reference data remain unavailable.

## Milestone 7 — Release the band-to-design workflow

**Work**

- Preserve a simple entry point, for example `just run-auto --target_f_low 500 --target_f_high 4000`; expose all default assumptions in the result.
- Export the resolved specification, ranked feasible designs, dimensions and actual areas, driver identity/provenance, response/impedance/coverage where supported, rejection reasons, validation status, search bounds/budget, and acoustic CAD.
- Export source/data/image hashes and commands needed to reproduce the recommendation. Use isolated run directories to prevent stale reports from being mixed with current results.
- Update README and design documentation to match supported behavior and clearly identify experimental capabilities.
- Add scheduled/release numerical validation plus a fast CI smoke test. A release verification job fails if required cases are skipped or incomplete.

**Exit gate**

From a clean checkout, run feasible narrow-band and wider-band fixtures, a size-constrained case, an impossible request, a missing-data case, an injected solver failure, and an interrupted/resumed run. Each produces the expected status and complete provenance. Resume agrees with an uninterrupted run within declared numerical tolerances. The user can distinguish validated predictions, experimental predictions, and insufficient evidence without reading code.

## Execution order and progress tracking

Execute milestones 1 → 2 → 3 before tuning automatic optimization. Curated data collection for milestone 4 and reference-data preparation for milestone 6 can begin earlier. Milestone 5 depends on validated metrics and eligibility; milestone 7 consolidates the demonstrated behavior.

Track each milestone with implementation links, test commands, result artifacts, unresolved failures, and explicit pass/fail against its gate. Do not use a single percentage-complete estimate: working software, correct acoustic predictions, successful search, and physical agreement are separate achievements.

The first implementation batch is milestone 1 plus the area/units schema from milestone 2. Its deliverable is a reproducible, honestly labelled single reference run and a reliable miniature automatic run. It is not a new large driver sweep or an expansion of the custom bullet design.

This document is a proposed execution plan. No implementation changes or new physical validations are claimed by its creation.
