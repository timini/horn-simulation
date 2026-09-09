# Single-horn status and delivery roadmap

Updated 9 September 2026 following the independent solver, run-management and worked-example qualification work. This is the current delivery sequence; older plans and issue bodies retain historical evidence. [Issue #81](https://github.com/timini/horn-simulation/issues/81) remains the completion tracker.

## Intended result

Give the pipeline an operating band and optional size/output constraints. Receive ranked driver-and-horn candidates, dimensions, response and impedance predictions, acoustic CAD, rejection reasons, explicit assumptions and reproduction information. The first qualified release targets rigid, axisymmetric, front-loaded, single-driver horns with characterized interfaces. It must be able to return no feasible design or insufficient evidence.

Multiple-entry horns, interacting sources, entry taps and crossover optimisation belong in [MEH Design Studio](https://github.com/timini/meh-design-studio). Removing MEH does not remove the need to validate a conventional horn's driver interface and radiated output.

## What works today

Band-driven geometry generation, analytical screening, production FEM evaluation, driver coupling, ranking, bounded refinement, seven profiles, axial STEP import and HTML reports are implemented. The launcher isolates runs and preserves source, inputs, container identity and resume evidence. Invalid/incomplete frequency sweeps cannot produce a valid recommendation.

The independent solver change expands the required production acoustic lane from 17 to 23 cases. Ordinary analysis CI also checks six archived independent motor references, and core CI checks ten independent Rayleigh integrals. Current merge/check links are recorded in [issue #81](https://github.com/timini/horn-simulation/issues/81). The numerical environment has a documented x86 OpenBLAS workaround and geometry safeguards; its original intermittent failure's upstream cause remains unconfirmed.

| Evidence | Result | Limit |
| --- | --- | --- |
| Exhaustive finite-grid search | Four bands, 20 geometries and three motors; all feasible top-ten results and winners retained | Numerical search quality, not measured driver ranking or a global optimum |
| Independent numerical references | 40/40 agree at supplied frequencies | Some supplied curves are sparse; no continuous-band conclusion for those cases |
| Production numerical health / mesh convergence | 15/15 runs and 5/5 geometries pass | Restricted reference geometries and tested quantities |
| Dense measured pipe impedance | 237/299 curves within frozen magnitude limits | Repeated measurements, not 299 assemblies; failures retained |
| Matched independent interior FEM and ideal motor | Six tube/cone meshes, 13 frequencies each; maximum complex relative discrepancy below 9e-13 | Shared P1 meshes and plane-wave termination; ideal characterized motor only |
| Independent cavity modes | Twelve modes converge through three P2 meshes; finest maximum error 0.033% | Production volume operator; rigid box |
| Uniform baffled-piston radiation equations | Ten independent Rayleigh-integral checks pass | Infinite baffle and uniform piston, not a full exterior horn solve |
| Complete horn-and-driver validation | Outstanding | Absolute output and recommendation ordering are not established |

See the [reference validation report](REFERENCE_VALIDATION_REPORT.md), [independent solver evidence](INDEPENDENT_SOLVER_VALIDATION.md), [physical reference audit](PHYSICAL_REFERENCE_AUDIT.md) and [acoustic contract](ACOUSTIC_CONTRACT.md). The measured pipe evidence applies to the opt-in loss model and matched terminations. The default lossless workflow does not inherit these passes. The three curated manufacturer-parameter motors still lack verified interfaces and separate diaphragm/rear-load data.

## Generate an experimental candidate now

With the documented Docker, Nextflow, Java and `just` prerequisites installed, build the images and run, for example:

```sh
just build
just run-auto --target_f_low 800 --target_f_high 1600 --drivers_db data/drivers-curated --lem_top_n 3 --refinement_budget 2 --num_bands 2 --num_intervals 101 --mesh_size 0.01
```

This band is an example request, not a promise that the curated motors can satisfy it. Use `just latest-run` to locate the completed run, then open `outputs/auto/report/auto_report.html` inside that directory. Inspect the resolved assumptions and evidence status. See the [worked example](WORKED_EXAMPLE.md) and [candidate resolution procedure](CANDIDATE_RESOLUTION.md). It may report no feasible design or insufficient evidence. Generated dimensions and air-volume CAD are useful for exploration; they are not a physically qualified or manufacturing-ready assembly.

## Remaining issues

| Issue | Current purpose | Release role |
| --- | --- | --- |
| [#81](https://github.com/timini/horn-simulation/issues/81) | Complete assembly evidence, driver/interface qualification, radiation, convergence and final acceptance | Main completion gate |
| [#78](https://github.com/timini/horn-simulation/issues/78) | Complete physical interface and full exterior comparisons beyond the delivered interior-FEM/ideal-motor references | Supports physics qualification |
| [#74](https://github.com/timini/horn-simulation/issues/74) | Verified LaVoce refresh and provenance/metadata audit | Catalogue maintenance; a small qualified set can ship first |
| [#45](https://github.com/timini/horn-simulation/issues/45) | Interior field export for diagnosis | Optional unless needed to investigate a failure |
| [#49](https://github.com/timini/horn-simulation/issues/49) | Non-axial boundary identification | Deferred unless a chosen reference needs it |
| [#51](https://github.com/timini/horn-simulation/issues/51) | Conventional complex-geometry extensions | Deferred; MEH is excluded |

Closed/superseded work includes #41 in favour of #81, #77 in favour of #78, #47 for the documented axial STEP conventions, and #79 for the sister project plus the conventional radiation gates in #81. Run-management #75 is delivered in #85. LaVoce classification #73 is corrected in #86; the larger catalogue audit remains #74. Do not interpret closing these issues as validation of arbitrary geometry or exterior radiation. Historical #34/#48 closures do not qualify the disabled legacy FEM–BEM path.

## Ordered completion work

### 1. Select a tractable physical reference and comparison

Search existing published/open DIY data first. For each proposed assembly, inventory exact air and interface geometry, driver variant and mass convention, rear load, voltage, microphone distance/axis, environment, processing and reuse terms. Two sufficiently characterized assemblies/configurations are required to evaluate recommendation ordering, not just agreement for one example. Existing acquired DIY curves remain supporting material where these facts are missing.

**Exit:** a versioned reference specification, raw-data identities, declared usable measurement band and uncertainty, and frozen comparison tolerances. If evidence is insufficient, explicitly record the missing fields and the measurement/build procedure needed; do not substitute another unmatched response curve. Purchasing or fabrication is a separate decision.

### 2. Qualify the driver/interface model

Model the chosen chamber/adapter and rear load at the fidelity the reference requires. Independently test voltage scaling, area conversion, impedance signs, energy accounting and limiting cases. Distinguish diaphragm-only mass from mass containing an existing air load; do not invent the missing decomposition. Qualify a small driver set before expanding the scraped catalogue.

**Exit:** an independently implemented circuit/solver reference agrees under matched conventions, and unsupported combinations are rejected or flagged. The FSN classification correction is delivered through #73; continue source audits through #74 without treating a database refresh as physical validation.

### 3. Establish numerical and radiated-output accuracy

Under #78, the matched interior-FEM and ideal-motor comparisons are now executed and archived in this repository. They agree on six shared tube/cone meshes under pinned runtime and acoustic conventions. The [complete-horn radiation study](EXTERIOR_RADIATION_VALIDATION.md) now retains a matched infinite-baffle comparison and finite-baffle pilots. Its reference refinement checks pass, but the strict throat-impedance comparison fails (up to 1.05 dB / 5.73 degrees); pressure is within 0.40 dB. Investigate the mouth boundary/field approximation before promoting its horn loading. The separate uniform-piston integral checks do not establish a uniform horn-mouth velocity.

The independent rectangular-cavity production-operator check is delivered. Use the candidate-resolution harness to check mesh, frequency-grid and loft-section changes for each profile and range admitted to the release; one worked candidate does not qualify all seven profiles. Investigate reference failures or explicitly exclude unsupported materials/terminations; retain frozen tolerances. Absolute on-axis output needs its own evidence. Only claim coverage/directivity after separate validation.

**Exit:** versioned comparison artifacts pass case-specific limits, and the supported geometry, interface, radiation and frequency domain is explicit. Broader exterior BEM implementation is necessary only if the selected release domain cannot be supported by the validated approximation.

### 4. Validate complete assemblies and selection

Freeze predictions before comparison with qualifying measurements. Check absolute response without a fitted level shift, electrical impedance, cutoff and ordering between the reference and comparison assembly. Preserve separate calibration and held-out evidence if any fitted correction is introduced. Investigate disagreements rather than changing gates to obtain a pass.

**Exit:** both assemblies meet the documented physical gates over their measurement-valid band. Rank differences smaller than combined uncertainty become ties. Only then can the tested domain acquire a physically validated recommendation status.

### 5. Release a repeatable design workflow

Re-run finite-grid search checks with the qualified metrics and drivers, plus feasible, constrained, impossible, missing-data, failed-band and interrupted/resumed workflows. Publish a clean-checkout example with driver identity, horn dimensions, interface/rear-load specification, raw response evidence, limitations and reproduction commands. Direct Nextflow requires an explicit output directory; standard launchers isolate output and support completed-run lookup. Run-management #75 is delivered.

For a buildable reference, provide or link verified mechanical drawings for walls, mounting and the characterized interface; keep these distinct from the acoustic air-volume STEP. Automatic manufacturing optimisation is not required for this release.

**Exit:** a user can reproduce a qualified example and understand whether their new request falls inside the tested domain. Issue #81 closes only when its software and physical gates have linked evidence. New bands, drivers or profiles outside that domain remain experimental.

## Decision

Continue the existing implementation. The major remaining uncertainty is complete-assembly evidence and model qualification, not the existence of the automation. A percentage-complete estimate or fixed finish date would hide that dependency. Focus effort on the reference/interface/radiation sequence before adding more profiles, a larger catalogue or complex geometry.
