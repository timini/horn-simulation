# Single-horn status and delivery roadmap

Reviewed 9 September 2026 against main `eacb273`. This is the current delivery sequence; older plans and issue bodies retain historical evidence. [Issue #81](https://github.com/timini/horn-simulation/issues/81) remains the completion tracker.

## Intended result

Give the pipeline an operating band and optional size/output constraints. Receive ranked driver-and-horn candidates, dimensions, response and impedance predictions, acoustic CAD, rejection reasons, explicit assumptions and reproduction information. The first qualified release targets rigid, axisymmetric, front-loaded, single-driver horns with characterized interfaces. It must be able to return no feasible design or insufficient evidence.

Multiple-entry horns, interacting sources, entry taps and crossover optimisation belong in [MEH Design Studio](https://github.com/timini/meh-design-studio). Removing MEH does not remove the need to validate a conventional horn's driver interface and radiated output.

## What works today

Band-driven geometry generation, analytical screening, production FEM evaluation, driver coupling, ranking, bounded refinement, seven profiles, axial STEP import and HTML reports are implemented. The launcher isolates runs and preserves source, inputs, container identity and resume evidence. Invalid/incomplete frequency sweeps cannot produce a valid recommendation.

The latest reviewed main revision passed [all seven CI jobs](https://github.com/timini/horn-simulation/actions/runs/34235435571) and [17 required production acoustic checks](https://github.com/timini/horn-simulation/actions/runs/34235435608). The numerical environment has a documented x86 OpenBLAS workaround and geometry safeguards; its original intermittent failure's upstream cause remains unconfirmed.

| Evidence | Result | Limit |
| --- | --- | --- |
| Exhaustive finite-grid search | Four bands, 20 geometries and three motors; all feasible top-ten results and winners retained | Numerical search quality, not measured driver ranking or a global optimum |
| Independent numerical references | 40/40 agree at supplied frequencies | Some supplied curves are sparse; no continuous-band conclusion for those cases |
| Production numerical health / mesh convergence | 15/15 runs and 5/5 geometries pass | Restricted reference geometries and tested quantities |
| Dense measured pipe impedance | 237/299 curves within frozen magnitude limits | Repeated measurements, not 299 assemblies; failures retained |
| Complete horn-and-driver validation | Outstanding | Absolute output and recommendation ordering are not established |

See the [reference validation report](REFERENCE_VALIDATION_REPORT.md) and [acoustic contract](ACOUSTIC_CONTRACT.md). The measured pipe evidence applies to the opt-in loss model and matched terminations. The default lossless workflow does not inherit these passes. The three curated manufacturer-parameter motors still lack verified interfaces and separate diaphragm/rear-load data.

## Generate an experimental candidate now

With the documented Docker, Nextflow, Java and `just` prerequisites installed, build the images and run, for example:

```sh
just build
just run-auto --target_f_low 500 --target_f_high 4000 --drivers_db data/drivers-curated
```

This band is an example request, not a promise that the curated motors can satisfy it. Inspect the new `results/<run-id>/` report, resolved assumptions and evidence status. It may report no feasible design or insufficient evidence. Generated dimensions and air-volume CAD are useful for exploration; they are not a physically qualified or manufacturing-ready assembly.

## Remaining issues

| Issue | Current purpose | Release role |
| --- | --- | --- |
| [#81](https://github.com/timini/horn-simulation/issues/81) | Complete assembly evidence, driver/interface qualification, radiation, convergence and final acceptance | Main completion gate |
| [#78](https://github.com/timini/horn-simulation/issues/78) | Execute matched independent solver/coupling comparisons; retain the earlier feature comparison as history | Supports physics qualification |
| [#73](https://github.com/timini/horn-simulation/issues/73) | Fix LaVoce FSN classification still wrong on reviewed main | Correctness fix for the wider catalogue |
| [#74](https://github.com/timini/horn-simulation/issues/74) | Verified LaVoce refresh and provenance/metadata audit | Catalogue maintenance; a small qualified set can ship first |
| [#75](https://github.com/timini/horn-simulation/issues/75) | Decide safe direct-Nextflow output behavior and finish latest-run convenience | Usability; standard launcher already isolates runs |
| [#45](https://github.com/timini/horn-simulation/issues/45) | Interior field export for diagnosis | Optional unless needed to investigate a failure |
| [#49](https://github.com/timini/horn-simulation/issues/49) | Non-axial boundary identification | Deferred unless a chosen reference needs it |
| [#51](https://github.com/timini/horn-simulation/issues/51) | Conventional complex-geometry extensions | Deferred; MEH is excluded |

Retire #41 in favour of #81, #77 in favour of #78, #47 as delivered within the documented axial STEP conventions, and #79 in favour of the sister project plus the conventional radiation gates in #81. Do not interpret closing these issues as validation of arbitrary geometry or exterior radiation. Historical #34/#48 closures do not qualify the disabled legacy FEM–BEM path.

## Ordered completion work

### 1. Select a tractable physical reference and comparison

Search existing published/open DIY data first. For each proposed assembly, inventory exact air and interface geometry, driver variant and mass convention, rear load, voltage, microphone distance/axis, environment, processing and reuse terms. Two sufficiently characterized assemblies/configurations are required to evaluate recommendation ordering, not just agreement for one example. Existing acquired DIY curves remain supporting material where these facts are missing.

**Exit:** a versioned reference specification, raw-data identities, declared usable measurement band and uncertainty, and frozen comparison tolerances. If evidence is insufficient, explicitly record the missing fields and the measurement/build procedure needed; do not substitute another unmatched response curve. Purchasing or fabrication is a separate decision.

### 2. Qualify the driver/interface model

Model the chosen chamber/adapter and rear load at the fidelity the reference requires. Independently test voltage scaling, area conversion, impedance signs, energy accounting and limiting cases. Distinguish diaphragm-only mass from mass containing an existing air load; do not invent the missing decomposition. Qualify a small driver set before expanding the scraped catalogue.

**Exit:** an independently implemented circuit/solver reference agrees under matched conventions, and unsupported combinations are rejected or flagged. Correct #73 and audit source records through #74 without treating a database refresh as physical validation.

### 3. Establish numerical and radiated-output accuracy

Under #78, run the same simple geometry in this solver and a pinned Boundary Lab installation with matching inlet, termination, air, frequency grid and phasor convention. Compare complex impedance and transfer quantities before adding exterior-model differences. Then compare the supported mouth/observer approximation with an independent exterior solution in its intended regime. The MEH repository's executed Boundary Lab adapter is useful integration evidence, not a completed comparison here.

Add an independent rectangular-cavity assembly check where appropriate. Complete mesh, frequency-grid and loft-section convergence for the profiles admitted to the release. Investigate reference failures or explicitly exclude unsupported materials/terminations; retain frozen tolerances. Absolute on-axis output needs its own evidence. Only claim coverage/directivity after separate validation.

**Exit:** versioned comparison artifacts pass case-specific limits, and the supported geometry, interface, radiation and frequency domain is explicit. Broader exterior BEM implementation is necessary only if the selected release domain cannot be supported by the validated approximation.

### 4. Validate complete assemblies and selection

Freeze predictions before comparison with qualifying measurements. Check absolute response without a fitted level shift, electrical impedance, cutoff and ordering between the reference and comparison assembly. Preserve separate calibration and held-out evidence if any fitted correction is introduced. Investigate disagreements rather than changing gates to obtain a pass.

**Exit:** both assemblies meet the documented physical gates over their measurement-valid band. Rank differences smaller than combined uncertainty become ties. Only then can the tested domain acquire a physically validated recommendation status.

### 5. Release a repeatable design workflow

Re-run finite-grid search checks with the qualified metrics and drivers, plus feasible, constrained, impossible, missing-data, failed-band and interrupted/resumed workflows. Publish a clean-checkout example with driver identity, horn dimensions, interface/rear-load specification, raw response evidence, limitations and reproduction commands. Finish #75's remaining output behavior explicitly.

For a buildable reference, provide or link verified mechanical drawings for walls, mounting and the characterized interface; keep these distinct from the acoustic air-volume STEP. Automatic manufacturing optimisation is not required for this release.

**Exit:** a user can reproduce a qualified example and understand whether their new request falls inside the tested domain. Issue #81 closes only when its software and physical gates have linked evidence. New bands, drivers or profiles outside that domain remain experimental.

## Decision

Continue the existing implementation. The major remaining uncertainty is complete-assembly evidence and model qualification, not the existence of the automation. A percentage-complete estimate or fixed finish date would hide that dependency. Focus effort on the reference/interface/radiation sequence before adding more profiles, a larger catalogue or complex geometry.
