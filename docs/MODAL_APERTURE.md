# Experimental nonlocal aperture radiation

`modal_baffled` couples the three-dimensional interior pressure field to 16 axisymmetric aperture-velocity modes in an infinite rigid baffle. It replaces the local piston Robin boundary with a modal Rayleigh radiation impedance matrix and integrates the resulting nonuniform mouth velocity for the on-axis observer. It is optional; the default remains `flanged_piston`.

## Why this model exists

The [matched independent reference](EXTERIOR_RADIATION_VALIDATION.md) found that the local boundary's throat impedance missed the frozen 0.5 dB / 5 degree comparison limits, despite on-axis pressure agreeing within 0.40 dB. Those failed results remain in the repository. The new operator [passes the same reference](MODAL_APERTURE_VALIDATION.md): maximum differences are 0.017 dB for throat impedance and 0.013 dB for on-axis pressure, with phase errors below 0.19 degrees. Separate mode and mesh refinements also pass. This establishes one ideal geometry and band.

## Supported implementation and limits

- One connected, rigid, axisymmetric acoustic volume with circular disk inlet and mouth, centered on the z axis. The inlet is at z=0 and the mouth at z=length.
- CAD checks sample the port boundaries and compare the complete volume with two rotated copies. STEP spline uncertainty is allowed at 0.1 micrometre radial tolerance, with a corresponding thin-shell volume tolerance. Annular ports, offset ports, rectangular ports and non-axisymmetric internal cuts are rejected.
- Lossless air, P1 elements and a single MPI rank. Infinite baffle only. Finite-baffle diffraction, non-axisymmetric modes, directivity, wall losses and fabrication detail are outside this implementation.
- The quadrature accepts mouth ka up to 30; this is a numerical domain limit, **not an experimentally validated bandwidth**. Every new horn/band needs mode, mesh and frequency convergence checks. Sixteen modes do not guarantee convergence for arbitrary dimensions.
- The final auto-ranking uses the full modal on-axis pressure. Webster screening remains a plane-mode approximation and is labelled accordingly. The [finite-grid search audit](MODAL_SEARCH_VALIDATION.md) retains all exhaustive top-ten results and winners in four bands; this does not establish a global optimum or performance outside that grid.
- Driver and interface evidence rules remain in force. Numerical radiation agreement does not characterize a cone-to-throat chamber, phase plug, rear load or moving-mass convention. Results remain experimental predictions.

## Running it

After rebuilding all images, add `--radiation_model modal_baffled` to an auto-mode request. The normal report records that model and the `modal_baffled_on_axis` output metric. The response CSV preserves all complex aperture coefficients, their mode count, the interface residual and acoustic power. Incomplete modal coefficients and inconsistent frequency-band joins are rejected.

```sh
just build
just run-auto --target_f_low 800 --target_f_high 1600 \
  --drivers_db data/drivers-curated --lem_top_n 3 --refinement_budget 2 \
  --num_bands 2 --num_intervals 101 --mesh_size 0.01 \
  --radiation_model modal_baffled
```

This command generates an experimental candidate, not a qualified physical assembly. Single-mode mouth-pressure diagnostics retain their existing meaning; an on-axis auto report uses the nonuniform aperture observer.

## Executed automatic design

A fresh 800–1600 Hz search completed all 18 processes using this model; a subsequent resume reused all 18 from cache. The leading candidate is an 18Sound 6NMB420 with a hyperbolic air passage: throat diameter 64.996 mm, mouth diameter 177.418 mm and length 117.906 mm. Its predicted target-band ripple is **1.485 dB**, with mean output **104.882 dB SPL at 1 m and 2.83 V RMS**. The report retains a near tie and insufficient-evidence status because the physical driver interface and moving-mass/rear-load separation are unverified.

The [worked-example archive](../data/validation/worked_example_modal_800_1600.tar.gz) contains the self-contained report, STEP geometry, rankings, raw and coupled curves, resolved specification, sealed manifest and fresh/resume logs. The [identity manifest](../data/validation/worked_example_modal_800_1600_manifest.json) names the executed revision. The selected STEP is `outputs/auto/refinement/refine_hyperbolic_0001.step`; it represents the acoustic air volume, with no wall or mounting design.

A second [archived workflow](../data/validation/modal_empty_workflow.tar.gz) requests a 2 m mouth radius outside this model's ka domain. It completes with an empty ranking and no FEM tasks, explaining the rejection instead of failing inside the solver. The older default-model [worked example](WORKED_EXAMPLE.md) remains separate: changing the radiation assumption changes its predicted ripple from about 0.99 dB to 1.48 dB even though the winning dimensions remain the same.

The actual ranked candidate also passes all six [mesh, loft and frequency comparisons](MODAL_CANDIDATE_RESOLUTION.md). The largest mesh change in coupled output is 0.063 dB. This numerical stability does not resolve the missing physical interface evidence.

## Mathematical and numerical checks

Modes are `J0(mu_n r/a)/J0(mu_n)`, where `mu_0=0` and subsequent roots satisfy `J1(mu_n)=0`. Each has area norm `S`. The Fourier–Bessel spectrum is integrated over propagating and evanescent plane waves with positive quadrature weights; this preserves reciprocity and nonnegative radiation resistance. The finite evanescent tail is explicit and tested by doubling its range.

The FEM pressure and scaled modal velocity `w=rho*c*v` are solved together:

```text
A p + i k B w = b
B^T p - S Z w = 0
```

Here `B` is the aperture projection and `Z` is the specific modal impedance divided by `rho*c`. Pressure and velocity use the `exp(+iwt)` convention and RMS amplitudes. The implementation checks linear residual, interface closure and radiated power. Core tests compare the plane-mode impedance with the analytic baffled piston and the observer with its exact near/far-axis solution. Solver tests compare pressure and velocity excitation and their energy accounting.

The independent reference is the separately executed [MMM toolbox](https://github.com/bkolbrek/MMM_toolbox), pinned and archived as documented in the reference report. The general modal-radiation method is described in [Zorumski's NASA record](https://ntrs.nasa.gov/citations/19740038818) and [Kolbrek's thesis](https://kolbrek.hornspeakersystems.info/images/misc/KolbrekPhDThesisSubmitted.pdf). No upstream implementation is vendored into this solver.
