# Experimental nonlocal aperture radiation

`modal_baffled` couples the three-dimensional interior pressure field to 16 axisymmetric aperture-velocity modes in an infinite rigid baffle. It replaces the local piston Robin boundary with a modal Rayleigh radiation impedance matrix and integrates the resulting nonuniform mouth velocity for the on-axis observer. It is optional; the default remains `flanged_piston`.

## Why this model exists

The [matched independent reference](EXTERIOR_RADIATION_VALIDATION.md) found that the local boundary's throat impedance missed the frozen 0.5 dB / 5 degree comparison limits, despite on-axis pressure agreeing within 0.40 dB. Those failed results remain in the repository. The new operator is being qualified against the same converged reference, with additional mode and mesh refinements; it does not inherit a pass merely because the equations differ.

## Supported implementation and limits

- One connected, rigid, axisymmetric acoustic volume with circular disk inlet and mouth, centered on the z axis. The inlet is at z=0 and the mouth at z=length.
- CAD checks sample the port boundaries and compare the complete volume with two rotated copies. STEP spline uncertainty is allowed at 0.1 micrometre radial tolerance, with a corresponding thin-shell volume tolerance. Annular ports, offset ports, rectangular ports and non-axisymmetric internal cuts are rejected.
- Lossless air, P1 elements and a single MPI rank. Infinite baffle only. Finite-baffle diffraction, non-axisymmetric modes, directivity, wall losses and fabrication detail are outside this implementation.
- The quadrature accepts mouth ka up to 30; this is a numerical domain limit, **not an experimentally validated bandwidth**. Every new horn/band needs mode, mesh and frequency convergence checks. Sixteen modes do not guarantee convergence for arbitrary dimensions.
- The final auto-ranking uses the full modal on-axis pressure. Webster screening remains a plane-mode approximation and is labelled accordingly; search recall must be rechecked before claiming a qualified optimiser using this mode.
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

## Mathematical and numerical checks

Modes are `J0(mu_n r/a)/J0(mu_n)`, where `mu_0=0` and subsequent roots satisfy `J1(mu_n)=0`. Each has area norm `S`. The Fourier–Bessel spectrum is integrated over propagating and evanescent plane waves with positive quadrature weights; this preserves reciprocity and nonnegative radiation resistance. The finite evanescent tail is explicit and tested by doubling its range.

The FEM pressure and scaled modal velocity `w=rho*c*v` are solved together:

```text
A p + i k B w = b
B^T p - S Z w = 0
```

Here `B` is the aperture projection and `Z` is the specific modal impedance divided by `rho*c`. Pressure and velocity use the `exp(+iwt)` convention and RMS amplitudes. The implementation checks linear residual, interface closure and radiated power. Core tests compare the plane-mode impedance with the analytic baffled piston and the observer with its exact near/far-axis solution. Solver tests compare pressure and velocity excitation and their energy accounting.

The independent reference is the separately executed [MMM toolbox](https://github.com/bkolbrek/MMM_toolbox), pinned and archived as documented in the reference report. The general modal-radiation method is described in [Zorumski's NASA record](https://ntrs.nasa.gov/citations/19740038818) and [Kolbrek's thesis](https://kolbrek.hornspeakersystems.info/images/misc/KolbrekPhDThesisSubmitted.pdf). No upstream implementation is vendored into this solver.
