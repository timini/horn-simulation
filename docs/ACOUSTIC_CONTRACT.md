# Acoustic contract, version 2

The supported automatic path is an experimental, linear, rigid-wall interior FEM model with a local mouth impedance and a uniform circular baffled-piston observation approximation. It does not yet establish the accuracy of a free-standing horn, its complete exterior field, or a particular commercial driver's phase plug.

## Quantities and conventions

- Complex amplitudes are RMS, with time dependence `exp(+iωt)`. Pressure: Pa; particle velocity: m/s; volume velocity: m³/s.
- Inlet velocity is positive into the horn; mouth velocity is positive outward. Specific acoustic impedance is `p/v`, in Pa·s/m. Acoustic impedance is `p/U = Z_specific / area`.
- Solver excitation is 1 Pa at the inlet. `spl` is the spatial RMS **mouth-plane pressure level**, not listening-distance sensitivity. `mouth_p_real/imag` is the complex area-mean mouth pressure; `mouth_u_real/imag` is integrated outward volume velocity.
- `inlet_area_m2` and `mouth_area_m2` are physical areas of the tagged CAD boundaries, including annular openings. Separate `mesh_*_area_m2` fields preserve the triangulated areas. FEM flux is integrated on the mesh; specific impedance is normalized by the physical inlet area so the conversion to `p/U` is consistent across meshes. Mesh-only API callers use their tagged mesh areas.
- Local normalized radiation impedance is `z = Z_specific/(ρc)`. The FEM Robin term is `+ik/z`, not `+ik*z`. With the adopted convention, pressure propagates as `exp(-ikz)` in a matched tube.
- Dirichlet inlet volume flux is recovered from the variational boundary reaction. Averaged linear-element boundary gradients produced large impedance errors near resonances and are not used for the local-model inlet transfer.

## Driver and observer

The motor equations are `V = Ze I + BL v` and `BL I = (Zm + Z_specific Sd²/Sthroat) v`. A force-balance test checks these independently and checks electrical input against resistive, mechanical and acoustic power. Defaults are 2.83 V RMS and 1 m from the mouth plane. **2.83 V is only 1 W for a purely resistive 8 Ω load.**

Both the analytical screen and final FEM ranking use mouth volume velocity, then the on-axis Rayleigh result for a uniformly vibrating circular piston:

`p(r) = ρ c (U/A) [exp(-ikr) − exp(-ik√(r²+A/π))]`.

This establishes a common, dimensionally defined observable. Applying it to a nonuniform horn aperture or an unbaffled assembly remains an approximation. It is not an exterior-wave solver or a directivity prediction. See the [COMSOL baffled-piston reference](https://doc.comsol.com/6.3/doc/com.comsol.help.models.aco.baffled_piston_radiation/baffled_piston_radiation.html) and [impedance boundary definition](https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.023.html).

When available, diaphragm-only mass and rear air-load mass are separate inputs. Legacy Mms coupling remains experimental: no guessed air-mass subtraction is applied. Missing driver operating band, parameter provenance, interface model, separated mass or maximum-output data remain evidence gaps. A feasible model response is not a verified driver recommendation.

## Acceptance and numerical evidence

Target bounds must be finite and increasing. All target-band samples, including interpolated endpoints, determine ripple and coverage; out-of-band peaks and internal notches cannot be hidden by outer cutoff estimates. Hard failures receive zero preference score and explicit rejection reasons. Missing evidence receives `insufficient_evidence`.

Schema-2 band merging requires every expected band, frequency sample, convention and area. It rejects missing/nonfinite data and mouth-level discontinuities greater than 0.5 dB at overlapping endpoints. The production tube tests additionally check complex impedance, phase, passive loading, three mesh resolutions and split/unsplit sweeps. A finite 20-geometry benchmark compares screening against exhaustive FEM using synthetic motors; it does not validate real drivers.

The optional BEM path remains experimental and is excluded from automatic optimization. Its current coupling code uses a whole-boundary trace rather than a verified aperture-only exterior model; backend sphere tests must not be mistaken for validation of that horn coupling. The [Bempp operator definitions](https://bempp.com/handbook/api/boundary_operators.html) define a mass matrix for the weak identity operator, not a coefficient-space identity.


## Losses and numerical acceptance

Optional boundary-layer wall loss, finite-flange pipe termination and quadratic elements are documented with their restricted domains in [loss validation](LOSS_PHYSICS_VALIDATION.md). Every production frequency solve rejects a nonfinite or relative residual above 1e-8 and a nonpositive solver convergence reason. Side-wall viscous/thermal and inlet/outlet RMS powers are exported separately. Direct motor coupling requires the default air properties. The legacy FEM–BEM entry points are disabled because their boundary mapping was invalid.


## Frequency-band joins

Adjacent endpoint samples must agree in mouth pressure level within 0.5 dB and in complex throat impedance, mouth pressure and mouth volume velocity. The complex tolerance is 5% of the larger endpoint magnitude plus 1% of a fixed unit-pressure reference scale (ρc for specific impedance, 1 Pa for mouth pressure, mouth area/(ρc) for volume velocity). The fixed floor keeps comparisons meaningful near zeros without allowing phase changes to hide behind equal SPL. A failed join requires mesh refinement; no smoothing or arbitrary choice of the duplicate bypasses it.

A known lower or upper driver-band limit is enforced even if the other endpoint is missing. Custom throat fractions must be finite and within (0, 1], and quantization cannot exceed the declared ka cap. Positive size caps that leave no geometry in the current heuristic range produce an explicit empty report; contradictory or nonpositive bounds remain invalid inputs.

Explicit throat radii obey the same high-frequency ka cap as screened radii. Driver fit is checked per candidate using effective diaphragm area versus actual mouth area. Single-mode ripple retains original response extrema; stitching errors must fail merge validation instead of being smoothed away. Both coupling paths use characterized diaphragm plus rear-load mass when available, with the same disclosed Mms fallback otherwise.

Imported STEP reports use tagged CAD inlet/mouth areas instead of parametric radius defaults. Variable-throat searches display the range and each result’s throat radius. Relative-shape validation has a separate shape gate and never passes the absolute-level gate. Run manifests record the resolved Nextflow launcher hash and engine version/build; the selected engine version is pinned for execution and an incompatible or older unrecorded manifest cannot resume.
