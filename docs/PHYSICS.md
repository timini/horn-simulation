# Physics

This page describes what the pipeline computes. [ACOUSTIC_CONTRACT](ACOUSTIC_CONTRACT.md) fixes the sign and unit conventions in detail, and [VALIDATION](VALIDATION.md) says which parts have been checked against what.

All phasors are RMS with time dependence exp(+iwt). Air is c = 343 m/s and rho = 1.225 kg/m3 unless a validation script states otherwise.

## Model in one paragraph

A three-dimensional finite element model solves the linear wave equation inside a rigid, axisymmetric air volume (the horn). The mouth radiates into an infinite rigid baffle. This gives the horn's throat impedance and the transfer from throat to mouth. A lumped driver model (Thiele-Small motor, optional front and rear chambers) is then attached to the throat impedance for each driver, which gives the diaphragm motion, electrical impedance and the on-axis pressure at the listening distance. Everything is linear. Nonlinear distortion, heating, wall vibration, finite baffles and cabinets are not modelled.

## Interior solve

The solver (`horn_solver.solver`) solves the Helmholtz equation for the complex pressure p,

```text
laplacian(p) + k^2 p = 0,   k = 2 pi f / c
```

on a tetrahedral mesh of the STEP volume with Lagrange elements of degree 1 (default) or 2, using FEniCSx/dolfinx. Each frequency is a direct LU solve through PETSc. A solve whose relative residual is above 1e-8 stops the run.

Boundaries are found from the CAD: the inlet at z = 0, the mouth at z = length, the rest is wall.

| Boundary | Condition |
| --- | --- |
| Inlet (throat) | p = 1 Pa. The inlet volume velocity comes from the boundary reaction, so the throat impedance is p / u. |
| Wall | Rigid (zero normal velocity). With `loss_model boundary_layer`, a thin viscous and thermal boundary layer term (Berggren, Bernland and Noreland 2018, eq. 31). |
| Mouth | The radiation model, below. |

Because the inlet is driven at 1 Pa, the solver output is a transfer function. The `spl` column is the RMS pressure over the mouth plane for 1 Pa at the throat. It is not the level at a listener. Listener levels come from the coupling step.

## Mouth radiation models

| Model | What it is | When it is used |
| --- | --- | --- |
| `modal_baffled` | Couples the interior field to 16 axisymmetric mouth velocity modes and a Rayleigh radiation impedance matrix for a circular aperture in an infinite rigid baffle. The mouth velocity can vary across the aperture. | Default for generated geometry with lossless P1 elements. See [MODAL_APERTURE](MODAL_APERTURE.md). |
| `flanged_piston` | Local Robin condition with the impedance of a rigid piston in an infinite baffle, assuming uniform mouth velocity. | Default for imported STEP, `boundary_layer` losses or P2 elements. Known throat impedance error up to 1.05 dB and 5.7 degrees against the matched reference; reports print this. |
| `unflanged_piston` | Local condition for an unflanged pipe end. | Valid only below mouth ka = 1.5; candidates beyond are rejected. |
| `finite_flange` | Local condition for a pipe with a flange of `flange_width`. | Opt-in. |
| `plane_wave` | rho c termination (no reflection for plane waves). | Opt-in, for checks. |
| `closed` | Rigid end. | Single mode only, for impedance checks. |

The local models enter the weak form as a Robin term `+ i k / z` on the mouth, with z the normalised specific radiation impedance. The old FEM-BEM coupling to bempp-cl has been removed from the solver: it coupled the whole boundary rather than the mouth.

## Meshing

Each frequency band gets its own mesh from gmsh, sized for that band's highest frequency (`horn_solver.meshing`):

- **Element size.** P1 uses lambda/10, P2 lambda/6. The phase speed error is about 1.6 % per wavelength for P1 at lambda/10, 4.6 % for P1 at lambda/6 and 0.08 % for P2 at lambda/6. The size used is the finer of this rule and `mesh_size`.
- **Throat.** Elements are capped at r_t/5 within one throat radius of the throat, grading to the bulk size over two more. The meshed throat area must be within 1 % of the CAD area, because the throat flux is integrated over the mesh. One retry at r_t/10, then the candidate is rejected (`throat_area_guard_failed`).
- **Node limit.** If a P1 mesh would exceed `max_mesh_nodes` (default 100,000), it is coarsened just enough to fit, never beyond lambda/6. Each row records `mesh_elements_per_wavelength` and `mesh_coarsened_for_node_limit`. If lambda/6 does not fit, an auto-mode candidate is rejected (`mesh_node_limit_exceeded`) and the run continues.

Example: at 8 kHz, lambda is 43 mm, so the P1 target is 4.3 mm. The default single-mode horn (conical, 50/200/500 mm) needs about 266,000 nodes at that size, so its top band is coarsened to about lambda/6.9; lower bands use lambda/10.

## Frequency sweep

Frequencies are spaced evenly in log frequency. The sweep is split into `num_bands` bands of equal log width, each with the same number of points: at least `num_intervals / num_bands`, raised to meet `min_points_per_octave`. Bands share their end points, so the merged grid is uniform. Bands are solved in parallel and merged; a missing band, missing sample or an inconsistent join stops the run.

Single mode sweeps `min_freq` to `max_freq`. Auto mode sweeps from one octave below `target_f_low` to half an octave above `target_f_high`, so that f3 points outside the band can be found.

## Driver coupling

`horn_core.driver_interface` and `horn_analysis.transfer_function` connect the driver to the solved throat impedance:

```text
V_g -- Re + i w Le -- Bl -- Rms + i w Mms + 1/(i w Cms)
     -- optional sealed rear chamber (air spring)
     -- area Sd
     -- optional front chamber (shunt compliance, or a short duct)
     -- horn throat impedance from the FEM solve
```

This gives the diaphragm velocity, excursion, electrical impedance, input power and the throat volume velocity at `voltage_rms`. Because the horn is linear, the FEM pressure field scales with that throat drive.

## Listener level

The level used for scoring, ranking and the coupled CSVs is the on-axis pressure at `observation_distance` from the mouth plane, in an infinite baffle:

- with `modal_baffled`, from the 16 mouth velocity modes by a Rayleigh integral (`modal_baffled_on_axis`);
- with the local models, from the mouth volume velocity treated as a uniform baffled piston (`uniform_baffled_piston_on_axis`).

Ranking rows record which was used in `output_metric`. The default distance is 1 m and the default drive 2.83 V RMS.

## Directivity

For `modal_baffled` horns the auto report adds an experimental far-field directivity, computed from the same 16 mouth modes in an infinite baffle (`horn_analysis.directivity`). It is relative to the axis, so it does not depend on the driver. It has not been compared with a measured polar. A bempp-cl comparison with a finite baffle differs, mostly from baffle diffraction; see [CROSS_SOFTWARE_VALIDATION](CROSS_SOFTWARE_VALIDATION.md).

## Analytical screen

The auto screen (`horn_core.webster`) slices the horn into thin cylinders and cascades plane-wave transfer matrices from throat to mouth, with a mouth radiation impedance. It takes milliseconds per geometry and is used only to choose which geometries to solve with FEM. It is compared with openwind in [CROSS_SOFTWARE_VALIDATION](CROSS_SOFTWARE_VALIDATION.md).

## Profiles

All profiles come from one module, `horn_core.profiles`, so the screen, the STEP, the render and the builder CSV use the same curve.

| Profile | Radius r(z) |
| --- | --- |
| `conical` | Straight line from throat to mouth radius. |
| `exponential` | r_t exp(z / x0). |
| `hyperbolic` | Salmon hypex with T = 0 (catenoidal), r_t cosh(z / x0). Other T values only through the Python API and `horn-geometry --hypex-t`. |
| `tractrix` | True tractrix. Its natural length has tractrix constant equal to the mouth radius. A longer length uses a larger constant; a shorter one is impossible and is rejected. |
| `os` | Oblate spheroidal (Geddes), sqrt(r_t^2 + z^2 tan^2(theta)). |
| `scaled_tractrix` | Tractrix contour with constant equal to the mouth radius, scaled axially to the length. Formerly called `lecleach`; it was never a Le Cleac'h horn. |
| `cd` | Exponential throat joined smoothly to a tangent cone. |

The STEP generator revolves the profile about the z axis and checks the result against the profile by dense sampling. It must stay within 1 %; the generator adds profile points beyond `num_sections` when needed.

## Limits of the model

- Rigid, axisymmetric air volume with circular inlet and mouth. Non-axisymmetric or rectangular horns are rejected by `modal_baffled` and are not generated.
- Infinite rigid baffle at the mouth. A free-standing horn or a cabinet radiates differently.
- Linear driver model. No distortion, thermal compression or large-signal behaviour. Excursion and power checks are linear rejections, not a maximum-output prediction.
- Driver interfaces (adapters, phase plugs) are lumped approximations or estimates; see [DRIVERS](DRIVERS.md).
- Wall losses are opt-in and checked only in a restricted pipe regime.
