# Cross-software validation

Until now the solver was checked mostly against itself, analytic formulas and one in-house run of the MMM toolbox. This page compares it with three established tools that were written independently of this project:

- [openwind](https://openwind.inria.fr/) 0.12.4 (Inria, GPL-3): one-dimensional finite elements for the input impedance of wind instruments.
- [bempp-cl](https://bempp.com/) 0.4.2 (MIT): three-dimensional boundary elements for exterior radiation, numba backend.
- [Hornresp](https://www.hornresp.net/) (freeware, Windows): a widely used horn loudspeaker simulator. It could not be run here (no Wine), so input sheets and a comparison script are provided for a manual run.

Neither openwind nor bempp-cl is a dependency of the project packages. They run in a separate environment from `validation/cross_software/requirements.txt`.

The openwind, bempp-cl and Hornresp tolerances were committed (799c0e2) before the first comparison was run, and none was changed afterwards. Checks added later each had their tolerance fixed in the script before that check first ran:

- the FEM tolerance (part A2) was written before `fem_compare.py` was first run, but committed (a7de55e, 16:08 UTC) after the conical FEM run had finished (16:01 UTC) and before the exponential run finished (16:14 UTC);
- the Rayleigh quadrature and the spherical-wave 1D reference;
- in the review fix pass: openwind's `planar_piston` low-ka check, the FEM mesh refinement criterion and the like-for-like Hornresp SPL reference.

Failures are listed as failures.

## Summary

| Comparison | What it tests | Result |
| --- | --- | --- |
| Webster TMM against openwind, same physics, 7 profiles, 5 cases each | our plane-wave propagation, Kirchhoff and thin-layer losses, piston load, closed-end resonances | Pass, 35 of 35 |
| Webster TMM (exact piston load) against openwind's flanged-pipe radiation fit, ka at most 3 | the radiation model choice | Fails 6 of 7 profiles; expected model difference, see below |
| horn-solver FEM (3D) against openwind and TMM (1D), conical and exponential | 3D against 1D throat impedance with the same piston mouth impedance | Fails the frozen low-ka limit (0.60 dB and 0.40 dB median). The offset below ka 0.5 is wavefront curvature; the larger differences between ka 0.5 and 1.5 are not explained |
| bempp-cl pulsating sphere | the BEM set-up | Pass |
| Project modal far field against direct Rayleigh quadrature | the far-field transform in the experimental directivity, for the solver's mode basis | Pass, below 1e-12 dB |
| bempp-cl piston in a finite baffle against the infinite-baffle piston | piston radiation impedance and directivity | Impedance passes at ka 2 and 4; directivity fails at ka 2 because of the finite baffle |
| bempp-cl horn in a finite baffle against the archived modal_baffled FEM | throat impedance, directivity and on-axis level of a complete horn | Throat impedance passes at 5 of 5 frequencies; directivity and on-axis level fail, dominated by finite-baffle diffraction |
| BEM mesh refinement | BEM discretisation | 7 of 9 pairs pass; two directivity changes of 0.10 and 0.16 dB |
| Hornresp | plane-wave horn impedance and the 800 to 1600 Hz driver example | Not run; kit in `validation/hornresp/` |

No bug in the project code was found. Wherever the two codes model the same physics the TMM matches openwind within 0.07 dB at the 95th percentile, and within 0.43 dB at worst on the flank of a sharp closed-end resonance. Most remaining differences come from modelling choices the project already documents: the uniform piston against a flanged pipe, plane against curved wavefronts at low ka, and an infinite against a finite baffle. One is not explained: the 3D FEM with the local `flanged_piston` mouth load differs from the 1D codes by up to about 1.5 dB between mouth ka 0.5 and 1.5 (part A2).

## Reproduce

```sh
just validate-cross-software          # creates .venv-cross-software, then openwind, bempp-cl and the Hornresp sheets
just validate-cross-software-fem results/<conical-run> results/<exponential-run>
```

The FEM runs come from the normal launcher:

```sh
just run --profile conical --throat_radius 0.025 --mouth_radius 0.15 --length 0.3 \
  --min_freq 50 --max_freq 5000 --radiation_model flanged_piston --loss_model lossless --num_bands 2
```

and the same with `--profile exponential`. The runs used for the results here were `results/20261116T154534Z-c3ed006e` (conical, 161 frequencies, up to 123,803 cells) and `results/20261116T160638Z-1f441625` (exponential, up to 85,504 cells), both with horn-solver image `sha256:cfa7a18f6f9c` (tagged `horn-solver:latest` on 2026-11-16), from source revisions `0303d43` (conical) and `799c0e2` (exponential); the solver code is the same in both. Their run directories are not committed; the comparison JSON keeps the numbers, and each result records the source revision, container image ids and mesh settings from the run manifest.

Run time on an Apple laptop: openwind about 90 minutes (1,981 frequencies, 7 profiles, 7 solves each, dense bores), bempp-cl about 40 minutes.

Results: [openwind JSON](../data/validation/cross_software_openwind.json), [FEM JSON](../data/validation/cross_software_fem.json), [bempp-cl JSON](../data/validation/cross_software_bempp.json), plots in [docs/validation/cross_software](validation/cross_software/).

## A. openwind: one-dimensional input impedance

### Set-up

Throat radius 25 mm, mouth radius 150 mm, length 300 mm, for all seven profiles from `horn_core.profiles`. The tractrix uses its natural length for these radii, 223.8 mm. Frequencies 50 Hz to 5 kHz in 2 Hz steps (1,981 points; mouth ka from 0.14 to 13.7).

openwind is given the conical and exponential horns as exact segments. The other profiles are sampled at 801 points and given as conical pieces, which openwind solves exactly piece by piece. openwind fixes its air from temperature (20 C here: c = 343.99 m/s, rho = 1.1993 kg/m3, with its viscosity, conductivity, heat capacity and gamma); those values are passed to our code so the gas model is not part of the difference. Both codes use exp(+i w t). The quantity compared is Z/Zc at the throat, Zc = rho c / S_throat.

Our side is `compute_horn_transfer_tmm` with 200 segments (the prescreen default) and, as a convergence check, 2,000.

| Case | Our model | openwind model |
| --- | --- | --- |
| matched_lossless | `plane_wave` load (Z = rho c), lossless | `total_transmission`, lossless |
| closed_bessel | `closed`, `kirchhoff` losses | `closed`, `bessel` (Zwikker-Kosten) losses |
| piston_lossless | `flanged_piston` (exact Struve piston), lossless | the same piston impedance supplied as `from_data`, lossless |
| piston_bessel | `flanged_piston`, `kirchhoff` | `from_data` piston, `bessel` |
| piston_thin_layer | `flanged_piston`, `boundary_layer` | `from_data` piston, `wl` (Webster-Lokshin thin layer) |
| flanged_fit_bessel | `flanged_piston`, `kirchhoff` | `flanged_non_causal` (Silva et al. fit to the exact flanged pipe), `bessel` |

The first five cases share the physics, so near agreement is expected. The last uses a different radiation model on each side.

### Tolerances

| Class | Magnitude (median / p95 / max) | Phase (median / p95 / max) | Peak frequency | Band |
| --- | --- | --- | --- | --- |
| Same physics | 0.1 / 0.5 / 2.0 dB | 1 / 5 / 20 degrees | 0.5 % | 50 Hz to 5 kHz |
| Radiation model | 0.25 / 1.0 / 2.0 dB | 2 / 8 / 20 degrees | 1 % | ka at most 3 |

Peaks are maxima of |Z| with at least 1 dB prominence, refined by a parabola. The same number of peaks is also required.

### Results (200-segment TMM)

| Profile | matched lossless, dB max | closed + losses, dB p95 / max | piston lossless, dB max | piston + Bessel, dB max | piston + thin layer, dB max | largest peak error, same physics |
| --- | --- | --- | --- | --- | --- | --- |
| conical | 0.002 | 0.054 / 0.429 | 0.002 | 0.002 | 0.058 | 0.08 % |
| exponential | 0.003 | 0.019 / 0.131 | 0.003 | 0.003 | 0.070 | 0.07 % |
| hyperbolic | 0.005 | 0.015 / 0.044 | 0.005 | 0.005 | 0.083 | 0.05 % |
| tractrix | 0.011 | 0.043 / 0.238 | 0.011 | 0.011 | 0.072 | 0.02 % |
| os | 0.002 | 0.034 / 0.123 | 0.003 | 0.002 | 0.072 | 0.08 % |
| scaled_tractrix | 0.017 | 0.070 / 0.291 | 0.018 | 0.017 | 0.082 | 0.07 % |
| cd | 0.002 | 0.031 / 0.251 | 0.002 | 0.002 | 0.066 | 0.08 % |

All 35 same-physics comparisons pass at 200 and at 2,000 segments. Phase p95 is at most 0.25 degrees. The closed-end case has the sharpest resonances; its largest dB differences sit on the resonance and anti-resonance flanks, where a 0.01 % frequency shift gives a few tenths of a dB. The thin-layer case differs slightly more than the Bessel case because our `boundary_layer` and openwind's `wl` are different first-order expansions of the same Kirchhoff theory. openwind's own FEM and TMM agree with each other to better than 0.0001 dB on every bore, so the reference is internally consistent.

With a baffled mouth this horn's input impedance changes by only about 0.06 dB when losses are added, so the loss formulas are tested mainly by the closed-end case, where they set the resonance heights.

### Radiation model case

| Profile | dB median (limit 0.25) | dB p95 | Phase p95 | Peaks, openwind and ours (Hz) | Result |
| --- | --- | --- | --- | --- | --- |
| conical | 0.256 | 0.652 | 3.9 | 345 / 338, 776 / 753 | Fail (median, peaks) |
| exponential | 0.295 | 0.760 | 5.2 | 434 / 423, 715 / 692 | Fail (median, phase median, peaks) |
| hyperbolic | 0.258 | 0.848 | 5.5 | 496 / 491 | Fail (median) |
| tractrix | 0.146 | 0.777 | 5.2 | 618 / 617 | Pass |
| os | 0.223 | 0.656 | 4.2 | 380 / 370, 799 / 776 | Fail (peaks) |
| scaled_tractrix | 0.286 | 0.853 | 5.7 | 447 / 440, 1002 / none | Fail (median, phase median, peaks) |
| cd | 0.240 | 0.670 | 4.2 | 368 / 359, 792 / 769 | Fail (peaks) |

Our `flanged_piston` and `modal_baffled` prescreen load is a rigid piston in an infinite baffle. openwind's `flanged_non_causal` is a fit to the exact flanged pipe, where the mouth velocity is not uniform. Their low-frequency end corrections differ (0.85 a for the piston, 0.82 a for the flanged pipe), so with our load the peaks sit 2 to 3 % lower. This is a known difference between two radiation models, not a coding error: when openwind is given our piston impedance, the curves agree to 0.02 dB. The project's main ranking uses the modal FEM aperture (`modal_baffled`), which already passes an independent flanged-mouth reference (see [MODAL_APERTURE.md](MODAL_APERTURE.md)). The result is a reminder that the piston prescreen moves Webster peaks by a few percent.

## A2. horn-solver FEM against openwind and the TMM

The production FEM solves the 3D air volume. It was run in single mode with the `flanged_piston` mouth load and no losses. All three codes use the exact baffled-piston impedance Z_piston(ka) at a flat mouth disc, but they apply it differently, and the throat is driven differently:

- The FEM applies Z_piston as a local Robin condition at every point of the mouth disc (`packages/horn-solver/src/horn_solver/solver.py`, the `flanged_piston` branch): u_n(r) = p(r) / Z_piston. The 1D codes apply it as a lumped load on the mean mouth pressure and velocity. The two agree only when the mouth pressure is uniform; once it is not, the FEM loads each part of the mouth with the impedance of the whole piston, which is not the true nonlocal radiation load.
- The FEM drives a flat throat disc at uniform pressure (`dirichlet`) and takes Z as that pressure over the mean throat velocity. The 1D codes assume a plane wave at the throat.
- The FEM interior field has curved wavefronts and, at higher ka, higher-order modes. Both 1D models assume plane wavefronts.

Tolerance: for mouth ka at most 1.5, median 0.2 dB, p95 0.5 dB, max 1.0 dB, phase p95 3 degrees. Above ka 1.5 the comparison is reported only, because a plane-wave model is not expected to follow the 3D field there. The tolerance was written before `fem_compare.py` was first run and committed after the conical FEM run had finished (see the note at the top of this page).

| Profile | ka at most 1.5, dB median / p95 / max | phase p95 | ka above 1.5, dB median / p95 / max | Result |
| --- | --- | --- | --- | --- |
| conical | 0.60 / 0.89 / 0.94 | 3.7 | 0.36 / 0.76 / 0.78 | Fail |
| exponential | 0.40 / 1.29 / 1.53 | 11.1 | 0.54 / 1.32 / 1.58 | Fail |

The TMM (2,000 segments) and openwind agree with each other to 0.00001 dB here, so the FEM differs from both 1D codes equally.

openwind can also model curved wavefronts in a flaring bore (`spherical_waves=True`). The FEM difference by mouth ka band, magnitude median / max in dB:

| ka band | conical, plane-wave 1D | conical, spherical-wave 1D | exponential, plane-wave 1D | exponential, spherical-wave 1D |
| --- | --- | --- | --- | --- |
| 0 to 0.5 | 0.60 / 0.64 | 0.09 / 0.16 | 0.35 / 0.42 | 0.01 / 0.04 |
| 0.5 to 1 | 0.74 / 0.94 | 0.32 / 1.38 | 0.57 / 1.10 | 0.13 / 0.31 |
| 1 to 1.5 | 0.33 / 0.69 | 0.46 / 1.60 | 1.13 / 1.53 | 0.99 / 2.16 |

What is explained: below ka 0.5 the difference is a nearly constant level offset with no phase error, 0.60 dB for the cone (22.6 degree half angle) and 0.35 dB for the exponential horn. The spherical-wave 1D model removes almost all of it (0.16 dB and 0.04 dB left). That offset is the plane-wavefront approximation of the 1D model.

What is not explained: between ka 0.5 and 1.5, where most of the frozen failures sit, the FEM differs from the spherical-wave 1D curve about as much as from the plane-wave curves (up to 1.4 and 1.6 dB for the cone, 2.2 dB for the exponential horn, with phase p95 of 6 to 12 degrees). Wavefront curvature alone does not account for it. Candidate causes, none of which has been isolated:

1. the local Robin application of the piston impedance in the FEM against the lumped load of the 1D codes (above), which matters exactly when the mouth pressure becomes nonuniform at moderate ka;
2. the throat conditions: a uniform-pressure flat disc in the FEM against a plane wave in the 1D codes;
3. wavefront curvature beyond what openwind's spherical-wave option represents, and the onset of higher-order modes.

The FEM volume mesh has been ruled out (see the refinement below).

The production ranking does not use this configuration: it uses the nonlocal `modal_baffled` mouth, which removes cause 1.

What this means for users: the Webster prescreen and the `flanged_piston` FEM differ by up to about 1.5 dB in throat impedance for these 300 mm horns between ka 0.5 and 1.5, and the reason is not yet established. The prescreen is only a prescreen; final predictions use the FEM with `modal_baffled`.

Plots: [conical](validation/cross_software/fem_conical.png), [exponential](validation/cross_software/fem_exponential.png).

### FEM mesh refinement

The conical horn was rerun from 50 to 667 Hz (mouth ka up to 1.83) on the same frequency grid with `--mesh_size 0.007` and one band: `results/20261116T181626Z-bb179610`, 116,690 cells, P1, same horn-solver image, source `d1a9d38` (no solver change). In the original run the frequencies up to 500 Hz (ka 1.37) were solved on 43,927 cells and those above on 123,803 cells, so the cell count grows by a factor of 2.7 below ka 1.37; above it the two meshes are of similar size. The throat element size, 5 mm, was the same in both runs.

Criterion, committed in `fem_compare.py` before the refined run started (the run's source revision `d1a9d38` contains it): refining must change |Z| by at most 0.1 dB and 1 degree for ka at most 1.5. Result: the largest change is 0.009 dB and 0.065 degrees (0.006 dB between ka 0.5 and 1). Pass. The refined FEM differs from openwind by the same amounts as the original (median 0.60 dB, max 0.94 dB for ka at most 1.5), so the volume mesh is not the cause. Finer throat meshing and P2 elements were not tried.

## B. bempp-cl: exterior radiation

### Set-up

The exterior Helmholtz problem is solved with P1 pressure and piecewise-constant normal velocity on a closed triangulated surface, using the Burton-Miller combined equation (the plain double-layer equation failed to converge at an interior resonance of the large baffle body at 1.2 kHz). Meshes come from gmsh with outward normals taken from the tetrahedron behind each face. bempp-cl uses exp(-i w t); its results are conjugated before comparison. Air: rho = 1.225 kg/m3, c = 343 m/s. Far field at 50 m, angles 0 to 90 degrees in 5 degree steps. Each case was run at two mesh densities (element size divided by 1.4).

A rigid piston or horn mouth "in an infinite baffle" cannot be modelled directly with closed-body boundary elements. The baffle is a thick disc, radius 5 times the radiator radius and 50 mm thick, with the radiator at the centre of its front face. A second run with radius 10 times the radiator shows how much of each difference is edge diffraction from that finite baffle.

### Tolerances

| Check | Limit |
| --- | --- |
| Pulsating sphere against the closed form | 0.1 dB, 1 degree |
| Modal far field (uniform velocity) against the analytic piston pattern | 1e-6 dB |
| Modal far field against direct Rayleigh quadrature, archived horn velocities | 0.01 dB, 0.1 degree |
| Piston radiation impedance against the exact baffled piston, ka at least 2 | 0.5 dB, 5 degrees |
| Piston directivity re axis, ka at least 2, angles up to 60 degrees, where the analytic level is above -20 dB | 1.0 dB |
| Horn throat impedance against the archived modal_baffled FEM | 0.5 dB, 5 degrees |
| Horn directivity re axis against the project's modal far field, up to 60 degrees, above -20 dB | 1.0 dB |
| Horn on-axis far-field level against the modal far field | 1.5 dB |
| Mesh refinement change (impedance and directivity) | 0.1 dB, 1 degree |

### Checks of the set-up and of the directivity formula

The pulsating sphere (4,098 unknowns) agrees with the closed form within 0.015 dB and 0.2 degrees at ka 0.5, 2 and 5. Pass.

The project's experimental directivity (`horn_core.modal_radiation.modal_baffled_far_field`) evaluates the Rayleigh far field of the 16 mouth-velocity modes through their mode spectrum. Integrating the same velocity profile directly (400-point Gauss quadrature of v(r) J0(k r sin theta) r) for all 101 archived frequencies of the worked example horn gives the same result within 3e-13 dB and 2e-12 degrees. With a uniform velocity it reproduces the analytic piston pattern exactly. Pass. The quadrature uses the same velocity basis, J0(mu_n r / a) / J0(mu_n), taken from `modal_boundary.aperture_projection`. It therefore confirms the Hankel-transform evaluation of the far field for that basis. That the solver's modal coefficients are expressed in the same basis was checked by reading the code, not by this test.

### Piston in a finite baffle

Piston radius 100 mm.

| ka | Unknowns | Z/rho c, BEM | Z/rho c, exact baffled piston | Z error dB / degrees | Directivity error, dB | On-axis, finite minus infinite baffle, dB | Result |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 1,057 | 0.139 + 0.386i | 0.120 + 0.397i | 0.08 / 3.0 | 5.13 | +4.17 | below ka 2: report only |
| 1 | 1,173 | 0.441 + 0.673i | 0.423 + 0.647i | 0.35 / 0.05 | 2.00 | -0.65 | below ka 2: report only |
| 2 | 1,981 | 1.013 + 0.520i | 1.033 + 0.535i | 0.18 / 0.20 | 3.10 | +1.97 | Fail (directivity) |
| 4 | 4,051 | 0.942 + 0.115i | 0.941 + 0.122i | 0.002 / 0.41 | 0.20 | -0.04 | Pass |

With the 10x baffle at ka 2 the directivity error falls from 3.10 to 1.16 dB, the on-axis difference from +1.97 to -0.15 dB, and the impedance error from 0.18 to 0.05 dB. At ka 1 the 10x baffle still gives +3.0 dB on axis. The remaining differences move with baffle size, so they are edge diffraction from the finite baffle, not an error in the piston formulas. The exact baffled-piston impedance used by the Webster prescreen and `flanged_piston` is confirmed within 0.2 dB from ka 2 upward, and the radiation resistance and reactance follow it at lower ka with the expected finite-baffle offset.

### Worked example horn in a finite baffle

The archived hyperbolic horn of the 800 to 1600 Hz worked example (throat radius 32.5 mm, mouth radius 88.7 mm, length 117.9 mm) is cut into the front of the baffle body and driven by a uniform throat velocity. This is a full 3D model of the horn interior and the exterior together, independent of the project's FEM and modal aperture. It is compared at the archived FEM frequencies nearest 800, 1000, 1200, 1400 and 1600 Hz with `outputs/auto/refinement/refine_hyperbolic_0001_results.csv` from `data/validation/worked_example_modal_800_1600.tar.gz` (unit throat pressure, `modal_baffled`, infinite baffle).

| f (Hz) | mouth ka | Unknowns | Z/rho c BEM | Z/rho c FEM modal | Z error dB / degrees | TMM piston vs BEM, dB | Directivity error, dB | On-axis BEM minus FEM, dB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 801 | 1.30 | 1,826 | 2.119 + 0.846i | 2.121 + 0.746i | 0.13 / 2.4 | 0.96 | 7.18 | -5.41 |
| 998 | 1.62 | 2,186 | 1.253 + 0.024i | 1.322 + 0.020i | 0.47 / 0.25 | 2.60 | 1.84 | +1.01 |
| 1199 | 1.95 | 2,534 | 1.143 + 0.169i | 1.145 + 0.201i | 0.05 / 1.6 | 2.66 | 4.04 | +2.48 |
| 1389 | 2.26 | 2,977 | 1.255 + 0.243i | 1.273 + 0.227i | 0.10 / 0.86 | 1.67 | 0.77 | -0.04 |
| 1598 | 2.60 | 3,420 | 1.292 + 0.004i | 1.338 + 0.021i | 0.30 / 0.72 | 0.44 | 3.61 | -2.59 |

Throat impedance passes at all five frequencies (largest 0.47 dB, 2.4 degrees). The two throat conditions are not identical: the BEM drives the throat with a uniform normal velocity, the archived FEM with a uniform pressure (`dirichlet`), and each takes Z as mean pressure over mean velocity. They differ slightly when the throat field is not plane; throat ka reaches 0.95 at 1.6 kHz. This is an independent 3D confirmation of the production `modal_baffled` throat loading for this horn, to within the finite-baffle effect. The 1D Webster prescreen with the piston load differs from the BEM by up to 2.7 dB here. That is a different quantity from, and larger than, the 1.05 dB throat-impedance error of the FEM with a local `flanged_piston` load recorded for another horn in [EXTERIOR_RADIATION_VALIDATION.md](EXTERIOR_RADIATION_VALIDATION.md); it is reported here, not explained.

Directivity fails at four of five frequencies and on-axis level at three. These cannot be separated from the finite baffle: the mouth ka is only 1.3 to 2.6, where the piston case also shows several dB of baffle diffraction. With the 10x baffle the on-axis difference at 801 Hz changes from -5.41 to -4.28 dB and at 1199 Hz from +2.48 to -2.41 dB, so it oscillates with baffle size as edge diffraction does. The project's directivity output is therefore correct for an infinite baffle (Rayleigh check above), but a real cabinet or a baffle of a few mouth diameters will change it by several dB in this ka range, as its label already says.

### Mesh refinement

Seven of nine refinement pairs meet 0.1 dB / 1 degree. The impedance changes by at most 0.06 dB and 0.46 degrees in every pair. Two directivity changes exceed 0.1 dB (0.10 dB for the piston at ka 4 and 0.16 dB for the horn at 1 kHz), both far smaller than the differences being reported. Fail as frozen, but immaterial to the conclusions.

### Full 3D horn BEM at larger size

The 300 mm horns of part A were not modelled with bempp-cl. Their 150 mm mouth reaches ka 13.7 at 5 kHz, and a dense BEM at six elements per wavelength on a baffle body several mouth diameters across needs tens of thousands of unknowns, more than fits in memory on this machine without a fast multipole backend. The worked example horn, at 1.6 kHz and below, needed at most 6,614 unknowns.

Plot: [directivity](validation/cross_software/bempp_directivity.png).

## C. Hornresp

Hornresp is a Windows program and Wine is not installed here, so no Hornresp result exists yet. [validation/hornresp](../validation/hornresp/README.md) contains:

- input sheets for the seven horns above and for the 800 to 1600 Hz worked example horn with the 18Sound 6NMB420 driver, with values in Hornresp units;
- our TMM throat impedance for each geometry exactly as entered;
- the same values as Hornresp import records (`records/import_unverified/`), written from the record layout published in user forums and not checked against the program;
- `compare_hornresp.py`, which reads a Hornresp text export and applies tolerances frozen before any Hornresp data: impedance median 0.2 dB, p95 1.0 dB, max 3.0 dB, phase p95 5 degrees, peaks within 1 % and the same peak count; SPL from 800 to 1600 Hz within median 0.5 dB and max 2.0 dB against two references.

The two SPL references differ on purpose. The first is the archived modal FEM on-axis response. Hornresp is a plane-wave model, and our own 1D model built the way Hornresp builds SPL already differs from that FEM curve by median 1.06 dB and max 1.41 dB, so this check is expected to fail even if Hornresp agrees exactly with our 1D model. The second reference, added in the review pass before any Hornresp data, is that like-for-like 1D curve: TMM horn, ideal transformer, mouth volume velocity into half space. Piston directivity does not affect either on-axis comparison: in the far field the on-axis level of a baffled piston equals that of a half-space point source with the same volume velocity, and the near-field correction at 1 m is 0.02 dB.

Conical, exponential, hyperbolic (Hornresp Hyp with T = 0), tractrix and CD (Exp then Con) are entered exactly. OS and the scaled tractrix have no Hornresp flare; they are entered as four Con or Exp segments with 1.30 % and 1.36 % largest radius error, and our reference curve uses the same segmented shape. The export parser is tested on synthetic files (`scripts/tests/test_hornresp_compare.py`); it has not seen a real Hornresp export, so column names may need the override options.

## Limits

- All comparisons are numerical. None is a measurement, and none qualifies a physical driver, adapter or enclosure.
- openwind and the TMM are both 1D plane-wave models. Their agreement checks our implementation of that model, not the model.
- The bempp-cl baffle is finite. An infinite baffle cannot be modelled with these closed-body boundary elements, so the directivity of the project's infinite-baffle output could only be bounded, not confirmed, below ka of about 4.
- Only two FEM runs were compared with the 1D codes, both with the `flanged_piston` load, plus one mesh refinement run of the cone. The cause of their difference between ka 0.5 and 1.5 is not established. The production `modal_baffled` FEM is compared here only through the archived worked example.
