# Duct losses and finite-flange validation

Latest expanded run: [8 September reference validation report](REFERENCE_VALIDATION_REPORT.md).

The earlier lossless brass-pipe diagnostic failed: first impedance peaks were 25–32 dB too high. It is preserved in `data/validation/pipe_diagnostic_baseline.json`. The implementation now includes a published thin-boundary-layer wall condition in the production three-dimensional FEM solver, matched complex propagation in TMM, and finite-flange pipe radiation. No measured curve was fitted.

## Model and supported domain

All phasors use RMS amplitudes and `exp(+iωt)`. The wall terms follow [Berggren, Bernland and Noreland, equation 31](https://arxiv.org/abs/1801.04177). A tangential-gradient term represents viscous dissipation; a pressure term represents thermal dissipation. The TMM uses the matching thin-layer density/compressibility; the core also implements the exact circular Kirchhoff reference. RMS input power equals radiation plus wall dissipation in the production reference test.

Boundary-layer mode is opt-in: `--loss_model boundary_layer --minimum_wall_scale <metres>`. The caller must declare a conservative minimum wall-curvature/gap scale. Both viscous and thermal layers must be at most 10% of that scale; narrow gaps beyond this approximation are rejected. This declaration is not an automatic CAD curvature analysis. Quadratic elements (`--element_degree 2`) also use quadratic CAD mesh geometry, which is necessary to avoid circular-wall faceting error. The lossless default remains explicitly experimental.

Finite-flange radiation uses [Silva et al.](https://arxiv.org/abs/0811.3625) and the Dalmont finite-flange relation used in the [benchmark's analysis](https://gitlab.inria.fr/aernoult/acoustic-impedance-benchmark). The exposed domain is circular outlets, `ka < 1.5`, and flange width/radius at most 1. It is a pipe termination model, not a solution for arbitrary horn rollback, baffle diffraction or exterior geometry. Zero flange width reduces to the unflanged limit. `closed` means a rigid outlet; loss terms are applied to the side wall, not a thermal boundary layer on an end cap.

## Independent data and frozen checks

Source: [Ernoult et al., Zenodo v2, record 20024938](https://zenodo.org/records/20024938), CC BY 4.0. The archived download is pinned by SHA-256 in `data/validation/references.json`. The importer retains source filenames, units, experimental operators and repeated measurements. The 299 measurement files are **not 299 different assemblies**.

All eight dense predictions were written and hashed before the comparison script opened measurement values. Brass open cylinders were development cases; the other 263 curves were held out. Models use the source's dry-air 25°C constants and nominal dimensions. No level shift, frequency shift, damping fit or threshold adjustment was applied. Per-curve magnitude gates are median absolute error ≤2 dB and 95th percentile ≤4 dB over 110–3900 Hz. Complex-error statistics are recorded separately; a magnitude pass does not establish a phase pass.

| Case | Curves | Magnitude passes | Median of per-curve 95th-percentile errors |
| --- | ---: | ---: | ---: |
| Brass, open, 2 mm flange | 36 | 36 | 1.14 dB |
| Wood, open, 7 mm flange | 36 | 29 | 3.41 dB |
| Printed cylinder, open | 37 | 36 | 1.63 dB |
| Brass, closed | 45 | 41 | 2.45 dB |
| Wood, closed | 45 | 21 | 4.08 dB |
| Printed cylinder, closed | 46 | 33 | 2.22 dB |
| Cone, open | 28 | 26 | 2.94 dB |
| Cone, closed | 26 | 15 | 3.52 dB |

Overall: **237/299 pass**, including **201/263 held-out curves**. These failures remain failures. Wall material, sealing, nominal versus actual dimensions and end-cap effects are plausible contributors; they have not been established as causes. General physical validation is therefore not claimed.

Four independently supplied one-dimensional FEM reference curves all pass, with 95th-percentile magnitude errors below 0.010 dB and maximum below 0.082 dB. Direct production FEM comparisons for brass-open, wood-open and cone-open cases pass 90 of 100 measured curves on a 121-frequency grid. The three-mesh production wall-loss test passes complex-impedance, ≤0.5 dB finest-mesh output change, passivity and energy-balance checks.

The machine-readable summary is `data/validation/loss_physics_validation.json`. Local complete comparisons, source/prediction hashes and a plot are under `results/validation-81/loss-physics/held-out-final/`.

## Reproduction and limits

After importing the pinned research archive with `horn_analysis.reference_data`, run inside the solver image with repository sources mounted:

```sh
python3 scripts/validate_duct_physics.py \
  --reference-dir results/validation-81/imported/ernoult-pipe-impedance-v2 \
  --output-dir results/validation-81/loss-physics/reproduction --with-fem
```

The production reference test is `packages/horn-solver/tests/validation/test_wall_losses.py`; CI requires it without skips. The full solver suite passed 48 cases with one unrelated optional straight-tube regime skipped; the required release lane permits no skips.

This evidence tests duct input impedance. It does **not** validate a commercial driver's chamber/phase plug, absolute listening-distance SPL, off-axis response, or recommendation ordering between physical horn assemblies. Those gates remain open in issue 81. The legacy FEM–BEM coupling is disabled at both public entry points because its whole-boundary trace did not implement a mouth-only exterior problem; standalone BEM operator tests are retained.
