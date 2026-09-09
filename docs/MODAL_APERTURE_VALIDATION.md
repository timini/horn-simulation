# Independent modal-aperture qualification

The new nonlocal mouth boundary **passes** the frozen independent comparison for the worked hyperbolic horn over 800–1600 Hz. This is numerical qualification of an ideal infinite-baffle problem, not measured qualification of a physical horn/driver assembly.

The earlier local Robin boundary comparison remains a failure in [the exterior-radiation report](EXTERIOR_RADIATION_VALIDATION.md). No tolerance, source amplitude, geometry or phase offset was fitted to obtain the new result.

## Fixed problem and reference

The throat radius is 32.498 mm, mouth radius 88.70898640584517 mm and length 117.90625 mm. The hyperbolic profile is the same 80-section STEP used by the earlier production comparison. Excitation is uniform throat velocity, 1 m/s RMS; air is rho=1.225 kg/m³ and c=343 m/s. The observer is on-axis, 1 m from the mouth. Every case uses the same 33 logarithmically spaced frequencies and `exp(+iwt)` convention. Integrated FEM inlet flow is normalized to the nominal circular throat area.

The reference is the separately executed, pinned MMM toolbox with 64 modes and 1,000 axial sections. Its own mode, section and radial-observer refinements passed before this boundary implementation was introduced. Its archive hashes are checked before preparation. The new boundary uses an independently implemented Fourier–Bessel radiation integral coupled to the three-dimensional P1 FEM; it does not call the reference toolbox.

## Executed result

| Comparison | Maximum impedance magnitude change | Impedance phase | Maximum pressure magnitude change | Pressure phase |
| --- | ---: | ---: | ---: | ---: |
| 8 → 16 aperture modes, 4 mm mesh | 0.001103 dB | 0.00813° | 0.000969 dB | 0.00461° |
| 4 → 3 mm mesh, 16 modes | 0.011350 dB | 0.06142° | 0.005775 dB | 0.11253° |
| Final FEM → independent MMM | 0.016993 dB | 0.15873° | 0.012094 dB | 0.18339° |

All three comparisons pass. Frozen refinement limits are 0.05 dB and 1°; independent comparison limits are 0.5 dB and 5°. The three sweeps contain 99 frequency solves. Meshes contain 74,329 and 172,448 tetrahedra. Maximum relative linear residual is 1.74e-12, interface residual 1.60e-14 and input/radiation power imbalance 5.78e-14. Every computed radiation power is positive.

## Evidence and reproduction

- `data/validation/modal_mouth_manifest.json` identifies the exact executed source revision and archive digests.
- `modal_mouth_artifacts.tar.gz` contains protocols, raw curves, reference inputs and provenance, the solve log, original solver image ID, comparison, output seal and an archive of the executed package/harness source.
- `modal_mouth_reference.json` is the machine-readable result. It retains `physical_validation_status: experimental_prediction`.

At the manifest's source revision, with a clean checkout:

```sh
python scripts/validate_modal_mouth.py prepare results/qualification/modal-mouth
# In a solver container mounted at /workspace, with PYTHONPATH=/usr/local/lib:
python3 scripts/validate_modal_mouth.py solve results/qualification/modal-mouth
# Back on the host, using the same source:
python scripts/validate_modal_mouth.py compare results/qualification/modal-mouth
```

The reference input is taken from the checked-in MMM archive; the toolbox does not need to be rerun for this comparison. To reproduce that independent reference itself, follow the exterior-radiation report. Solves and comparisons refuse changed source, changed inputs, missing frequency rows or changed output seals.

This result covers one profile, geometry and target band under uniform inlet velocity. Production pressure-transfer integration is now exercised by the [complete automatic workflow](MODAL_APERTURE.md), its [candidate-resolution study](MODAL_CANDIDATE_RESOLUTION.md), pressure/velocity solver tests and the [four-band search audit](MODAL_SEARCH_VALIDATION.md). Additional profiles/bands, finite-baffle behavior and actual driver/interface data still need their own evidence. Reports continue to show insufficient evidence for an uncharacterized assembly.
