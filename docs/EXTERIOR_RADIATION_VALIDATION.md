# Complete-horn radiation validation

A correct interior FEM solve and a correct uniform-piston formula do not establish the radiation of a complete horn. The mouth velocity can vary across the aperture, and a real baffle/body can change both loading and the listening-position response. This work separates those effects for the hyperbolic horn in the [800–1600 Hz worked example](WORKED_EXAMPLE.md).

## Matched infinite-baffle comparison

The independently maintained [MMM Toolbox](https://github.com/bkolbrek/MMM_toolbox) implements mode matching for axisymmetric horns, including modal radiation into an infinite baffle. The reference is pinned to `24d533b191f299aafe1c6aae28cba294abc8f8e2` and executed as a separate Octave process. Its acoustic code is not included in this repository.

The frozen study uses the worked candidate's mathematical hyperbolic profile, 1 m/s RMS uniform throat velocity, 1.225 kg/m³ air density, 343 m/s sound speed and an observation point 1 m on axis from the mouth. Both sides use `exp(+iwt)`. The reference computes the full modal mouth loading and a spatial Rayleigh integral; production uses its local flanged-piston boundary and uniform-mouth observer approximation. Polygonal inlet areas are corrected to the same nominal volume velocity before comparing output. No fitted level shift is allowed.

Thirty-three logarithmically spaced frequencies cover 800–1600 Hz. Reference refinements use 16, 32 and 64 radial modes; 250, 500 and 1,000 axial sections; and 201 versus 401 radial integration points. Production uses a 4 mm P1 mesh and an 80-section loft, whose resolution is separately studied in [the candidate report](CANDIDATE_RESOLUTION.md). The reference uses directly integrated modal radiation matrices, without the toolbox's precomputed interpolation tables or high-frequency approximation.

Before the dense comparison, the protocol fixes final-reference refinement limits of 0.05 dB and 1 degree for complex pressure and throat impedance, and model-comparison targets of 0.5 dB in pressure/impedance magnitude and 5 degrees in phase. These strict engineering targets assess this approximation; they are not measurement uncertainty or a universal standard. Preliminary three-frequency studies informed the method, not a physical calibration. Any failed target remains a failure and must be investigated.

Reproduction requires a clean source checkout and a pristine separate clone of the pinned toolbox:

```sh
docker build -f scripts/mmm_reference.Dockerfile -t horn-octave-reference:local .
python scripts/validate_mmm_radiation.py prepare results/mmm-study
python scripts/validate_mmm_radiation.py reference results/mmm-study \
  --checkout /path/to/pinned/MMM_toolbox

docker run --rm -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
  -v "$PWD:/workspace" -w /workspace -e PYTHONPATH=/usr/local/lib \
  horn-solver:latest python3 scripts/validate_mmm_radiation.py production results/mmm-study

python scripts/validate_mmm_radiation.py compare results/mmm-study
```

The host needs NumPy, SciPy and pandas. The adapter generates and preserves Bessel roots as reference inputs, provides only a scalar `contains` compatibility helper for Octave 7, and verifies the upstream entry functions' locations. Reference source is mounted read-only. The protocol records source/input hashes; solve evidence seals all output files. The reference image ID and runtime version are retained. A Git worktree should be prepared/compared on its host; only the production stage runs in the solver container. Preserve failed runs and use a fresh directory after any source change.

### Executed result

All five reference sweeps (165 frequency solves) and the 33-point production sweep completed. The reference passes its declared refinement checks: 32→64 modes changes pressure by at most 0.0257 dB and impedance by 0.0399 dB; 500→1,000 axial sections changes them by 0.0073 and 0.0112 dB. Increasing radial integration from 201 to 401 points changes pressure by less than 0.000006 dB. These are observed refinement changes, not an asymptotic error bound.

**The strict model comparison fails.** Production versus the final reference differs by at most 0.39194 dB / 3.1725 degrees in on-axis pressure, and 1.04964 dB / 5.7307 degrees in throat impedance. Output pressure meets its frozen targets; throat impedance exceeds both its 0.5 dB and 5-degree limits. The approximation must not be promoted to a high-accuracy horn-loading model on the strength of its closer pressure result. The discrepancy can affect driver coupling and recommendation ordering.

The [comparison JSON](../data/validation/mmm_radiation_reference.json), [raw archive](../data/validation/mmm_radiation_artifacts.tar.gz) and [identity manifest](../data/validation/mmm_radiation_manifest.json) preserve the failure, exact reproduction source, pristine upstream identity, numerical inputs, response arrays, console logs and runtime/image identities. Reproduce the historical study from the manifest's source commit; do not overwrite its inputs to apply newer source. The earlier v1 attempt was stopped after the production stage rejected harmless frequency-grid roundoff between host and container NumPy versions; v2 fixes that check and reruns both sides. Neither result is a physical assembly measurement.

## Finite-baffle diagnostic pilot

The separate Boundary Lab pilot uses the exact worked-example air-volume STEP inside a rigid cylindrical body. Its front face is flush with the mouth, the back is 10 mm behind the throat, and the throat floor supplies uniform velocity into the horn. All other surfaces are rigid. The outward-oriented surface is closed and checked for watertight edges and source orientation. This is a defined finite mounting geometry, not an infinite-baffle reference or a proposed fabrication drawing.

Pinned Boundary Lab `8cb166226e412877d3f71f2845918e479b97aa85` solved the complete exterior Burton–Miller Neumann problem. The source derivative is checked against unit velocity, and integrated source pressure against the reported mechanical force/velocity impedance. Conversion to our specific throat impedance uses both the source area and the opposite phasor convention. These closure checks do not replace a linear-system residual; this backend does not export an exterior residual.

The table gives **reference minus production** on-axis output differences, with no fitted gain:

| Baffle/body radius | Interior / outer surface mesh target | 800 Hz | 1200 Hz | 1600 Hz |
| --- | --- | --- | --- | --- |
| 250 mm | 8 / 40 mm, 1,814 triangles | +3.022 dB | −3.675 dB | +0.362 dB |
| 250 mm | 4 / 20 mm, 6,694 triangles | +2.972 dB | −3.835 dB | +0.472 dB |
| 500 mm | 8 / 40 mm | +2.799 dB | −2.080 dB | −0.046 dB |

These rows use double precision and fixed fourth-order regular quadrature. Increasing that quadrature to sixth order on the coarse 250 mm case changes output by less than 0.006 dB. The initial single-precision, automatically selected second-order calculation differs from the double-precision fourth-order result by less than 0.007 dB. Refining the 250 mm surface changes output by up to 0.161 dB at the three points. Changing baffle radius produces a much larger response/phase change, so the finite mounting cannot be ignored.

This pilot **does not pass an accuracy or physical-validation gate**. Three frequencies do not bound continuous-band error, the larger baffle is not separately mesh-qualified, and neither finite body is the same physical boundary as an infinite baffle. Do not attribute the entire discrepancy to the interior FEM or apply it as a fitted correction to other horns.

The [comparison and raw-file hashes](../data/validation/exterior_baffle_pilot_reference.json), [raw pilot archive](../data/validation/exterior_baffle_pilot_artifacts.tar.gz) and [identity manifest](../data/validation/exterior_baffle_pilot_manifest.json) retain the 15 exterior solves, source/geometry inputs and authored pilot scripts. Those historical scripts contain the original local paths; the archived inputs permit re-execution with equivalent path mappings. Run `scripts/audit_exterior_pilot.py` against the extracted root to repeat the source/units checks and numerical comparison. No physical measurements are included in this archive.

## Release implication

The default prediction assumes infinite-baffle radiation and a simplified mouth field. A free-standing horn or a finite mounting needs its own radiation evidence. A driver also needs its characterized interface and rear load: an ideal source comparison cannot supply them. [Issue #81](https://github.com/timini/horn-simulation/issues/81) remains open for the complete physical design gate, and the [physical-reference audit](PHYSICAL_REFERENCE_AUDIT.md) records the usable measurements and missing information for each current lead.
