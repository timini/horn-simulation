# Independent numerical solver validation

This suite supports [#78](https://github.com/timini/horn-simulation/issues/78). It tests the implemented equations with matched inputs; it does not qualify a commercial driver, its chamber/adapter, an exterior radiation model or a physical assembly.

## Frozen comparison

[Protocol](../data/validation/boundary_lab_protocol.json): a straight tube and a cone, each on three shared P1 tetrahedral meshes, with 13 logarithmic frequencies from 250 to 2000 Hz. Both solvers read the **same Gmsh 4.1 mesh**, including inlet/mouth/wall groups. The walls are rigid, the mouth has a plane-wave termination, and air properties match exactly.

The independent implementation is Boundary Lab commit `8cb166226e412877d3f71f2845918e479b97aa85`, Python 3.11 and Julia 1.12.6. Its native 2.83 V RMS electrodynamic piston has a synthetic dry diaphragm mass, known compliance/damping/inductance and zero rear load. Its effective piston area is the actual polygonal mesh inlet area, not an assumed CAD circle. The reference's `exp(-i omega t)` quantities are conjugated into this project's `exp(+i omega t)` convention before comparison.

The production solver uses prescribed inward velocity. The comparison checks complex throat impedance, mouth pressure per inlet velocity, motor current and velocity, throat pressure and coupled mouth pressure. It applies the existing production driver operating-point calculation to the production FEM load, then compares against Boundary Lab's independently coupled solution. The frozen maximum relative complex error is `1e-4`; both linear solvers must meet a `1e-8` residual limit. Every expected frequency must be present and healthy. There is no fitted phase or level shift.

Agreement on shared meshes separates implementation errors from discretization differences. It does **not** establish continuum convergence: the polygonal inlet area itself changes with mesh resolution. Nor does ideal zero-rear-load coupling validate the scraped manufacturer's Mms, phase plug or rear enclosure.

## Recorded result

All **6/6 cases, 13 frequencies each**, pass. The largest complex relative discrepancy across all six compared quantities is **8.68e-13**, below the frozen 1e-4 limit; the largest independent residual is **1.74e-14**, below 1e-8. This near-machine-precision agreement is expected when two correct implementations solve the same discrete equations. It is not a physical accuracy estimate.

| Shared mesh case | Maximum complex relative discrepancy |
| --- | ---: |
| Tube, 8 mm | 2.57e-13 |
| Tube, 6 mm | 1.49e-13 |
| Tube, 4 mm | 3.26e-13 |
| Cone, 8 mm | 1.95e-13 |
| Cone, 6 mm | 1.03e-13 |
| Cone, 4 mm | 8.68e-13 |

The [numeric reference](../data/validation/boundary_lab_reference.json), [complete 2 MB artifact archive](../data/validation/boundary_lab_reference_artifacts.tar.gz) and [identity manifest](../data/validation/boundary_lab_reference_manifest.json) retain original meshes, STEP geometry, both solver outputs, logs, input/source hashes and independent runtime details. The manifest records the exact source commit used for reproduction. Six offline regression tests use the independent complex loads and motor results in ordinary analysis CI.

The earlier Gmsh 2.2 pilot failed the independent reader; a later attempt timed out under shared host load and exposed worker cleanup. Neither counts as a pass. The accepted run starts fresh, uses 4.1 input and the corrected process-group cleanup, and completes every case. Local failed-attempt logs are retained; the checked-in archive contains the accepted complete run.

## Reproduction

Use a clean pinned Boundary Lab checkout with its Python dependencies and Julia project installed according to that checkout's instructions. The runner verifies the imported module comes from that checkout, checks its revision/cleanliness and runtime versions, and records installed Python packages. Do not reuse a mutable checkout being edited by another task.

The preparation/comparison Python environment needs NumPy, pandas, SciPy, Gmsh, meshio and this repository's `horn-core`, `horn-drivers` and `horn-analysis` packages. The production stage uses the repository's solver Docker image. Set `JULIA_DEPOT_PATH` if the independent installation has a dedicated depot.

```sh
python scripts/validate_boundary_lab.py prepare results/qualification/new-comparison

docker run --rm -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
  -v "$PWD:/workspace" -w /workspace \
  -e PYTHONPATH=/usr/local/lib:/workspace/packages/horn-core/src:/workspace/packages/horn-drivers/src:/workspace/packages/horn-analysis/src:/workspace/packages/horn-solver/src:/workspace/packages/horn-geometry/src \
  horn-solver:latest python3 scripts/validate_boundary_lab.py horn results/qualification/new-comparison

python scripts/validate_boundary_lab.py reference results/qualification/new-comparison \
  --checkout /path/to/pinned/boundary-lab \
  --python /path/to/boundary-lab-runtime/bin/python \
  --julia /path/to/julia-1.12.6/bin/julia

python scripts/validate_boundary_lab.py compare results/qualification/new-comparison
```

Preparation refuses an existing directory. It freezes geometry, input and source hashes. Each solver stage seals its output hashes. Comparison verifies these hashes, units, shapes, mesh identity, frequency completeness, conventions and residuals before publishing `comparison.json`; changed evidence fails. Keep the whole run directory, not just its summary. A timeout or partial solve has no passing final evidence. Start a fresh directory after a source change or failed attempt. The external runner bounds each case to ten minutes and kills its own remaining process group, including workers surviving the CLI parent.

## Required production checks

The velocity-inlet test solves an analytically matched travelling-wave tube with complex RMS velocity. It checks complex inlet/mouth pressure, integrated flow, absolute acoustic power and power conservation; invalid prescribed velocities are rejected.

The cavity check uses the **production volume operator** to assemble a rigid rectangular cavity, then independently compares its first twelve nonzero eigenfrequencies with analytical box modes. All twelve errors decrease through three P2 mesh refinements; the finest maximum error is approximately **0.033%**, below the frozen 1% bound. The constant-pressure zero mode is also checked.

These six additional cases join the existing seventeen in the required acoustic CI lane. That lane rejects skips as well as failures. The standalone Python CI lane runs the comparison-input/evidence safeguards and subprocess cleanup tests without requiring Julia.

Physical-reference suitability is recorded separately in [the physical audit](PHYSICAL_REFERENCE_AUDIT.md). No numerical pass changes a candidate's physical validation status.

## Independent radiation-integral check

Ten pure-Python cases integrate the [Rayleigh monopole kernel](https://euphonics.org/4-3-2-the-rayleigh-integral-and-the-baffled-piston/) independently of the production closed forms. They cover three radii, near-to-far axial distances and `ka` from 0.01 to 20, including complex RMS amplitude and phase. A double surface integral, reduced using disk-overlap area, checks surface-averaged radiation load without Bessel or Struve functions. Two quadrature resolutions must agree before comparison with production; all ten cases pass the fixed tolerances.

This verifies the implemented uniform circular piston in an infinite baffle. It is not a COMSOL run, a full exterior horn solve, finite-baffle validation or evidence that an arbitrary horn mouth has uniform velocity. The [COMSOL verification model](https://doc.comsol.com/6.3/doc/com.comsol.help.models.aco.baffled_piston_radiation/baffled_piston_radiation.html) describes this distinction between the surface integral and the general exterior-field calculation.
