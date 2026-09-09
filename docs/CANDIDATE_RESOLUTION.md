# Numerical resolution of a generated candidate

A completed optimisation is a search result, not evidence that its mesh or frequency grid is adequate. `scripts/validate_candidate_resolution.py` freezes a selected candidate, its driver record and its target band, then checks the production geometry/solver/coupling at increasing resolutions. It never changes the physical evidence status or certifies a driver interface.

The initial protocol supports the default lossless, flanged-piston, P1 model. It uses the candidate's exact profile, dimensions, voltage and observation distance. Its simulation range extends half an octave beyond each target edge. Three mesh sizes (10, 6 and 4 mm), three loft counts (20, 40 and 80 sections), and two nested grids (101 and 201 frequencies) isolate mesh, geometry and sampling changes. The solver may impose a smaller wavelength-based mesh limit; each output records the actual cell count.

Before solving, the study freezes these engineering gates:

| Change | Maximum coupled output change | Maximum throat impedance magnitude change | Maximum throat impedance phase change |
| --- | --- | --- | --- |
| Each successive mesh pair | 0.5 dB | 0.5 dB | 5 degrees |
| Each successive loft pair | 0.5 dB | 0.5 dB | 5 degrees |
| Frequency grid | 0.2 dB | 0.2 dB | 2 degrees |

Frequency refinement also limits the change in target-band ripple to 0.2 dB. Comparisons use the entire finer grid, so an additional peak or notch cannot be discarded by downsampling. Target-band ripple includes both band edges. Each solve must have a complete frequency grid, finite values, the specified model/phase contract, residual at most 1e-8, positive convergence reason, passive input load and relative acoustic power imbalance at most 1e-7.

These are resolution-change bounds for one fixed model, not an asymptotic error estimate, uncertainty in a physical loudspeaker, universal acoustic standards or proof that no still-narrower feature exists. They do not establish the ranking of nearly tied candidates. Preserve any failing study and investigate it before freezing a new protocol.

## Run a study

Use the normal prerequisites and build the images. In a clean checkout, select an actual ranking entry (index zero is the leading result) and its single driver JSON record:

```sh
python scripts/validate_candidate_resolution.py prepare results/my-resolution-study \
  --ranking results/my-design/outputs/auto/report/auto_ranking.json \
  --driver data/drivers-curated/18Sound/18sound-6nmb420.json \
  --f-low 800 --f-high 1600 --candidate-index 0

docker run --rm -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
  -v "$PWD:/workspace" -w /workspace -e PYTHONPATH=/usr/local/lib \
  horn-solver:latest python3 scripts/validate_candidate_resolution.py solve results/my-resolution-study

python scripts/validate_candidate_resolution.py compare results/my-resolution-study
```

The host needs NumPy and the comparison stage additionally needs pandas and the dependencies of `horn-core`, `horn-drivers` and `horn-analysis`. The script loads production Python modules from the checkout it hashes. Retain the solver image identity and full console log alongside the study. Each directory is new; there is no partial-resume shortcut. Source changes or edited inputs invalidate comparison. A completed solve seals every STEP and response file; the comparison rejects missing or modified evidence and writes explicit passes/failures to `comparison.json`.

A pass supports numerical exploration of that candidate within the stated model and band. It does not turn an acoustic air-volume STEP into a mounting drawing, establish rear-load or phase-plug behaviour, or replace the two complete assembly references required by [the roadmap](SINGLE_HORN_ROADMAP.md).

## Executed 800–1600 Hz example

All six solves and all five comparisons passed for the hyperbolic/18Sound 6NMB420 candidate in the [worked example](WORKED_EXAMPLE.md). The full comparison band was approximately 566–2263 Hz.

| Refinement | Maximum output change (dB) | Maximum throat impedance change (dB) | Maximum phase change (degrees) |
| --- | --- | --- | --- |
| Mesh 10 → 6 mm | 0.03964 | 0.11677 | 0.52536 |
| Mesh 6 → 4 mm | 0.01298 | 0.03757 | 0.16794 |
| Loft 20 → 40 sections | 0.0006692 | 0.001197 | 0.005445 |
| Loft 40 → 80 sections | 0.00001859 | 0.00001896 | 0.0007732 |
| Frequency 101 → 201 points | 0.002646 | 0.008821 | 0.04889 |

Target-band ripple changed by 0.001501 dB when doubling frequency sampling. Meshes contained 5,333 to 74,334 tetrahedra. All residuals were below 1e-15 and relative power imbalance below 1e-13. These changes establish stability of the evaluated model at these resolutions; they do not bound physical model error.

The [comparison JSON](../data/validation/candidate_resolution_reference.json), [complete raw archive](../data/validation/candidate_resolution_artifacts.tar.gz) and [identity manifest](../data/validation/candidate_resolution_manifest.json) retain every geometry, response, frozen input/source identity, solver image ID and console log. The manifest records the exact reproduction commit; use that commit to reproduce this historical study, or prepare a new study for newer source.
