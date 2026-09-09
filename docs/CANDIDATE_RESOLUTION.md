# Numerical resolution of a generated candidate

A completed optimisation is a search result, not evidence that its mesh or frequency grid is adequate. `scripts/validate_candidate_resolution.py` freezes a selected candidate, its driver record and its target band, then checks the production geometry/solver/coupling at increasing resolutions. It never changes the physical evidence status or certifies a driver interface.

The initial protocol supports the default lossless, flanged-piston, P1 model. It uses the candidate's exact profile, dimensions, voltage and observation distance. Its simulation range extends half an octave beyond each target edge. Three mesh sizes (10, 6 and 4 mm), three loft counts (20, 40 and 80 sections), and two nested grids (101 and 201 frequencies) isolate mesh, geometry and sampling changes. Every spatial refinement runs at 201 frequencies; the finest geometry also runs at 101 points for the sampling comparison. The actual response that produced the ranking is compared with the first refined solve, and the three mesh cases use the original ranked STEP. Preparation rejects bands whose wavelength-based cap would collapse any mesh refinement. Comparison also requires strictly increasing cell counts across the three mesh cases.

Before solving, the study freezes these engineering gates:

| Change | Maximum coupled output change | Maximum throat impedance magnitude change | Maximum throat impedance phase change |
| --- | --- | --- | --- |
| Each successive mesh pair | 0.5 dB | 0.5 dB | 5 degrees |
| Each successive loft pair | 0.5 dB | 0.5 dB | 5 degrees |
| Frequency grid, including original ranking → baseline | 0.2 dB | 0.2 dB | 2 degrees |

Frequency refinement also limits the change in target-band ripple to 0.2 dB. Comparisons interpolate SPL and complex impedance components against logarithmic frequency and use the entire finer grid, so an additional peak or notch cannot be discarded by downsampling. Target-band ripple includes both band edges with the same logarithmic-frequency interpolation used by the production ranking. Each solve must have a complete frequency grid, finite values, the specified model/phase contract, residual at most 1e-8, positive convergence reason, passive input load and relative acoustic power imbalance at most 1e-7.

These are resolution-change bounds for one fixed model, not an asymptotic error estimate, uncertainty in a physical loudspeaker, universal acoustic standards or proof that no still-narrower feature exists. They do not establish the ranking of nearly tied candidates. Preserve any failing study and investigate it before freezing a new protocol.

## Run a study

Use the normal prerequisites and build the images. In a clean checkout, select an actual ranking entry (index zero is the leading result) and its single driver JSON record:

```sh
python scripts/validate_candidate_resolution.py prepare results/my-resolution-study \
  --run-dir results/my-design \
  --ranking results/my-design/outputs/auto/report/auto_ranking.json \
  --driver data/drivers-curated/18Sound/18sound-6nmb420.json \
  --f-low 800 --f-high 1600 --candidate-index 0

python scripts/validate_candidate_resolution.py solve results/my-resolution-study --jobs 3

python scripts/validate_candidate_resolution.py compare results/my-resolution-study
```

The host needs NumPy and the comparison stage additionally needs pandas and the dependencies of `horn-core`, `horn-drivers` and `horn-analysis`. The host launches the immutable solver image recorded by the originating run; a mutable latest tag is not used for execution. It mounts the frozen source read-only, captures each worker’s Python, numerical-library, MPI, PETSc and BLAS identity, and seals the successful host execution. Comparison verifies that execution seal, image identity and runtime record. The script loads production Python modules from the checkout it hashes. Only originating runs with the protocol baseline (10 mm mesh, 20 loft sections and 101 requested frequency points) are admitted. Preparation requires the completed run manifest, resolved settings, original ranking/STEP/response and a complete individual-driver database in the source snapshot. The originating package source hashes must match the current study. Ranking, resolved settings, STEP and response must match their completion-time output digests. Driver bytes must match that snapshot, and comparison must reproduce the stored ripple and average level before assessing resolution. Other originating settings require a separately designed protocol. Use `solve --jobs 3` to run three independent cases concurrently when resources permit. Workers shut down normally through a process executor; successful workers are not sent SIGTERM. Retain the full console log alongside the study. The host enforces a two-hour execution timeout and removes its own container on failure. Each directory is new; there is no partial-resume shortcut. Preparation and comparison require a clean Git checkout, including no untracked files or ignored Python source. Container solves verify the frozen source hashes before and after each case; the host comparison repeats the Git check after solving. Source changes or edited inputs invalidate comparison. A completed solve seals every STEP and response file; the comparison rejects missing or modified evidence and writes explicit passes/failures to `comparison.json`.

A pass supports numerical exploration of that candidate within the stated model and band. It does not turn an acoustic air-volume STEP into a mounting drawing, establish rear-load or phase-plug behaviour, or replace the two complete assembly references required by [the roadmap](SINGLE_HORN_ROADMAP.md).

## Executed 800–1600 Hz example

All six new solves, the originating response and all six comparisons passed for the hyperbolic/18Sound 6NMB420 candidate in the [worked example](WORKED_EXAMPLE.md). The full comparison band was approximately 566–2263 Hz.

| Refinement | Maximum output change (dB) | Maximum throat impedance change (dB) | Maximum phase change (degrees) |
| --- | --- | --- | --- |
| Original ranking → 201-point baseline | 0.004455 | 0.014793 | 0.08263 |
| Mesh 10 → 6 mm | 0.03969 | 0.11677 | 0.52536 |
| Mesh 6 → 4 mm | 0.01299 | 0.03757 | 0.16794 |
| Loft 20 → 40 sections | 0.0006692 | 0.001197 | 0.005445 |
| Loft 40 → 80 sections | 0.00001859 | 0.00001896 | 0.0007732 |
| Frequency 101 → 201 points | 0.002646 | 0.008821 | 0.04889 |

Target-band ripple changed by 0.001501 dB when doubling frequency sampling. Meshes contained 5,333 to 74,334 tetrahedra. All residuals were below 1e-15 and relative power imbalance below 1e-13. These changes establish stability of the evaluated model at these resolutions; they do not bound physical model error.

The [comparison JSON](../data/validation/candidate_resolution_reference.json), [complete raw archive](../data/validation/candidate_resolution_artifacts.tar.gz) and [identity manifest](../data/validation/candidate_resolution_manifest.json) retain every geometry, response, frozen input/source identity, solver image ID and console log. The manifest records the exact reproduction commit; use that commit to reproduce this historical study, or prepare a new study for newer source.

The final archived study completed with normal worker shutdown and parent exit code zero. Earlier studies are superseded in the manifest, including a numerical pass that required forced cleanup of a completed MPI worker. The originating run and its cached resume completed all 18 processes with output digests recorded.
