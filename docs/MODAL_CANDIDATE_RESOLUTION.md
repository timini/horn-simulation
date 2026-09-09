# Resolution of the automatic modal candidate

The complete 800–1600 Hz [modal workflow](MODAL_APERTURE.md#executed-automatic-design) passes all six fixed resolution comparisons. This checks the actual ranked driver, geometry and response, followed by mesh, loft and frequency refinement. It does not qualify the physical driver interface.

The protocol and integrity checks are the same as the [default-model study](CANDIDATE_RESOLUTION.md). This study uses `modal_baffled` throughout, including all 16 aperture coefficients in the coupled on-axis observer. The target band is 800–1600 Hz and the complete comparison sweep is approximately 566–2263 Hz.

| Refinement | Maximum output change (dB) | Maximum impedance change (dB) | Maximum phase change (degrees) |
| --- | ---: | ---: | ---: |
| Original ranking → 201-point baseline | 0.004801 | 0.016741 | 0.08760 |
| Mesh 10 → 6 mm | 0.062262 | 0.145505 | 0.53024 |
| Mesh 6 → 4 mm | 0.020956 | 0.048957 | 0.17421 |
| Loft 20 → 40 sections | 0.000738 | 0.001257 | 0.00455 |
| Loft 40 → 80 sections | 0.0000240 | 0.0000274 | 0.000778 |
| Frequency 101 → 201 points | 0.002834 | 0.009952 | 0.05184 |

Target-band ripple changes are 0.002571 dB from the originating grid to the baseline, and 0.001700 dB from 101 to 201 frequencies on the finest geometry. Both use the stricter 0.2 dB / 2 degree frequency limits and 0.2 dB ripple limit. Spatial limits remain 0.5 dB / 5 degrees. All seven health records pass residual, grid, passivity and acoustic power checks; meshes span 5,333–74,334 cells.

## Evidence

The [identity manifest](../data/validation/modal_candidate_resolution_manifest.json), [comparison](../data/validation/modal_candidate_resolution_reference.json) and [raw archive](../data/validation/modal_candidate_resolution_artifacts.tar.gz) preserve original ranking/STEP/response, driver bytes, complete originating source snapshot, protocol, six new solves, exact image/runtime, clean host execution and all source/output seals. The earlier partial unpinned attempt is excluded.

The final union-grid comparison is a separately identified reanalysis of these sealed curves. The original comparison and both comparator source snapshots remain in the archive. No numerical result or tolerance was replaced. Reproduce the original solve at `reproduction_source_commit`, or the union-grid assessment at `analysis_revision`, following the [reanalysis procedure](CANDIDATE_RESOLUTION.md#executed-8001600-hz-example).

Mode-count refinement has its own [independent qualification](MODAL_APERTURE_VALIDATION.md) for this geometry over 800–1600 Hz. This study does not extend that independent reference comparison to the extra sweep margins, every profile, arbitrary ka or a finite baffle. New candidates need their own convergence checks.
