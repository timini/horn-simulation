# Modal-radiation search qualification

The fast screening stage retains all exhaustive top-ten results and every winner in the same fixed 20-horn, three-motor, four-band benchmark used for the earlier local-radiation model. Final rankings use the complete 16-mode aperture observer. Screening uses the labelled plane-mode approximation. No score or shortlist setting was tuned using the new FEM results.

| Target band | Role | Model-feasible pairs | Top-ten recall | Winner retained | Score regret |
| --- | --- | ---: | ---: | --- | ---: |
| 800–1600 Hz | Development | 32 | 100% | Yes | 0 |
| 1000–1200 Hz | Held out | 40 | 100% | Yes | 0 |
| 1200–2000 Hz | Held out | 38 | 100% | Yes | 0 |
| 600–1000 Hz | Held out | 27 | 100% | Yes | 0 |

All four cases pass the frozen gates: at least 90% top-ten recall, retained winner, score regret no more than 0.02 and at least ten model-feasible pairs. The grid contains conical/exponential profiles, 22.5 mm throat radius, 45/70 mm mouth radii and lengths 40/55/70/85/100 mm. Every horn is solved at 100 logarithmic frequencies from approximately 424 to 2828 Hz with an 8 mm mesh; all curves pass finite-data, grid, residual, convergence and power-balance checks before ranking. Twenty horn simulations generate all 240 driver/horn/band comparisons.

This tests acceleration against exhaustive FEM **at the same fixed resolution**. It does not prove mesh convergence for this whole grid, a global optimum, seven-profile coverage or physical driver suitability. The three manufacturer-parameter motors retain unverified interface and moving-mass/rear-load status. A held-out band here is held out from screening development, not an independent physical measurement.

## Evidence and reproduction

The [identity manifest](../data/validation/modal_search_manifest.json) records the executed revision and immutable solver image. The [archive](../data/validation/modal_search_artifacts.tar.gz) contains all 20 STEP files, 20 raw responses, all screening/exhaustive rankings, runtime identity, source snapshot, host execution/seals, exact host runner and console log. The [summary](../data/validation/modal_search_reference.json) records the gates and outcomes. An earlier attempt failed before any simulation because a diagnostic helper was absent; it is explicitly excluded.

At the recorded source revision, use the pinned image and archived host runner, adjusting only the workspace/output paths. The runner checks clean source, freezes package/driver/harness hashes, mounts source read-only, captures actual runtime and seals outputs after execution. The benchmark entry point inside that environment is:

```sh
python scripts/benchmark_search.py --output-dir /study --radiation-model modal_baffled
```

The archive regression test reconstructs all 240 rankings from the sealed raw FEM responses and rechecks every selection gate. It does not depend on the stored pass flag alone.
