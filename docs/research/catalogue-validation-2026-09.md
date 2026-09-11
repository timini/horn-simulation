# Catalogue and annular-screen validation

11 September 2026. Scope: the catalogue update and 500–6,500 Hz reduced-order study in PR #94, stacked on the earlier annular-model work in PR #93.

## Local checks

- 144 catalogue, loader, scraper, validator and prescreen tests passed before review. Two additional regression tests exercise missing inductance rejection and valid zero-inductance resume; all 43 scraper tests pass after that fix.
- 314 horn-core and scripts/report tests passed.
- A live selected-size REDCATT import produced three records whose parameters exactly match the saved catalogue.
- 12,960 grid cases completed. Across 27 winning scenarios, maximum section-doubling change was 0.002355 dB, maximum frequency-grid interpolation change was 0.007012 dB, and relative real-power discrepancy was below 7.3e-14.
- Selected and Celestion STEP geometries passed the export/re-import checks. Response and section plots were visually inspected. Source/output hashes matched the generated manifest.
- `git diff --check` passed.

The repeated quarantine warnings during database tests are expected: suspect historical records remain in the repository and are explicitly refused by the simulation loader. Numerical plausibility tests for active records and separate quarantine regression tests serve different purposes. The 1,485 legacy-unverified records are not certified by a passing test suite.

## One review round

The first review attempt did not run because the installed CLI was incompatible with its configured model. The single allowed retry completed using GPT-5.5. It found one substantive issue: a scraped record could omit Le, be marked current, and subsequently fail the stricter loader.

The fix adds Le to the scraper's essential fields and resume checks, permits an explicit zero, rejects a missing value before writing, and tests both outcomes. All 43 scraper regressions pass. No second full review was requested.

## Scope of the predictions

The calculation is a plane-wave acoustic network with a rigid-piston motor. It does not validate cone breakup, phase equalisation, annular higher modes, HF output or summation, or the outer 15-inch horn. Factory-sealed drivers receive no additional rear-box compliance. The STEP files describe acoustic/housing concepts rather than a fabrication-ready driver adapter.

Remote workflow status is available on [PR #94](https://github.com/timini/horn-simulation/pull/94); this document does not substitute local results for GitHub checks. No merge has been authorised or performed.
