# Search and scoring

Auto mode turns a target band into ranked driver and horn pairs. It evaluates a fixed number of designs and reports the best ones it evaluated. It does not prove a global optimum, and close scores are not a reliable physical ordering.

## Stages

1. **Driver screen** (`horn_analysis.prescreen`). A driver is kept when its Fs is below 1.5 times `target_f_low`, its Fs/Qes (EBP) is at least 50, its nominal diameter is inside `min_diameter` and `max_diameter` when given (a diameter that cannot be parsed passes), and its diaphragm, or for a driver with an exit area (a compression driver) its exit, fits inside the mouth (the fixed mouth, or a ka = 1 mouth at `target_f_low`). Drivers without a Qes are dropped.
2. **Throat radii.** Unless `throat_radius` is fixed, throats are sampled up to ka = 3.83 at `target_f_high`, the first radial duct mode, above which a plane-wave screen means little. The samples are 0.3, 0.65 and 1.0 of that limit (or `throat_fractions`), plus throats matching each driver's diaphragm area and throats that put each driver's front chamber corner at `target_f_high`. Radii within 1 % are merged and at most five are kept. Throats above ka = 1.84 (the first asymmetric mode) are flagged `throat_above_first_asymmetric_mode`, and drivers with a diaphragm smaller than a throat are flagged `diaphragm_smaller_than_throat`. A throat you fix may go up to ka = 2 pi.
3. **Geometry grid** (`horn_core.geometry_designer`). Unfixed dimensions are derived from `target_f_low`: mouth radius around the ka = 1 radius c/(2 pi f_low), length from a quarter to a half wavelength. `num_mouth_radii` and `num_lengths` set the grid. `min_*` and `max_*` bounds replace that side of the derived range. With only a maximum, the search runs from 0.4 times the maximum up to it. When a cap is below the ka = 1 mouth or quarter-wave length, `design.json`, the ranking and the report state the expected low-frequency shortfall in octaves. Every grid point is crossed with every profile in `profiles`.
4. **Webster screen** (`horn_analysis.lem_prescreen`). Every geometry is coupled to every driver with a one-dimensional transfer-matrix model of the horn and scored in milliseconds. Geometries outside the radiation model's valid range, or impossible for their profile (for example a tractrix shorter than its natural length), are rejected here.
5. **FEM shortlist** (`horn_analysis.search.shortlist_geometries`). `lem_top_n` geometries go to FEM. The first half are taken by screen score alone. The best design of each profile not yet chosen then fills slots up to 80 % of the budget, and no profile takes more than half the budget while others have comparable candidates. The last 20 % prefer geometrically different near-ties. If every screened design is infeasible, the least-infeasible ones are sent.
6. **FEM solve and ranking.** Each shortlisted STEP is solved in `num_bands` parallel bands over the auto sweep, one octave below `target_f_low` to half an octave above `target_f_high`. Every pre-screened driver is coupled to every solved horn and scored.
7. **Refinement** (`horn_analysis.refine`). A deterministic coordinate search adjusts throat, mouth and length of the best design within `refinement_budget` extra FEM solves, keeping its profile. It starts from the best feasible FEM design or, if every FEM design failed, from the least infeasible one (fewest failed limits, then smallest ripple excess). When it ends on a bound that came from a heuristic, it widens that bound by 25 % and continues, at most twice. Bounds you set are never moved. `refinement/search_audit.json` records every step.

The report ranks the final rows and shows the top five in detail. `top_n` rows are kept.

## Hard limits

A pair is feasible only when all of these hold. Infeasible pairs score 0 and appear in the rejection table with their margins.

| Limit | Reason code |
| --- | --- |
| Target-band ripple at most `max_ripple_db` | `target_band_ripple_exceeded` |
| Lowest in-band level at least `min_output_db`, when set | `minimum_output_not_met` |
| Length and mouth within any `max_length` and `max_mouth_radius` | `maximum_length_exceeded`, `maximum_mouth_radius_exceeded` |
| Sd / throat area at most `max_compression_ratio` | `compression_ratio_exceeded` |
| Diaphragm no larger than the mouth | `driver_larger_than_mouth` |
| Target band inside the driver's usable band, when known | `outside_driver_usable_band` |
| Peak excursion at most Xmax at `voltage_rms` | `driver_excursion_exceeded` |
| Input power at most the rated power at `voltage_rms` | `driver_nominal_power_exceeded` |
| In-band excursion at most Xmax at `test_power_w`, when set | `test_power_excursion_exceeded` |
| In-band input power at most the driver rating at `test_power_w`, when set | `test_power_nominal_power_exceeded` |

A feasible pair with missing driver or interface data is labelled `insufficient_evidence` rather than rejected. The gaps are listed in `evidence_gaps`.

## Front chamber in the search

Levels and scores include the driver's front chamber, so the search prefers throats and horns that work with it. A cone driver without chamber data gets a geometric estimate. If a limit fails only because of that estimate, the design is kept, the reason moves to `estimated_chamber_reasons` and the gap `response_limited_by_estimated_front_chamber` is added. See [DRIVERS](DRIVERS.md#front-and-rear-chambers).

## Score

Feasible pairs are ordered by a weighted sum of four terms, each mapped linearly onto 0 to 1 on a fixed scale:

```text
score = w_sensitivity * clip((L_mean - 60) / 60, 0, 1)
      + w_ripple      * clip(1 - ripple / 12, 0, 1)
      + w_extension   * clip(0.5 + E / 2, 0, 1)
      + w_size        * clip(0.5 - log2(V / V_ref) / 8, 0, 1)
```

- `L_mean`: log-frequency mean observer level over the target band, dB.
- `ripple`: its peak-to-peak range, dB.
- `E = log2(f_low / f3_low) + log2(f3_high / f_high)`, in octaves. The f3 points are 3 dB below `L_mean`, found over the full sweep. E is negative when the response is already 3 dB down inside the band. When an f3 is a sweep bound, E is a lower bound and `extension_is_lower_bound` is set.
- `V = pi r_mouth^2 L`. `V_ref` is the cylinder of a ka = 1 mouth and a quarter-wave length at `f_low`.

Default weights are sensitivity 0.625, ripple 0.25, extension 0.125, size 0. Set others with `--score_weights`, for example `sensitivity=0.4,ripple=0.4,extension=0.2,size=0`. Weights are normalised to sum to 1.

Because each term is linear, the weights set fixed exchange rates. With the defaults, 1 dB less ripple is worth 2 dB of sensitivity, one octave more extension is worth 6 dB, and size does not matter. So a design 15 dB more sensitive is not beaten by one that is 3 dB flatter. The report and `ranked_results.json` (`scoring`) print the formula, weights, exchange rates and the origin of every search bound.

Rows within 0.02 of the best score are marked as near ties. This margin is a search policy, not a measured uncertainty.

## Search quality evidence

The shortlist was checked against exhaustive FEM over a fixed grid of 20 horns, three motors and four bands, for both the `flanged_piston` and `modal_baffled` models. In every case the screen kept all top-ten pairs and the winner, with zero score regret. See [VALIDATION](VALIDATION.md#search-and-workflow). This shows the screen does not lose good designs on that grid. It says nothing about physical ranking.
