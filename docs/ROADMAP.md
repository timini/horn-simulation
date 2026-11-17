# Roadmap

Written 2026-11-17, after [PR #99](https://github.com/timini/horn-simulation/pull/99) merged. This page says where the project is, what comes next and in what order, and gives the design of the planned Pareto design database and its web page. [STATUS](STATUS.md) keeps the short list of what works and the physical qualification plan. [DESCOPED](DESCOPED.md) lists what was dropped from the original brief.

Each item below is marked decided or open. Decided means the design is agreed and only needs building. Open means a question still needs an answer, usually from a measurement or a pilot run.

## Where the project is now

### Recent merges

| PR | Merged | What it added |
| --- | --- | --- |
| [#95](https://github.com/timini/horn-simulation/pull/95) | 2026-11-17 | The remaining gaps in the original brief, comparisons with openwind and bempp-cl, rewritten docs, a refreshed site |
| [#96](https://github.com/timini/horn-simulation/pull/96) | 2026-11-17 | A test fix for the plain range separator in the auto report |
| [#97](https://github.com/timini/horn-simulation/pull/97) | 2026-11-17 | Band join failures reject the candidate instead of stopping the run; widest flat band fields in every ranked row; the site rewritten for DIY builders |
| [#98](https://github.com/timini/horn-simulation/pull/98) | 2026-11-17 | The chosen driver and a to-scale horn profile for each example on the site |
| [#99](https://github.com/timini/horn-simulation/pull/99) | 2026-11-17 | A priced driver database and a price filter for auto mode |

Earlier work, up to September 2026, is summarised in [history](history/README.md).

### What works and has been checked

- The 3D finite element solver agrees with an independent FEM and motor model on shared meshes, and gives correct cavity modes. See [INDEPENDENT_SOLVER_VALIDATION](INDEPENDENT_SOLVER_VALIDATION.md).
- The `modal_baffled` mouth model matches an independent mode-matching code in an infinite baffle. See [MODAL_APERTURE](MODAL_APERTURE.md).
- The one-dimensional Webster screen (`horn_core.webster`) agrees with openwind on 35 of 35 cases across all seven profiles. See [CROSS_SOFTWARE_VALIDATION](CROSS_SOFTWARE_VALIDATION.md).
- The screen keeps the best FEM designs on a fixed test grid, for four bands and both radiation models.
- Refinement and the two worked candidates are converged in mesh, loft and frequency grid. See [CANDIDATE_RESOLUTION](CANDIDATE_RESOLUTION.md).
- bempp-cl confirms throat impedance of a complete horn at 5 of 5 frequencies.
- One measured horn (Post and Hixson exponential horn) matches in input impedance from ka 1 to 5. See [POST_HIXSON_VALIDATION](POST_HIXSON_VALIDATION.md).

### What is not validated

- Directivity. bempp-cl in a finite baffle disagrees with the infinite-baffle model, mostly because of baffle diffraction.
- Absolute level for a real driver on a real horn, and the order drivers are ranked in. No complete assembly has been measured.
- The Hajicek and Neely conical horn standard fails its frozen gates, and the 3D FEM differs from the 1D codes by up to about 1.5 dB between mouth ka 0.5 and 1.5 for reasons not yet explained. Both are kept as failures.

[VALIDATION](VALIDATION.md) has every record.

### Drivers and prices

From [DRIVER_COVERAGE](DRIVER_COVERAGE.md), generated on 2026-11-17:

| Item | Count |
| --- | --- |
| Usable driver records | 2,122 (plus 7 curated) |
| Records held back as incomplete | 396 |
| Records with Xmax | 1,989 (93 %) |
| Records with rated power and sensitivity | 833 (39 %) |
| Drivers with at least one exact or manual price | 367 |
| Priced drivers per region | EU 316, UK 189, CA 70, US 12 |

Prices come from six retailer and list sources, stored as dated JSON lines snapshots in `data/prices/`. Offers older than 30 days are treated as stale. Auto mode takes `--price_region` and `--max_price`. Price is shown with its source and date and the cheapest feasible design is marked, but price never enters the score. See [DRIVERS, prices](DRIVERS.md#prices).

### Search outputs added this month

- Widest flat band. Every ranked row reports the widest contiguous band within +-1, +-2 and +-3 dB (`horn_analysis.flat_band`, fields `flat_band_1db_*`, `flat_band_2db_*`, `flat_band_3db_*`). These are reported only; they do not change the score. See [SEARCH_AND_SCORING](SEARCH_AND_SCORING.md#widest-flat-band).
- Band joins. When two frequency bands disagree where they meet, the candidate is rejected as `band_boundary_mismatch` and the run carries on. The summary says the mesh is the cause. See [ACOUSTIC_CONTRACT](ACOUSTIC_CONTRACT.md#frequency-band-joins).

### Site

The [example site](https://timini.github.io/horn-simulation/) explains horns for DIY builders and shows five example briefs (PA midrange, hi-fi upper, low-mid, loud vocal, flat midrange) with their reports, chosen drivers and to-scale profiles. Three of them needed a finer `--mesh_size` to get a design. The low-mid run hit a Java memory limit.

## Near-term work

In this order. Each item says when it is done.

### 1. Finer mesh in the upper bands automatically

Today the user has to pass a smaller `--mesh_size` when band joins fail. In a local run of the 500 to 4000 Hz PA midrange brief at the default mesh, every candidate was rejected at the join near 3830 Hz. One candidate, `auto_tra_0150`, also showed a false spike of 123.75 dB at the join. (This run folder is not committed.)

The plan: when a candidate fails a join, re-mesh the two bands on each side of the join with a finer element size and solve those bands again, up to a fixed number of tries. Record the mesh used per band in the candidate's provenance.

Done when the five example briefs complete at default settings with no `--mesh_size`, the join check is unchanged, and a test covers the retry.

### 2. Widest-band objective

`--objective widest_band` with `--flat_band_tolerance_db`, as designed in [SEARCH_AND_SCORING, planned](SEARCH_AND_SCORING.md#planned-widest-band-objective). The search finds the widest +-N dB band inside an outer window, under a size cap.

Done when the objective runs end to end in auto mode, the screen, FEM ranking and refinement all use the same score, a size cap is required or warned about, and one example on the site uses it.

### 3. Maximum SPL metric

Report the highest level each design reaches before the driver hits Xmax or its rated power, inside the flat band. Excursion and power are already computed per row, so this is a new output, not new physics. It needs rated power, which 63 % of records have; sensitivity and power together are present for 39 %.

Done when ranked rows carry `max_spl_db` and which limit set it, rows without the data say so, and the report shows it.

### 4. More price sources

The US has 12 priced drivers, all Dayton list prices.

- Parts Express pages need JavaScript. A product feed from the retailer is needed. Open.
- loudspeakerdatabase.com currently returns the same page for every URL, so its health check fails and no driver pages are fetched. The owner has been emailed. Open.
- OpenISD data carries no licence. Ask its author before using it. Open.

Done when the US has prices for at least the drivers that the examples use, each new source has a recorded robots and terms check in `data/sources.json`, and [DRIVER_COVERAGE](DRIVER_COVERAGE.md) is regenerated.

### 5. Trust in the estimated front chamber

Most cone drivers have no front chamber data, so the search uses a geometric estimate from Sd (`horn_core.driver_interface.estimate_cone_front_chamber_m3`). Designs whose result depends on the estimate are flagged, but the estimate itself is unchecked.

Done when the estimate is compared with measured or drawn chamber volumes for a handful of real drivers, the error is stated in [DRIVERS](DRIVERS.md), and the high-frequency corner it implies is printed per design.

## Pareto design database and web page

### Goal

For every driver in the database, find the horns that give the best trade-offs, and publish them as a page you can search. You state what you need, for example "500 to 4000 Hz within +-2 dB, 110 dB, mouth under 300 mm, cheapest". The page returns the cheapest driver and horn pairs that meet it, or the closest misses and which limit each one fails.

Nothing like this exists publicly. Today a builder picks a driver by instinct and then iterates in Hornresp or ATH.

### Who it is for

- Traditional horn builders, passive crossovers: flatness matters most.
- DSP builders: maximum SPL and sensitivity matter more; coverage later.
- Small PA and DIY makers: price and size.

### Definitions (decided)

| Term | Definition |
| --- | --- |
| Flat band | The widest single continuous range whose peak-to-peak level is within the window. Two patches either side of a dip do not count. Computed by `horn_analysis.flat_band.widest_flat_band`. |
| Tolerance | +-1, +-2 or +-3 dB, meaning a 2, 4 or 6 dB window at any level. |
| Bandwidth | Octaves, log2(f_high / f_low). |
| Size class | Small: mouth diameter and length up to 150 mm. Medium: up to 300 mm. Large: up to 500 mm. Without a cap the best horn is always the biggest. |
| Drive and observer | 2.83 V, 1 m, on axis, infinite baffle. |
| Objective | Keep every design on the multi-metric Pareto front. No single weighted score. The widest band per tolerance and size class is one view of the front. |

### Metrics (decided, except where noted)

| Metric | Unit | Better | Source |
| --- | --- | --- | --- |
| Flat bandwidth at +-1, +-2, +-3 dB | octaves | higher | 1D model; FEM for confirmed rows |
| Sensitivity | dB at 2.83 V, 1 m | higher | model, mean over the flat band |
| Maximum SPL before Xmax or rated power | dB | higher | model (near-term item 3) |
| Throat air distortion at a reference SPL | % second harmonic | lower | standard throat intensity formula |
| Excursion at a reference SPL | % of Xmax | lower | model; a proxy for driver distortion |
| Minimum electrical impedance | ohm | higher | model |
| Mouth diameter and length | mm | lower | geometry |
| Driver price | currency of the region | lower | dated retailer offer |
| Coverage smoothness | to be defined | better | later, after directivity is validated |

The reference SPL for distortion and excursion is open; 105 dB at 1 m is the working value.

### Many metrics and ties (decided)

With nine metrics almost every design is on the front. Two rules keep the result useful.

1. Epsilon dominance. Design A beats design B only if A is at least as good as B on every chosen metric, within a tie margin, and better by more than the margin on at least one. The margins are 1 dB for SPL and sensitivity, 0.1 octave for bandwidth, 0.5 percentage points for distortion, and 5 % for price.
2. Requirements first. The user sets hard limits (band, tolerance, SPL at a distortion limit, size). Front members that meet every limit are listed, cheapest first. If none meets them, the closest misses are listed with the limit each one fails and by how much.

### Choose your own metrics (decided)

The page lets you pick any subset of metrics and set each one to optimise, ignore, or hard limit. It then shows the front for just that subset, worked out in the browser. No server is needed.

This works because a design on the front for a subset of metrics is either on the front for all metrics, or exactly tied on that subset with a design that is. So the build stores the all-metrics front with exact ties kept, and the page filters it. Two consequences:

- The stored front uses exact dominance, not epsilon dominance. Epsilon margins are applied only in the browser, so no design the browser might need is dropped at build time.
- Hard limits on stored metrics are safe, because a design that dominates another also meets any limit the other meets. Filters on things that are not metrics (maker, profile, nominal size) are not safe in that sense: a design hidden only because another maker's driver beat it would stay hidden. Open: either add those as group keys and keep a front per group, or show a note when such a filter is used.

### Method

1. Horn library. About 50,000 horn shapes per run: the seven profiles crossed with throat, mouth and length grids for each size class, from `horn_core.geometry_designer` and `horn_core.profiles`. Each shape is solved once with the transfer matrix model (`horn_core.webster.compute_horn_transfer_tmm`). The horn alone does not depend on the driver. Measured at about 7 ms per shape for 200 frequencies, so about 6 minutes on one core.
2. Driver pairing. Each driver's electro-mechanical model, with front and rear chamber, is coupled to every compatible horn using the same code as auto mode (`horn_core.driver_interface.driver_operating_point` through `horn_analysis.evaluation.evaluate_response`). Compatible means the diaphragm fits the mouth and the compression ratio is within the auto-mode limit, so most pairs are never computed.
3. All tolerances at once. `horn_analysis.flat_band.flat_band_fields` reads +-1, +-2 and +-3 dB from the same response.
4. Keep the front per driver. For each driver and size class, keep only the horns on that driver's own front. A design on the global front is always on its driver's front, so this loses nothing and keeps memory small.
5. Refine. A short local search around each front member, reusing the coordinate search in `horn_analysis.refine` with the Webster model instead of FEM.
6. Confirm with FEM. The solver runs on front horns, about 40 s per horn locally. The 1D model is about 1 dB off near the mouth, so +-1 dB rows are published only after FEM. Every row carries its check level.

Compute estimates (open until the pilot measures them):

| Step | Estimate |
| --- | --- |
| Horn library | about 6 minutes on one core |
| Driver pairing, coupling only | the agreed design said well under an hour; it was estimated with fewer drivers and needs re-measuring at 2,122 |
| Flat band search | about 1.5 ms per row today, which is too slow for tens of millions of pairs; needs a vectorised version over many responses, or a cheap pre-check (for example peak-to-peak over the target window) before the full search |
| FEM confirmation | a few hours for the top 100 drivers; 1 to 2 days for all front horns |

### Data files (decided, format details open)

One file per size class, versioned, under `data/designs/`:

```text
data/designs/
  v1/
    manifest.json             code version, inputs hash, build date, row counts
    horns-medium.jsonl.gz     horn library: one row per horn shape
    designs-medium.jsonl.gz   one row per driver x horn on the per-driver front
    front-medium-UK.json.gz   site file: global front for one region
```

Plain gzipped JSON lines need no new dependency. Parquet or Arrow would be smaller and faster to load; adding `pyarrow` is open.

Columns in `designs-<class>.jsonl.gz`:

| Column | Unit | Notes |
| --- | --- | --- |
| `design_id` | | hash of driver_id, horn_id and code version |
| `driver_id` | | key into `data/drivers` |
| `horn_id` | | key into the horn library |
| `size_class` | | small, medium, large |
| `profile` | | one of the seven profiles |
| `throat_radius_mm`, `mouth_radius_mm`, `length_mm` | mm | geometry |
| `flat_band_1db_low_hz`, `_high_hz`, `_octaves`, `_mean_level_db`, `_is_lower_bound` | Hz, octaves, dB | same for 2 dB and 3 dB; names match the ranking rows |
| `sensitivity_db` | dB at 2.83 V, 1 m | |
| `max_spl_db`, `max_spl_limit` | dB, text | `xmax` or `power`; null when driver data is missing |
| `throat_distortion_pct` | % | at the reference SPL |
| `excursion_pct_xmax` | % | at the reference SPL |
| `min_impedance_ohm` | ohm | |
| `front_chamber` | | `supplied`, `zero_phase_plug` or `estimated` |
| `check_level` | | `analytic` or `fem` |
| `fem_flat_band_*` | | FEM values, when `check_level` is `fem` |
| `evidence_gaps` | | list, from the existing evidence gap code |
| `driver_data_source`, `driver_data_date` | | from the driver record provenance |
| `code_version`, `inputs_hash` | | git commit and the hash of the driver, price and grid inputs |

Price is not stored in this file. It is joined at build time of the site file (next section).

The site file holds only the global front for one region, with the columns the page needs, in a columnar layout (one array per column) to make it small. Each row there adds `price`, `currency`, `price_source`, `price_date` and `buy_url`. The target is tens of thousands of rows and a few MB per file after compression.

### How it is built (decided)

A new script, `scripts/build_design_db.py`, with a `just design-db` recipe. A Nextflow mode is not needed for the 1D stages, which run in minutes on one machine. FEM confirmation reuses the existing single-mode Nextflow solve, one task per front horn, so it can run in parallel and resume.

Stages:

1. Read drivers through the same loader and validator as auto mode. Drivers that fail validation are skipped and counted.
2. Build or reuse the horn library (cached by its inputs hash).
3. Pair, compute metrics, keep each driver's front, write `designs-<class>.jsonl.gz`.
4. Join prices and write the site files.
5. Optionally queue FEM for front horns and merge the results back as `check_level` fem.

The manifest records the code version, the inputs hash and the counts at each stage, so one command reproduces the table.

### Prices refresh without new acoustics (decided)

Acoustic rows do not change when prices do. Stage 4 alone joins the newest offers by `driver_id`, using `horn_drivers.prices.price_snapshot` with the same rules as auto mode: exact or manual matches only, list prices only when no shop offer exists, offers older than 30 days dropped. It then recomputes the front for each region, which takes seconds. Every price on the page shows its source and date. Drivers with no price in a region stay on the page with price shown as unknown and are treated as worst on price.

### Web page (decided, details open)

A static page on the existing GitHub Pages site, `site/designs.html`. No server, no build step, no framework. Plain JavaScript, with at most one small pinned library for the parallel coordinates plot if writing it by hand proves slow. Data is fetched as the gzipped columnar JSON above; the browser decompresses it.

Front recompute in the browser:

- Apply hard limits first, which usually removes most rows.
- Sort the rest by the first optimised metric, then run a block nested loop: compare each row only with the current front, and drop it if a front member beats it by more than the margins. Cost is rows times front size times metrics. For 30,000 rows and a front of a few thousand this is around 100 million comparisons in the worst case, which is under a second in a Web Worker on a laptop. The pilot measures it on a phone.
- Rows within the margins of each other are ties. Ties are all kept and shown together, sorted by price.

The whole query (size class, region, metrics, their modes and limits) is stored in the URL, so a query can be shared or bookmarked.

Accessibility and phones: the results list is a real table that works with a keyboard and a screen reader, and it is the main view. The parallel coordinates plot is an extra; every filter it offers is also a form field. On a phone the plot sits below the list and scrolls sideways inside its own box; the page itself does not.

### Page flow

1. Pick a size class and a price region.
2. For each metric choose optimise, ignore, or hard limit with a value. Defaults: bandwidth at +-2 dB optimise, price optimise, everything else ignore.
3. The results list shows front members that meet the limits, cheapest first. If none do, it shows the closest misses with the failing limit and the shortfall.
4. The parallel coordinates plot shows the same rows, one axis per chosen metric. Dragging a range on an axis filters the list.
5. Clicking a row opens its detail page.

```text
+------------------------------------------------------------------+
| Size [medium v]   Region [UK v]   Tolerance [+-2 dB v]           |
+------------------------------------------------------------------+
| Metric            Mode          Limit                            |
| Bandwidth         (o) optimise  ( ) limit [   ] ( ) ignore       |
| Max SPL           ( ) optimise  (o) limit [110] dB ( ) ignore    |
| Mouth diameter    ( ) optimise  (o) limit [300] mm ( ) ignore    |
| Price             (o) optimise  ( ) limit [   ] ( ) ignore       |
| Distortion        ( ) optimise  ( ) limit [   ] (o) ignore       |
+------------------------------------------------------------------+
| 14 designs meet your limits. Cheapest first.                     |
| Price   Driver          Band (+-2 dB)    SPL   Mouth  Check      |
| 000 GBP Maker Model     000 to 0000 Hz   000   000    FEM        |
| ...                                                              |
+------------------------------------------------------------------+
| Parallel coordinates: Band | SPL | Mouth | Price                 |
|   one line per design, drag on an axis to filter                 |
+------------------------------------------------------------------+
```

The row in the wireframe shows the layout only; the values are placeholders.

Detail page per design:

- driver, horn profile and dimensions, with the to-scale horn drawing used on the [examples page](https://timini.github.io/horn-simulation/examples.html);
- predicted response with the flat bands shaded, and electrical impedance;
- downloads: STEP of the air volume, profile CSV, response CSV;
- buy links with retailer, price and the date it was read;
- check level, front chamber source and evidence gaps.

STEP files are generated on demand by the existing geometry code from the stored dimensions, or pre-built for front rows only. Which one is open; pre-building all of them could make the site too large.

### Keeping it honest (decided)

- The page carries the same experimental label as the reports.
- Each row shows its check level: Analytic or FEM. +-1 dB rows appear only with FEM.
- Rows whose result depends on an estimated front chamber are flagged, with the implied corner frequency.
- Driver data source and date are shown on the detail page.
- Closest misses are shown instead of an empty page.
- An invitation to send measurements, and failures published next to passes.

### Phases

| Phase | Scope | Done when |
| --- | --- | --- |
| 1. Pilot | Top 50 to 100 priced drivers, three size classes, FEM confirmed | Page published; every row's FEM band within its stated tolerance; analytic against FEM band difference reported; posted on diyAudio for feedback |
| 2. Full database | Every driver that passes the validator; analytic rows for all, FEM for as many as compute allows | Every row has a check level; one command rebuilds the table with code version and inputs recorded |
| 3. Coverage | Coverage column | Directivity check against an independent exterior solver passes |
| 4. Measured spot checks | At least two rows built and measured | Agreement reported, including failures |

Acceptance for the feature as a whole: one command rebuilds the table and records the code version and inputs; a unit test shows the band finder picks one continuous band (already true of `widest_flat_band`); for the pilot, the analytic and FEM bands are both reported and the FEM value is the one published; each row links STEP, profile CSV, response and impedance.

### Limits to state on the page

- Driver data is mostly datasheet or scraped Thiele-Small values, with unchecked accuracy.
- Few compression drivers have exit and diaphragm data, so the pilot is cone drivers.
- No measured validation of a complete assembly yet.
- Round horns in an infinite baffle only.

### Risks

| Risk | Mitigation |
| --- | --- |
| Wrong driver data | Validator in CI, provenance per field, consistency checks, conflicts kept |
| Estimated front chamber | Flag it, let builders enter real volumes, show the corner frequency |
| 1D error near +-1 dB | +-1 dB rows only after FEM |
| People build without checking | Experimental label, invite measurements, publish failures |
| FEM compute | Pilot first; top drivers only, or rented machines for the full run |
| Price data goes stale or a retailer blocks access | Prices dated, 30 day stale rule, price refresh separate from acoustics |

### Open questions

- Reference SPL for distortion and excursion.
- Parquet or Arrow instead of gzipped JSON.
- Non-metric filters (maker, profile, nominal size): front per group, or a warning.
- STEP on demand or pre-built.
- Flat band search speed at full scale.
- Driver pairing time at 2,122 drivers.

## Later

- Coverage column, once directivity is checked against an independent exterior solver.
- Measured spot checks of database rows, and the complete-assembly comparison in [STATUS](STATUS.md#remaining-work).
- Compression drivers, once exit geometry and phase plug data exist for enough of them.
