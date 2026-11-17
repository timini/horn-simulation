# Drivers

Drivers are JSON files, one per driver, under version control. Choose a database with `--drivers_db`.

| Folder | Contents |
| --- | --- |
| `data/drivers` | 1,454 drivers (1,342 labelled cone, 112 of unknown type) scraped from loudspeakerdatabase.com across nine manufacturers (18Sound, B&C, Beyma, Celestion, Ciare, Dayton, Eminence, LaVoce, SICA). The default. |
| `data/drivers-curated` | Seven hand-checked records from manufacturer pages and datasheets: three cone drivers (18Sound 6NMB420, Celestion AF3010 and TF0512HE) and four Peerless compression drivers. See [its README](../data/drivers-curated/README.md). |
| `data/drivers-incomplete` | Not loadable. Compression drivers whose public data lack a full Thiele-Small set, and scraped records removed from `data/drivers` with the reason. See [its README](../data/drivers-incomplete/README.md). |

None of these records is a physically qualified driver. Every one lacks a verified throat interface and a separate diaphragm mass and rear load, so rankings label them `insufficient_evidence`.

## Record format

```json
{
  "driver_id": "celestion-af3010-16",
  "manufacturer": "Celestion",
  "model_name": "AF3010 (16 ohm)",
  "nominal_diameter": "3in",
  "driver_type": "cone",
  "parameters": {
    "fs_hz": 128.9, "re_ohm": 10.48, "bl_tm": 5.98, "sd_m2": 0.003318,
    "mms_kg": 0.00355, "le_h": 0.00026, "qms": 6.096, "qes": 0.843, "qts": 0.74,
    "xmax_m": 0.00275, "power_w": 35.0, "nominal_impedance_ohm": 16.0
  },
  "usable_f_low_hz": 120,
  "usable_f_high_hz": 18000,
  "parameter_source": "https://celestion.com/product/af3010-2/",
  "interface_model": null,
  "record_schema_version": 1,
  "provenance": {"retrieved_on": "2026-09-07", "supplied_fields": ["..."], "derived_fields": []}
}
```

Parameters are SI. `power_w` is the nominal (AES, RMS or continuous) rating and `peak_power_w` the peak or program rating; the loader never uses one for the other. `sensitivity_db_1w1m` is the datasheet sensitivity in dB SPL at 1 W and 1 m; it is shown in the report but not used in the score. The validator warns when a cone driver's stated sensitivity is more than 4 dB from the Thiele-Small reference efficiency (sheets mix the 1 W and 2.83 V conventions). The usable band may sit at the top level or inside `parameters`.

The loader also accepts these unit-suffixed keys and converts them: `sd_cm2`, `vas_l`, `mms_g`, `le_mh`, `xmax_mm`, `exit_area_cm2`, `exit_diameter_mm`, `front_chamber_volume_cm3`, `front_chamber_area_cm2`, `rear_chamber_volume_cm3` and `power_aes_watts`. A value orders of magnitude off is rejected with a hint about the unit. Missing Mms is not replaced with the dry Mmd, and program power is never used as continuous power.

`driver_type` is `cone`, `compression` or `unknown`. A compression driver's throat is its exit (`exit_diameter_mm` or `exit_area_cm2`).

### Record schema version 2

New records from the source adapters use `"record_schema_version": 2`. The loader reads versions 1 and 2. Version 2 keeps every version 1 field and adds:

```json
{
  "identity": {"manufacturer": "BC", "model": "8NDL51",
               "variant": {"impedance_ohm": 8.0, "coil": null}, "aliases": ["8NDL51-8"]},
  "field_provenance": {
    "fs_hz": {"source_id": "bc", "url": "https://www.bcspeakers.com/...", "retrieved_on": "2026-11-17", "method": "page"}
  }
}
```

- `identity.manufacturer` must be a canonical name from `horn_drivers.identity.MANUFACTURERS`, which matches the folder names (`BC` for B&C, `18Sound` for Eighteen Sound). Aliases there cover retailer spellings. `aliases` holds retailer SKUs and model spellings for matching.
- `field_provenance` gives the source of each stored field. `method` is `page`, `pdf`, `derived` or `manual`. Every stored Fs, Re, Bl, Sd, Mms and Le needs an entry (an error otherwise); other fields without one get a warning. An entry may list `conflicts`, the values that lost to a better source.
- `horn_drivers.records.make_record` builds a version 2 record and checks its structure.

## Sources and polite fetching

`data/sources.json` lists every site the project may fetch from: id, kind (`manufacturer`, `aggregator`, `retailer`, `reference`), base URL, robots.txt snapshot date, terms URL, approval status, minimum delay and cache window. Only `approved` sources run. `horn-sources robots` refreshes the robots.txt copies in `data/sources/robots/`.

All fetching goes through `horn_drivers.fetch.PoliteFetcher`:

- robots.txt is read for each host and disallowed paths are refused. A URL must pass both Python's `urllib.robotparser` and an RFC 9309 matcher that understands `*` and `$` (robotparser ignores wildcard rules such as `Disallow: *?filter=`);
- one request at a time per host, at least the source delay (5 s for retailers, 10 s where robots.txt sets Crawl-delay 10, 30 s for loudspeakerdatabase.com) plus jitter;
- every kept response is cached under `results/cache/<source>/` and not refetched inside the cache window (7 days for prices, 90 for Thiele-Small data);
- two consecutive 429 or 403 responses stop the source for the run;
- the User-Agent names the project and its public repository.

It does not work around bot protection. A site that blocks polite access is marked `blocked` in the registry. Adapters store facts only (numbers, units, URL, date, SKU), never descriptions or images.

```sh
.venv/bin/horn-sources list
.venv/bin/horn-sources fx                       # ECB rates into data/prices/fx/
.venv/bin/horn-sources run thomann-uk --limit 20
```

Adapters live in `horn_drivers/sources/`, one module per source, and register with `@register_adapter`. A `DriverSource` writes version 2 records to `data/drivers` (or `data/drivers-incomplete` with the reason); a `PriceSource` writes exact matches to `data/prices` and everything else (fuzzy or unmatched) only to a review CSV under `results/price-review/`.

## Prices

### Retailer adapters

| Source | Region, currency | How it reads prices |
| --- | --- | --- |
| `thomann-uk`, `thomann-de` | UK GBP, EU EUR, VAT included | Brand-and-size loudspeaker category pages from the sitemap; each holds a JSON block per product. A product page is fetched only when the category entry has no impedance and the model is in the database. The German site is used for EUR because the English pages on thomann.de sit under the disallowed `/intl/`. |
| `lautsprechershop` | EU EUR, VAT included | Maker category pages from the English sitemap; one block per product with its order number and gross price. |
| `bax-shop-uk` | UK GBP, VAT included | Loudspeaker product pages from the pro audio sitemap. |
| `solen` | CA CAD, tax excluded | The shop's listing endpoint, which answers JSON when asked (`Accept: application/json`). Pair prices go to review. |
| `madisound` | US USD, tax excluded | Product pages from the sitemap (JSON-LD offer). Stocks almost none of the pilot makers. |
| `dayton-msrp` | US USD | MSRP on Dayton Audio's own product pages (pro audio category); a list price, not a shop offer. |

Not scraped: Parts Express (pages need JavaScript; ask for a feed), the Eminence store (its `/agents.md` asks automated clients to sign requests with Shopify Web Bot Auth, which this project does not do), Loudspeakerfreaks (robots.txt disallows everything), Simply Speakers (robots.txt answers 403).

Matching (`DriverIndex.match`): an offer is `exact` when the canonical manufacturer, the normalised model and the impedance agree with a database record and the listing's own title or product block names both maker and model. Impedance comes from the model name ("12NW76 8 Ohms"), a labelled spec line ("Impedance: 8 Ohm", never Re) or the product name. When the listing's model is not in the database, a database model named in the title is used, unless a revision suffix follows it ("Beta-12A-2"). Repair kits, diaphragms and other parts never match. A different impedance variant from the one in the database is not matched.

`horn-sources coverage --manufacturers 18Sound,BC,Celestion,Beyma,Eminence --driver-type cone` counts drivers with at least one exact price, per region, source and maker.

### Snapshots

Prices are kept out of the driver records, in `data/prices/<source>/<YYYY-MM-DD>.jsonl`, one offer per line: `driver_id`, `source_id`, `url`, `sku`, `price`, `currency`, `vat_included`, `region`, `in_stock`, `retrieved_on`, `match_confidence`. The price is stored as listed, shipping excluded.

With `--price_region`, a run takes for each driver the latest offer per source in that region, ignores offers older than 30 days and matches that are only `fuzzy`, converts other currencies with the stored ECB rate file (`data/prices/fx/`, newest on or before the run date) and keeps the cheapest, preferring offers in stock. `--max_price` then drops drivers above the limit or without a current offer. Price is shown, never scored. `horn-validate-drivers --db data/drivers --prices data/prices` checks the snapshot files and warns about offers older than 30 days, unknown drivers and fuzzy matches.

## Add a driver from a datasheet

Give each value in the unit the datasheet uses:

```sh
.venv/bin/horn-add-driver --db data/drivers-curated --manufacturer Acme --model "AC10 8ohm" --type cone \
  --fs-hz 45 --re-ohm 5.6 --bl-tm 15.2 --sd-cm2 346 --vas-l 52 --le-mh 0.9 \
  --qms 5.1 --qes 0.38 --qts 0.354 --xmax-mm 4 --power-w 300 \
  --f-low-hz 50 --f-high-hz 3000 --source-url https://example.com/ac10.pdf
```

Mms is derived from Vas, Sd and Fs when not given. Compression drivers take `--exit-diameter-mm`. Adapter or phase-plug air goes in `--front-chamber-cm3` (0 means none), with `--front-chamber-area-cm2` when it acts as a short duct. A sealed rear volume goes in `--rear-chamber-cm3`. The command converts to SI, runs the loader and validator, records the source URL and date, and writes nothing if either check fails. `--dry-run` prints the record instead.

## Validate a database

```sh
just validate-drivers
.venv/bin/horn-validate-drivers --db data/drivers
```

Errors include Qms, Qes and Qts that disagree by more than 10 %, records the loader would reject, garbled model names and the same parameter set across unrelated models (a copied page). CI requires zero errors on both shipped databases. Warnings flag missing optional data.

## Scrape drivers

```sh
uv sync --all-packages --extra test      # the test extra includes requests and beautifulsoup4
.venv/bin/horn-scrape-drivers --db data/drivers --manufacturers Eminence --state results/scrape-state.json
```

The scraper resumes from its state file and skips valid current records; `--refresh` refetches them. Each file is replaced atomically, so a failed refresh keeps the last good record. Discovery failures, incomplete pages and failed fetches return a failure status. `--patience-hours` waits through site outages for a bounded time.

Pages must name the requested driver, or they are refused. Each record stores the source URL, retrieval date and schema versions, plus the site's Rms, Cms, power, sensitivity and nominal frequency range (as `usable_f_low_hz` and `usable_f_high_hz`, labelled nominal). The driver category comes from the site's label and is recorded in `type_source`; diaphragm area alone cannot tell a compression driver from a small cone. Values entered by hand survive a refresh unless the source now supplies them.

## Front and rear chambers

The coupling (`horn_core.driver_interface`) puts these between the motor and the throat:

- **Front chamber**, the air between diaphragm and throat. Set `front_chamber_volume_m3` (or `_cm3`). Alone it acts as a shunt compliance, with a corner near f = c S_t / (2 pi V). Add `front_chamber_area_m2` (or `_cm2`) to model it as a short duct instead. `0` means an ideal area transformer with no chamber.
- **Rear chamber**, a sealed volume behind the diaphragm. Set `rear_chamber_volume_m3` (or `_cm3`). It adds the air-spring stiffness in series with the suspension. Unset keeps the free-air rear load already in the Thiele-Small parameters.

A cone driver with no front chamber value gets an estimate: the air in its own cone, with depth a quarter of the effective diameter, closed by a flat adapter plate, as a short duct of area Sd. Rows using it carry `front_chamber_estimated: true` and the gap `front_chamber_volume_estimated`. This is not a bound: a spacer adds volume, and a phase plug or cone filler removes it. An estimated chamber never rejects a design by itself. A limit that fails only with the estimate moves to `estimated_chamber_reasons`, the gap `response_limited_by_estimated_front_chamber` is added, and `without_front_chamber` holds the ripple and level with no chamber. Other driver types without a value use the ideal transformer, and get the gap `front_chamber_volume_unknown` when the diaphragm is larger than the throat. Supply the real adapter volume to replace the estimate.

## Evidence gaps

A feasible pair is still labelled `insufficient_evidence` when any of these apply. They are listed in `evidence_gaps`:

| Gap | Meaning |
| --- | --- |
| `driver_usable_band_unknown` | No usable frequency range recorded. |
| `mechanical_damping_unknown` | No Rms. |
| `driver_parameter_provenance_missing` | No source URL. |
| `driver_interface_unverified` | No verified throat adapter or phase plug model. |
| `moving_mass_and_rear_load_not_separated` | Mms not split into diaphragm mass and rear air load. |
| `maximum_output_limits_unknown` | No Xmax or rated power, so those checks could not run. |
| `front_chamber_volume_estimated`, `front_chamber_volume_unknown`, `response_limited_by_estimated_front_chamber` | Front chamber policy above. |
| `diaphragm_smaller_than_throat` | Sd below the throat area: the interface expands, which the plane-wave transformer does not describe. |
| `front_chamber_not_lumped_in_band` | The band reaches ka = 1 for the diaphragm radius, where a single chamber pressure is a poor approximation. |
| `throat_above_first_asymmetric_mode` | Throat ka above 1.84 at `target_f_high`. |

Every shipped record has at least the interface and mass gaps.
