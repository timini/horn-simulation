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

Parameters are SI. The loader also accepts these unit-suffixed keys and converts them: `sd_cm2`, `vas_l`, `mms_g`, `le_mh`, `xmax_mm`, `exit_area_cm2`, `exit_diameter_mm`, `front_chamber_volume_cm3`, `front_chamber_area_cm2`, `rear_chamber_volume_cm3` and `power_aes_watts`. A value orders of magnitude off is rejected with a hint about the unit. Missing Mms is not replaced with the dry Mmd, and program power is never used as continuous power.

`driver_type` is `cone`, `compression` or `unknown`. A compression driver's throat is its exit (`exit_diameter_mm` or `exit_area_cm2`).

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
