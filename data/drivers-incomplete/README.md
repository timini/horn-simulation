# Compression drivers with incomplete public data

`compression-drivers.json` lists compression drivers whose manufacturer pages or datasheets give the exit diameter, frequency range, power rating and sensitivity, but not the full Thiele-Small set (Fs, Re, Bl, Sd, Mms, Le and the Qs) that the driver model needs. Each record lists the values that were found, with the source URL and retrieval date, and names the missing fields in `missing_ts_fields`. No value was estimated.

These records are not in a driver database directory and the pipeline cannot load them. They record what is still needed: measure the missing parameters (an impedance sweep with an added-mass or closed-box step gives Fs, Re, Le, Qms, Qes, Bl and Mms; Sd needs the diaphragm dimensions) or obtain them from the manufacturer, then enter the driver with `horn-add-driver`.

The manufacturer folders (`18Sound/`, `BC/`, ...) hold version 2 records from the driver-base import, mostly compression drivers, with `incomplete_reason` naming the missing fields. They carry the power, sensitivity and band driver-base read from the manufacturer page.

Four compression drivers with complete published parameters are in `data/drivers-curated/Peerless`.

## Scraped records removed from data/drivers

`removed-scraped-records.json` lists scraped records taken out of `data/drivers`, so that the models are not silently lost:

- `copied_pages_not_restored`: 376 Faital, B&C and 18Sound records whose stored page was a copy of another driver's page. They were removed in 396f2d2. A re-scrape was attempted on 2026-11-15, but the site answered HTTP 429 to every request, so none has been restored. To restore them, run `horn-scrape-drivers --manufacturers Faital,BC,18Sound` when the site responds; the scraper refuses a page that names a different driver.
- `inconsistent_source_values`: 7 records whose own values contradict each other (Qts far from Qms*Qes/(Qms+Qes), or a page shared by two unrelated models). Their scraped values and the validator's reason are kept. Each needs the manufacturer datasheet.
