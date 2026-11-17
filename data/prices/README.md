# Driver price snapshots

Prices live here, not in the driver records.

- `<source_id>/<YYYY-MM-DD>.jsonl`: one offer per line, appended by a retailer adapter
  (`horn-sources run <source_id>`). Fields: `driver_id`, `source_id`, `url`, `sku`, `price`,
  `currency`, `vat_included`, `region`, `in_stock`, `retrieved_on`, `match_confidence`
  (`exact`, `manual` or `fuzzy`). The price is stored as listed, shipping excluded.
- `fx/ecb-<rate date>.json`: European Central Bank euro reference rates, fetched with
  `horn-sources fx`. Prices are converted only at ranking time, with the newest rate file
  dated on or before the ranking date.

Only `exact` and `manual` matches are used in a run. An offer older than 30 days on the run
date is not used, and `horn-validate-drivers --prices data/prices` warns about it.
Field rules and the loader are in `packages/horn-drivers/src/horn_drivers/prices.py`;
see [docs/DRIVERS.md](../../docs/DRIVERS.md).
