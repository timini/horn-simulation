# driver-base import

Data from driver-base (https://github.com/JonathanDCohen/driver-base), licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).

Imported from commit `f0fc602b2fb4154b20b90d4cae03fa82296b249a` (`web/drivers.json`, generated 2026-11-06) on 2026-11-17 by `horn-sources run driver-base`.

Changes made in this project: units converted to SI; values that driver-base derived or inferred (for example program power as twice AES power, or 1 W sensitivity computed from 2.83 V) left out; records matched against `data/drivers` and `data/drivers-curated`; power, sensitivity and band added to matched scraped records that lacked them; records whose Fs, Qes or Vas contradict their other values sent to `data/drivers-incomplete` (or, for a Cms that alone disagrees, the Cms left out); driver-base warn flags copied into the record notes. No existing value was changed.

Files written by the import:

- `summary.json`: counts from the last run.
- `matches.json`: each driver-base record matched to an existing record, how it was matched, the fields that disagree and the fields added.
- `prefix_review.csv`: driver-base records whose model name extends an existing one by a suffix that is not an impedance code (for example "EVO"), with matching Fs, Re, Bl, Sd and Mms. They may be a revision or the same driver, so they are neither filled nor imported until someone decides.
- `conflicts.csv`: per-field disagreements beyond 5 % (10 % for Le and Xmax). The existing value is kept. Where the existing record notes that `mms_kg` may be the dry mass Mmd, the row says so.
