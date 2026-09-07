# Manufacturer-sourced numerical comparison set

These three identified impedance variants were transcribed from the manufacturers' pages on 7 September 2026. Every JSON records its source, units, supplied fields and unresolved assumptions. The original large scraped database remains available; this set can be selected with `--drivers_db data/drivers-curated`.

| Driver | Variant | Primary source |
| --- | --- | --- |
| 18Sound 6NMB420 | 8 ohm | https://www.eighteensound.it/es/products/lf-driver/6-5/8/6NMB420 |
| Celestion AF3010 | 16 ohm | https://celestion.com/product/af3010-2/ |
| Celestion TF0512HE | 8 ohm | https://celestion.com/product/tf0512he/ |

This is a traceable parameter set, **not a physically validated horn-driver shortlist**. Manufacturer nominal frequency limits are only an outer screening envelope. None of these files supplies a verified throat chamber/adapter, diaphragm-only Mmd or a separately characterized rear acoustic load. They correctly produce `insufficient_evidence` even when the linear model meets the requested band.

The calculation checks sinusoidal peak excursion and real electrical input against supplied Xmax and nominal/AES power at the requested RMS voltage. These are rejection checks, not a prediction of broadband thermal capacity, distortion, or maximum SPL. Program-power marketing values are retained separately and never replace nominal power. Raw web pages were not archived; source URLs and the transcription date are provided explicitly.
