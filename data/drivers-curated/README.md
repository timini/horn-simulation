# Manufacturer-sourced numerical comparison set

The three cone drivers below were transcribed from the manufacturers' pages on 7 September 2026. Every JSON records its source, units, supplied fields and unresolved assumptions. The original large scraped database remains available; this set can be selected with `--drivers_db data/drivers-curated`.

| Driver | Variant | Primary source |
| --- | --- | --- |
| 18Sound 6NMB420 | 8 ohm | https://www.eighteensound.it/es/products/lf-driver/6-5/8/6NMB420 |
| Celestion AF3010 | 16 ohm | https://celestion.com/product/af3010-2/ |
| Celestion TF0512HE | 8 ohm | https://celestion.com/product/tf0512he/ |

Four compression drivers were added on 15 November 2026 with `horn-add-driver` from Tymphany "Driver Specification Sheets" (watermarked preliminary) hosted by distributors. They are the only compression drivers found with a complete published Thiele-Small set. The exit area is the nominal 1 inch exit; the usable band is the sheet's rated-noise-power test bandwidth, because the sheets give no other frequency range.

| Driver | Variant | Source |
| --- | --- | --- |
| Peerless DFM-2544R00-08 | 8 ohm, 1 in exit | https://www.parts-express.com/pedocs/specs/264-1466--peerless-dfm-2544r00-08-compression-dfm-driver-specifications.pdf |
| Peerless DFM-2535R00-08 | 8 ohm, 1 in exit | https://mm.digikey.com/Volume0/opasdata/d220001/medias/docus/7466/dfmspecsheet.pdf |
| Peerless DFM-2535R01-08 | 8 ohm, 1 in exit | https://mm.digikey.com/Volume0/opasdata/d220001/medias/docus/7466/DFM-2535R01-08.pdf |
| Peerless DFL-2525R00-08 | 8 ohm, 1 in exit | https://www.parts-express.com/pedocs/specs/264-1464--peerless-dfl-2525r00-08-titanium-dome-driver-specifications.pdf |

B&C, 18Sound, Celestion, BMS, RCF and Faital publish exit size, band and power for their compression drivers but not Fs, Bl, Mms and Sd together; those are recorded in `data/drivers-incomplete`.

This is a traceable parameter set, **not a physically validated horn-driver shortlist**. Manufacturer nominal frequency limits are only an outer screening envelope. None of these files supplies a verified throat chamber/adapter, diaphragm-only Mmd or a separately characterized rear acoustic load. They correctly produce `insufficient_evidence` even when the linear model meets the requested band.

The calculation checks sinusoidal peak excursion and real electrical input against supplied Xmax and nominal/AES power at the requested RMS voltage. These are rejection checks, not a prediction of broadband thermal capacity, distortion, or maximum SPL. Program-power marketing values are retained separately and never replace nominal power. Raw web pages were not archived; source URLs and the transcription date are provided explicitly.
