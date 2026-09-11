# Driver sources for a compact three-way coaxial horn

**Research date: 11 September 2026.** Scope: nominal 6.5-inch cone drivers for an experimental 500–6,500 Hz mid horn, with a cylindrical HF unit in its centre and the mid assembly inside a 15-inch LF horn. Source discovery includes adjacent 6-inch products, but the automated exact-size run admits only verified 6.5-inch records. Prices, UK stock, lead times and sample-to-sample measurements are outside this assessment.

## Findings

The catalogue needed both broader coverage and a source-integrity repair. The initial 1,867 records included three large groups totalling 376 entries with repeated motor parameters across unrelated models. A numerical plausibility check could pass these entries because the substituted parameter sets described plausible *different* drivers. On 11 September, the shared secondary source also returned the same B&C 12FL64 page for unrelated requested models. Its existing origin-identity guard correctly refused a refresh; the historic records predated that protection. This finding comes from the repository audit and direct source responses, rather than a manufacturer claim. The source should remain disabled for bulk refresh until routing is healthy.[^1]

The best new **data-access source** is REDCATT's documented engineering API. The best clearly recent **product lead** is SB Audience NERO-6MRN151D: the retrieved datasheet bears revision R.0 / 03.08.2026. These are separate kinds of novelty. A model newly added to this project's catalogue is not necessarily newly released. PHL's substantial missing range, for example, is useful coverage rather than evidence of a 2026 launch.[^2][^3][^4]

For the requested band, the primary-source shortlist includes Beyma 6MCF200Nd, Beyma 6MI100/6MI90, B&C 6PEV13 and Celestion CF0617M. Published frequency ranges justify investigating them; they do not establish smooth horn operation to 6.5 kHz. The automated study selects the **Beyma 6MCF200Nd within its simplified model**. REDCATT 6NPM deserves further investigation, but its imported formal engineering fields do not establish a usable upper band. SB Audience's recent driver is retained in the shared catalogue as a 6-inch product, excluded from the exact 6.5-inch screen.[^3][^5][^6][^7][^8][^9]

## Sources worth adding and maintaining

| Source | What is available | Assessment for this project |
|---|---|---|
| REDCATT | Structured public JSON, unit/lifecycle contract, drawings, STEP and response assets | Strong automated source. Three 6.5-inch products imported; selected-size refresh command added. |
| PHL Audio | Product index and consistent two-page manufacturer PDFs | Seventeen 6.5-inch references found; fifteen single-cone variants imported. Two 1365 coaxials excluded from single-cone ingestion. Most relevant mids specify 5 kHz upper operation. |
| SB Audience | Product pages and downloadable specification sheets | NERO-6MRN151D is a recent, useful adjacent-size lead. Preserve the manufacturer's 6-inch designation. |
| Beyma | Detailed PDF sheets and catalogue | 6MCF200Nd adds a sealed-back, high-resonance mid candidate. 6MI100 and 6MI90 refreshed from primary specifications. |
| Celestion | Structured product specifications and published test-bench evidence | CF0617M has attractive phase-plug clearance features, but inductance evidence needs care. |
| B&C | Model/impedance-specific pages and PDFs | Correct 6MDN44 impedance variants explicitly; retain 6PEV13 as a wide-band comparison with lifecycle/availability still to check. |
| LaVoce | Current product pages, PDF catalogue, comparative XLSX | Current pages take precedence over older catalogue values. MAF061.50 and MAN061.80 refreshed; both specify an upper limit below 6.5 kHz. |
| PRV Audio | Detailed technical tables, drawings, PDFs and crossover guidance | Two neo mids imported with explicit size conflicts. Marketing power and headline bandwidth need careful interpretation. |
| Ground Zero | Product pages and owner manuals | GZCM 6.5N-PROX is a useful missing brand/model lead. The manual examined lacks explicit Sd; no complete simulation record fabricated. |
| Dayton Audio | Manufacturer product family pages; supplier-hosted manufacturer specifications and FRD/ZMA downloads | ODEUM 6.5N merits a dedicated primary-file refresh. Not imported from the currently broken secondary database. |
| Oberton | Detailed primary product pages and drawings | 6NMB200V is explicitly 6.5 inches; examined 16-ohm data specifies 200–5,000 Hz, so not a nominal full-band candidate. |
| BMS / Kartesian | Manufacturer specifications | Useful wider source coverage; examined 6N160 and Wom165_vHE-e2 are lower-band designs, not evidence for the requested upper crossover. |

Sources for this comparison are the linked manufacturer records and source inventory below.[^2][^3][^4][^5][^6][^7][^8][^9][^10][^11][^12][^13][^14][^15][^16][^17][^18]

The landscape also includes the catalogue's established 18Sound, Faital, SICA, Ciare, Eminence, JBL, RCF and Precision Devices sources. This refresh does **not** claim their complete ranges were independently rechecked. In particular, the mass Faital and 18Sound anomalies remain quarantined. New-source research should expand discovery without turning unreviewed records into apparently verified simulation inputs.

## Product evidence and what it permits

### Exact 6.5-inch candidates with published upper coverage

| Driver, nominal 8 ohms | Published range, Hz | Fs, Hz | Le used, mH | Maximum frame span × depth, mm | Key qualification |
|---|---:|---:|---:|---:|---|
| Beyma 6MCF200Nd | 400–12,000 | 406 | 0.10 | 174 × 75 | Factory-sealed back; published compliance includes that load. |
| Beyma 6MI100 | 100–8,000 | 128 | 0.20 | 174 × 85 | Open basket; horn adapter and rear volume needed. |
| Beyma 6MI90 | 140–8,000 | 134 | 0.40 | 174 × 84 | Open basket; upper nominal range alone does not qualify distortion. |
| B&C 6PEV13 | 150–8,000 | 126 | 0.60 | 187 × 78 | Low moving mass; primary sheet found on legacy manufacturer host. Check production and supply. |
| Celestion CF0617M | 300–7,000 | 116.6 | 1.73 | 189 × 78.5 | Published Le materially changes simplified high-frequency prediction. |

The table uses manufacturer data, with exact URLs attached to each shared catalogue record. Size refers to nominal product class; frame span includes the much larger mounting envelope where known.[^5][^6][^7][^8][^9]

The Beyma 6MCF200Nd deserves attention because its published band reaches well past the desired crossover and its motor inductance is low. Its Fs is high because this is a specialised sealed-back mid. Treating it as a normal open driver and adding a second modelled rear chamber changes its compliance incorrectly. Its quoted excursion also uses a manufacturer definition; it should not be translated directly into a guaranteed clean horn SPL.[^5]

Celestion remains interesting for the physical phase-plug concept. Its inverted dustcap helps close placement, but the current product page's 1.73 mH differs substantially from the approximately 0.28–0.30 mH figures in Vance Dickason's earlier inductance-versus-displacement test. Those measurements have different contexts. The study retains the current nominal value and labels 0.29 mH as a sensitivity case, rather than replacing a source value with a more favourable one. A measured complex impedance curve would be a better input than either constant.[^9][^19]

### New and missing sources that need qualification

**REDCATT 6NPM.** The API identifies a 6.5-inch, 8-ohm mid with Fs 150 Hz, Mms 11.1 g, Le 0.16 mH and 180 W AES. Its numerical model is promising, but the API's free-text summary is not treated as measured bandwidth. The listed 162 mm overall diameter is smaller than its 172 mm bolt circle, indicating that it cannot be accepted as the maximum frame envelope. The record therefore leaves the maximum width unknown until the drawing is resolved. The 0.4 mm published Xmax is another reason to verify output requirements with a sample.[^20]

**SB Audience NERO-6MRN151D.** The revision dated 3 August 2026 specifies 230–10,500 Hz, 180 W AES, Fs 230 Hz and Le 0.17 mH. This is the clearest recent product found in this research. Its specification calls it 6 inches, with a 185 mm overall frame. It is included in the shared catalogue with those original labels. The manufacturer's site still lists NERO-6MRN150D as well, so an explicit production/supersession relationship should be confirmed before treating the old model as discontinued. The larger published bandwidth does not remove the need for cone-contour and horn-response measurements.[^3]

**PHL.** The manufacturer index offers seventeen 6.5-inch references. The fifteen imported single-cone variants span bass-mid and dedicated-mid roles. For example, 1120, 1663 and 1660NdM specify 300–5,000 Hz; 1660NdM-SQ2 specifies 450–5,000 Hz and includes its own chamber arrangement. These are credible sources for a lower crossover study, but none of the imported PHL records claims coverage to 6.5 kHz. They remain comparisons, not nominal winners. Their separate Le values at 1 and 10 kHz also illustrate why a single inductance number is an imperfect high-frequency model.[^4][^21][^22][^23]

**PRV.** 6MR500-NDY and 6MR600X-NDY advertise 6.5-inch products, yet both technical tables say 6 inches / 152 mm. Their 164.5 mm outer frames do not resolve the nominal-class contradiction automatically. The shared records retain the conflict and cannot satisfy exact-size filtering. Their continuous-program headline powers are twice their nominal ratings: the imported nominal values are 250 and 300 W. PRV recommends low-pass points of 5 and 6 kHz respectively when pairing with HF drivers. Their broad advertised response limits should therefore not be read as unconditional crossover recommendations.[^11][^12]

**Ground Zero.** GZCM 6.5N-PROX is explicitly 6.5 inches, 4 ohms and 180 Hz–12 kHz in the manufacturer manual. Its supplied central aluminium phase plug is part of the mid driver, not an independent tweeter. Replacing or surrounding that feature with a bullet HF driver is a mechanical/acoustic redesign. The manual includes motor parameters but no explicit Sd in the inspected specification table, so it remains a discovery lead rather than receiving a guessed piston area. Its Xmax is stated peak-to-peak and must not be stored as one-way excursion without conversion.[^13][^24]

**LaVoce.** Current MAF061.50 specifications differ from the 2021 catalogue in Re, Fs, Sd and other fields. For example, the current page gives Sd 132 cm², while the older sheet gives 143.1 cm². Mixing a new motor parameter with the old area would create a parameter set that belongs to neither revision. The current-page set is retained together. The comparative XLSX is useful for discovery, but its filename says 2023 Rev00 even though its hosting path contains 2026; upload location is not a product revision date.[^10][^25][^26]

## Catalogue repair and source policy

The shared catalogue now contains **1,889 records: 29 manufacturer-verified parameter records, 375 quarantined records and 1,485 legacy-unverified records**. Of the 29 refreshed/added records, 26 are verified nominal 6.5-inch cone drivers, one is a verified 6-inch SB Audience driver, and two PRV records have unresolved nominal-size conflicts. Twenty-two records were added and seven existing records replaced with sourced sets. One of the original 376 anomalous records, B&C 6MDN44 8 ohms, has been recovered from its correct primary page. These counts are reproducible in [the audit](../../data/catalogue-audit.json).

“Manufacturer verified” here means the stored facts were checked against the identified manufacturer's publication. It does not mean independently measured, available from stock, or validated in this horn. The remaining legacy records are deliberately labelled unverified. Quarantine is retained on disk for traceability and rejected by the simulation loader. Legacy nominal-size estimates cannot satisfy size-constrained searches.

Source records should preserve manufacturer, exact model, impedance variant, revision or retrieval date, units, power definition, and whether the driver includes a rear chamber. Unknown Le must remain unknown rather than becoming zero. Dry Mmd and air-loaded Mms are different quantities. Nominal size, piston area, outer frame, bolt circle and horn exit area are also different dimensions; none is a safe universal replacement for another.

The changes address concrete failure modes: repeated-parameter detection and quarantine; strict nominal-size filtering; incorrect recommendation-card extraction; preservation of stale unsourced numeric fields; Unicode page parsing and URL equivalence; invalid/nonfinite numbers; and shell wrappers that previously lost a failing command's exit status. The source audit is independent of numerical plausibility so a physically possible wrong product cannot receive a clean source verdict.

REDCATT's API contract is especially useful because it defines nulls, mixed LF/HF ratings, stable ordering codes and lifecycle behaviour. Its full-list disappearance rule cannot be safely applied to a filtered or delta response. The importer therefore converts only a validated selected-size response and performs no inferred deletions. The API permits attributed research/comparison use but prohibits bulk republication as a competing catalogue; this repository imports the three relevant engineering records rather than reproducing the full product marketing database.[^2][^27]

## Implications for the coaxial horn

For an unobstructed area equivalent to a circle of diameter D and a central cylinder of diameter B, the annular outer diameter is **sqrt(D² + B²)**. A 64 mm central housing and a 76 mm equivalent-area throat therefore require a 99.36 mm outer opening. Calling that throat “76 mm” without identifying it as an *area-equivalent diameter* would conceal the packaging problem.

The completed reduced-order study evaluated 480 horn geometries for each of 26 verified 6.5-inch drivers, plus one separately labelled Celestion inductance scenario: **12,960 calculations**. Profiles, dimensions, inputs, response curves, numerical checks and acoustic STEP files are preserved in [the study report](../../examples/6p5-mid-500-6500/README.md).

The nominal winner is Beyma 6MCF200Nd with a conical 190 mm mouth, 125 mm acoustic length and 76 mm area-equivalent throat around the 64 mm housing. The model predicts 3.55 dB peak-to-peak variation and 104.44 dB mean level at 2.83 V / 1 m across 500–6,500 Hz. These are model outputs, not manufacturer measurements. The chosen throat lies at the largest tested equivalent diameter, so this is a finite-grid result rather than proof of a global optimum or of the smallest possible throat.

REDCATT 6NPM gives a competitive exploratory result but cannot be a nominal winner while its formal usable band remains unknown. Celestion's ranking changes markedly with assumed inductance. These observations favour measuring shortlisted samples before treating a computed winner as a purchase or build decision.

The model represents a rigid-piston motor and a plane-wave acoustic network. It includes the central body's blocked area but does not resolve cone breakup, equal-length phase-plug passages, annular higher modes, HF diffraction, distortion, or the response of the surrounding 15-inch horn. Its blunt central-body termination is a packaging placeholder. The production modal FEM path currently rejects annular inlet ports; these results are explicitly reduced-order, not full FEM.

The next physical qualification requires the actual cone/dustcap profile, a feasible front chamber and equalising passages, mounting ears and rear packaging, then horn-loaded response/impedance and crossover measurements. The 15-inch outer horn dimensions remain necessary to assess the assembled fit. A T90A-sized housing is reserved, but the T90A's stated minimum crossover is 7 kHz, so that tweeter itself is not qualified for the 6.5 kHz target.[^28]

## Sources

All undated web pages were retrieved or checked on 11 September 2026. Product launch dates are not inferred from crawl dates. PDF source hashes are stored in the PHL records; current facts and source URLs are stored with every imported record.

[^1]: Loudspeaker Database, [B&C 6MDN44 8-ohm route](https://loudspeakerdatabase.com/BC/6MDN44_8%CE%A9), inspected 11 September 2026; wrong-page response observed locally. Repository audit: [catalogue-audit.json](../../data/catalogue-audit.json).
[^2]: REDCATT, [Engineering Data](https://www.redcatt.net/engineering-data), current engineering-data and API overview.
[^3]: SB Audience, [NERO-6MRN151D product page](https://www.sbaudience.com/index.php/products/woofers/nero-6mrn151d/) and [manufacturer datasheet](https://www.sbaudience.com/index.php/download_file/-/view/3645/), R.0 / 03.08.2026, one page. PDF retrieved directly from its download endpoint.
[^4]: PHL Audio, [Product index](https://phlaudio.com/products/index.html), 6.5-inch references.
[^5]: Beyma, [6MCF200Nd datasheet](https://www.beyma.com/speakers/Fichas_Tecnicas/beyma-speakers-data-sheet-low-mid-frequency-6MCF200Nd.pdf), pp. 1–2.
[^6]: Beyma, [6MI100 product specifications](https://www.beyma.com/en/products/c/low-mid-frequency/106MI108/loudspeaker-6mi100-8-oh/) and [datasheet](https://www.beyma.com/speakers/Fichas_Tecnicas/beyma-speakers-data-sheet-low-mid-frequency-6MI100.pdf).
[^7]: Beyma, [6MI90 datasheet](https://www.beyma.com/speakers/Fichas_Tecnicas/beyma-speakers-data-sheet-low-mid-frequency-6MI90.pdf).
[^8]: B&C Speakers, [6PEV13 8-ohm datasheet](https://old.bcspeakers.com/en/products/lf-driver/6-5/8/6pev13.pdf) and [6MDN44 8-ohm specifications](https://www.bcspeakers.com/en/products/lf-driver/6.5/8/6MDN44).
[^9]: Celestion, [CF0617M](https://celestion.com/product/cf0617m/), current product specifications.
[^10]: LaVoce, [MAF061.50](https://lavocespeakers.com/product/maf061-50/) and [MAN061.80](https://lavocespeakers.com/product/man061-80/), current specifications.
[^11]: PRV Audio, [6MR500-NDY](https://prvaudio.com/products/prv-6-5-inch-neo-midrange-speaker-6mr500-ndy/), technical specifications and crossover guidance.
[^12]: PRV Audio, [6MR600X-NDY](https://prvaudio.com/products/6mr600x-ndy-6-5-inch-neodymium-loudspeaker/), technical specifications and crossover guidance.
[^13]: Ground Zero, [GZCM 6.5N-PROX](https://www.ground-zero-audio.com/en/gzcm-6-5n-prox/).
[^14]: Dayton Audio, [ODEUM family](https://www.daytonaudio.com/topic/Odeum) and [Pro Series index](https://www.daytonaudio.com/category/283/pro-series); primary manufacturer discovery sources. Supplier resource lead: [Parts Express 295-638](https://www.parts-express.com/Dayton-Audio-6MB200N-8-6.5-Professional-Neodymium-Mid-Bass-295-638).
[^15]: Oberton, [6NMB200V specifications](https://www.oberton.com/en/products/neodymium-loudspeakers/458-6nmb200v.html?showall=1), examined table is the 16-ohm variant.
[^16]: BMS, [6N160 specifications](https://www.bmsspeakers.com/index.php-226.html?id=6n160_specification).
[^17]: Kartesian, [Wom165_vHE](https://www.kartesian-acoustic.com/product/wom165_vhe).
[^18]: LaVoce, [6.5-inch LF index](https://lavocespeakers.com/product/?size=6.5-6.5&type=LF+Products), current discovery inventory.
[^19]: Vance Dickason, Voice Coil/audioXpress, [Test Bench: Celestion CF0617M ProSound Midrange Driver](https://audioxpress.com/article/test-bench-celestion-cf0617m-prosound-midrange-driver), original independent bench measurements, 2015.
[^20]: REDCATT, [6NPM](https://www.redcatt.net/products/dr-65-006-8r-b2-c0001), ordering code DR-6.5-006-8R-B2-C0001, API updated_at 4 September 2026.
[^21]: PHL Audio, [1120 datasheet](https://phlaudio.com/fileadmin/user_upload/phl_audio/1120_SpecSheet.pdf), pp. 1–2.
[^22]: PHL Audio, [1663 datasheet](https://phlaudio.com/fileadmin/user_upload/phl_audio/1663_SpecSheet.pdf), pp. 1–2.
[^23]: PHL Audio, [1660NdM](https://phlaudio.com/fileadmin/user_upload/phl_audio/1660NdM_SpecSheet.pdf) and [1660NdM-SQ2](https://phlaudio.com/fileadmin/user_upload/phl_audio/1660NdM-SQ2_SpecSheet.pdf), pp. 1–2. Remaining imported variant sources are linked directly in their [catalogue records](../../data/drivers/PHL/).
[^24]: Ground Zero, [GZCM N-PROX owner manual](https://www.ground-zero-audio.com/wp-content/uploads/GZCM-N-PROX_OM_AWV2.1.pdf), AWV2.1, p. 3.
[^25]: LaVoce, [Products Catalogue 2021](https://www.lavocespeakers.com/wp-content/uploads/2021/06/LAVOCE_Products-Catalogue_2021-low.pdf), MAF061.50 sheet, p. 44.
[^26]: LaVoce, [Product Comparison Table 2023 Rev00](https://lavocespeakers.com/wp-content/uploads/2026/01/LAVOCE-PRODUCT-COMPARISON-TABLE-2023-Rev00.xlsx), discovery source; filename revision differs from hosting-path year.
[^27]: REDCATT, [API schema](https://www.redcatt.net/api/v1/products/schema), v1.1, unit definitions, lifecycle and reuse terms. [Product feed](https://www.redcatt.net/api/v1/products?per_page=500).
[^28]: Fostex, [T90A](https://www.fostex.jp/products/t90a/), manufacturer specifications; recommended crossover above 7 kHz with at least 12 dB/octave filtering.
