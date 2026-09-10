# 800–1600 Hz: generate an experimental horn candidate

Open [report.html](report.html) locally in a browser to inspect the predicted response, shortlisted designs, assumptions and evidence gaps. The report is self-contained; Docker is only required to generate a fresh design.

| Selected numerical candidate | Value |
| --- | --- |
| Driver | 18Sound 6NMB420 |
| Horn profile | Hyperbolic |
| Throat diameter | 64.996 mm |
| Mouth diameter | 177.418 mm |
| Acoustic length | 117.906 mm |
| Target-band ripple | 1.485 dB |
| Mean target-band prediction | 104.882 dB SPL, 2.83 V RMS, 1 m on axis |
| Radiation setting | Nonuniform axisymmetric mouth in an infinite baffle (`modal_baffled`) |
| Physical evidence status | Insufficient evidence; experimental prediction |

The selected candidate is near-tied with other designs. Its throat dimensions describe an ideal acoustic connection. The real cone driver's chamber/adapter, rear load and moving-mass decomposition have not been qualified.

## Inspect the design

- [horn.step](horn.step): winning acoustic air volume, readable in a STEP-capable CAD viewer. Walls, mounting and the physical driver adapter require separate mechanical design.
- [response.csv](response.csv): predicted coupled response for the selected driver and horn.
- [horn-impedance.csv](horn-impedance.csv): raw horn FEM results, including throat acoustic impedance and modal mouth coefficients.
- [ranking.json](ranking.json): the full exported shortlist and evidence flags.
- [resolved-specification.json](resolved-specification.json): all resolved settings from the original run, including its historical output path.
- [original-run-manifest.json](original-run-manifest.json): original source, input, container and completion identities.
- [example-manifest.json](example-manifest.json): maps these conveniently named files to their original archive members and hashes.

The files are copied unchanged from the [completed example archive](../../data/validation/worked_example_modal_800_1600.tar.gz), executed at source `77d915b`. This folder contains selected outputs, while the archive preserves the complete output collection. The original run manifest describes that original collection and its paths.

## Generate a fresh design

From the repository root, after installing the [README prerequisites](../../README.md):

```sh
just build
just run-auto --target_f_low 800 --target_f_high 1600 \
  --drivers_db data/drivers-curated --lem_top_n 3 --refinement_budget 2 \
  --num_bands 2 --num_intervals 101 --mesh_size 0.01 \
  --radiation_model modal_baffled
just latest-run
```

The last command prints the new run directory. Open `outputs/auto/report/auto_report.html` inside it. Change `target_f_low` and `target_f_high` for a different band; optional `--max_mouth_radius` and `--max_length` limits are in metres. The search can return no feasible design or insufficient evidence. New requests do not inherit qualification from this example.

The sample uses a small search budget: three screened geometries and two refinement proposals. It demonstrates a repeatable workflow, not a global optimum. The [candidate mesh/loft/frequency study](../../docs/MODAL_CANDIDATE_RESOLUTION.md) and [four-band search comparison](../../docs/MODAL_SEARCH_VALIDATION.md) document its numerical evidence. [Issue #81](https://github.com/timini/horn-simulation/issues/81) tracks the remaining complete-assembly evidence.
