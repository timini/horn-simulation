# Worked example: an 800 to 1600 Hz horn

This walks through one auto-mode run from command to result. The result is an experimental prediction, not a qualified or buildable assembly.

## Run it

After the [README](../README.md#quick-start) prerequisites:

```sh
just build
just run-auto --target_f_low 800 --target_f_high 1600 --drivers_db data/drivers-curated \
  --lem_top_n 3 --refinement_budget 2 --num_bands 2 --num_intervals 101 --mesh_size 0.01
just latest-run
```

`just latest-run` prints the run folder. Open `outputs/auto/report/auto_report.html` inside it.

The settings keep the run short: three screened geometries go to FEM and refinement gets two extra solves. Drive is 2.83 V RMS, the observer is 1 m on axis, the mouth model is `modal_baffled` and the ripple limit is 6 dB, all defaults. A larger budget may find a different design. This run does not show a global optimum.

## What happened

The numbers below are from run `20261116T223748Z-45ab2da6`, at source revision `4307d67` (the pipeline code is the same as at `ba3cb10`). It took under two minutes on an Apple laptop. Run folders are not committed; rerun the command to reproduce them. Small differences can come from a different machine.

1. **Driver screen.** Five of the seven curated drivers passed: the three cone drivers and two Peerless compression drivers.
2. **Throats.** Five throat radii from 12.7 mm to 130.7 mm (the ka = 3.83 limit at 1600 Hz), including one that puts the 6NMB420's estimated front chamber corner at 1600 Hz.
3. **Grid and screen.** 230 geometries across the seven profiles, 1,150 driver pairs scored by the Webster screen. Three went to FEM: one `cd`, one exponential and one scaled tractrix.
4. **FEM and ranking.** 25 driver and horn pairs were scored. Both Peerless drivers were rejected: the band is outside their rated usable band, and one also exceeds Xmax at 2.83 V.
5. **Refinement.** Two more `cd` geometries, starting from the best FEM design.

## Result

| Item | Value |
| --- | --- |
| Driver | 18Sound 6NMB420 (8 ohm), curated record |
| Profile | `cd` (exponential throat into a cone) |
| Throat diameter | 72.12 mm |
| Mouth diameter | 177.42 mm |
| Length | 117.91 mm |
| Target-band ripple | 2.08 dB |
| Mean target-band level | 104.13 dB at 1 m, 2.83 V RMS |
| f3 | 485 Hz and 2024 Hz (3 dB below the band mean) |
| Score | 0.791 |
| Status | `insufficient_evidence` |

The next two designs (another `cd` and an exponential horn of the same throat and mouth) score within 0.02 and are marked near ties. Treat the three as equivalent.

The level includes the 6NMB420's estimated front chamber: 139 cm3, the air in its cone closed by a flat adapter plate, with a corner near 1600 Hz on this throat (`front_chamber_estimated: true`). Without the chamber the same horn gives 3.23 dB ripple and 102.33 dB (`without_front_chamber`). The real adapter will differ: a spacer adds volume, a phase plug or cone filler removes it. Put the real volume in the driver record as `front_chamber_volume_cm3` and rerun.

Evidence gaps on the winner: `driver_interface_unverified`, `moving_mass_and_rear_load_not_separated`, `front_chamber_volume_estimated` and `front_chamber_not_lumped_in_band`.

## Files

Inside `outputs/auto/`:

- `refinement/refine_cd_0001.step`: the winning air volume, true size in millimetres.
- `report/profile_01_refine_cd_0001.csv`: its wall profile, `z_mm` and `r_mm`.
- `report/coupled_01_18sound-6nmb420_refine_cd_0001.csv`: level, electrical impedance, excursion, power and group delay.
- `refinement/refine_cd_0001_results.csv`: the raw FEM output.
- `refinement/ranked_results.json`, `report/auto_rejections.json`, `refinement/search_audit.json`: ranking, rejections and the refinement trail.

[OUTPUTS](OUTPUTS.md) explains every file. The STEP describes the air passage only: no walls, mounting or adapter.

## Resume

```sh
python3 scripts/run_pipeline.py --run-dir results/<run-id> -resume
```

The source, inputs and images must be unchanged. See [USAGE](USAGE.md#resume).

## Earlier archived runs of this band

Two earlier runs of the same command are archived with full provenance. Both used an ideal area transformer at the throat, before the front chamber model, and chose a hyperbolic horn with a 65.0 mm throat, 177.4 mm mouth and 117.9 mm length:

- `flanged_piston` mouth: [archive](../data/validation/worked_example_800_1600.tar.gz) and [manifest](../data/validation/worked_example_800_1600_manifest.json). About 0.99 dB ripple and 104.94 dB at the time. Its [resolution study](CANDIDATE_RESOLUTION.md#executed-800-to-1600-hz-example) passes all six comparisons.
- `modal_baffled` mouth: [examples/800-1600-modal](../examples/800-1600-modal/README.md), 1.49 dB and 104.88 dB at the time. Its [resolution study](MODAL_CANDIDATE_RESOLUTION.md) passes too.

With the current front chamber estimate, the archived modal hyperbolic horn gives about 6.4 dB ripple, which is why the current run picks a larger throat. The archives are kept unchanged as evidence of what the code did then.

The [complete-horn radiation comparison](EXTERIOR_RADIATION_VALIDATION.md) and [STATUS](STATUS.md#remaining-work) explain what still separates this from a physically validated recommendation.
