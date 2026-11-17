# Outputs

Every launcher run writes to `results/<run-id>/`. Find the newest completed one with `just latest-run`.

```text
results/<run-id>/
  manifest.json          provenance: source revision and hashes, inputs, images, Nextflow, command, status
  source.tar.gz, source.patch
  containers.config, nextflow.log, work/
  outputs/
    resolved_specification.json   every parameter value actually used, including defaults
    ...                           single mode files, or auto/ for auto mode
```

All levels are dB re 20 uPa. All complex quantities are RMS phasors with time dependence exp(+iwt). Lengths in CSV and JSON files are metres unless the column name says `_mm`.

## Single mode

| File | Contents |
| --- | --- |
| `horn.step` | Generated acoustic air volume (absent when `step_file` is given). |
| `final_results.csv` | Solver output for the whole sweep, one row per frequency. Columns below. |
| `kpis.json` | Peak-lobe metrics of the horn-only response: peak level and frequency, f3 points, bandwidth, ripple and mean level, with `*_is_bound` flags. |
| `frequency_response.png`, `impedance.png`, `phase_response.png`, `dashboard.png` | Plots of the horn-only solve. |
| `horn_3d.png` | Render of the profile, or the image given with `horn_3d_png`. |
| `single_report.html` | Self-contained report with all of the above. |
| `coupled_spl.csv` | With `driver_id` only. `frequency`, `coupled_spl` (observer level with the driver), `horn_only_spl` (solver mouth-plane level) and the electrical columns listed under auto mode. |
| `coupled_spl.png`, `driver_horn_kpis.json` | With `driver_id` only. Coupled plot, and the same assessment fields as an auto ranking row, evaluated over `target_f_low` to `target_f_high` clipped to the sweep. |

### Solver CSV columns

`final_results.csv` and every `*_results.csv` in auto mode have these columns:

| Column | Meaning |
| --- | --- |
| `frequency` | Hz. |
| `spl` | Spatial RMS pressure level over the mouth plane for 1 Pa at the throat. This is a transfer quantity, not the level at a listener. |
| `z_real`, `z_imag` | Specific acoustic impedance at the throat, Pa s/m. |
| `inlet_p_*`, `inlet_u_*` | Complex throat pressure and volume velocity. |
| `mouth_p_*`, `mouth_u_*` | Complex area-mean mouth pressure and outward volume velocity. |
| `modal_v_<n>_real/imag`, `modal_mode_count` | The 16 mouth velocity mode amplitudes (`modal_baffled` only). Used for the on-axis level and directivity. |
| `inlet_area_m2`, `mouth_area_m2`, `mesh_*_area_m2` | CAD port areas and the meshed areas. |
| `input_acoustic_power_w`, `mouth_acoustic_power_w`, `*_wall_power_w`, `power_balance_available` | Power accounting. Wall terms are nonzero only with `boundary_layer`. |
| `radiation_model`, `loss_model`, `element_degree`, `flange_width_m`, `bc_mode` | Model settings. |
| `mesh_nodes`, `mesh_cells`, `mesh_size_m`, `mesh_elements_per_wavelength`, `mesh_coarsened_for_node_limit`, `mesh_throat_*` | Mesh used for that band. |
| `relative_residual`, `ksp_type`, `pc_type`, `converged_reason` | Linear solve health. |
| `air_*`, `schema_version`, `phasor_convention` | Air properties and conventions. |

## Auto mode

Everything is under `outputs/auto/`.

| File | Contents |
| --- | --- |
| `prescreen_result.json` | Drivers that passed the screen, the throat radii sampled and their flags (for example `throat_above_first_asymmetric_mode`, `diaphragm_smaller_than_throat`). |
| `design.json` | Derived search ranges for throat, mouth and length, where each bound came from (`bound_sources`: user, heuristic, driver prescreen or size cap), the solver sweep (`sim_freq_range`) and any size compromise in octaves. |
| `candidates.csv` | The full geometry grid: `candidate_id`, `profile`, `throat_radius`, `mouth_radius`, `length`. |
| `lem_results.json`, `lem_filtered_candidates.csv` | Webster screen of every grid geometry with every driver, the shortlist policy and the geometries sent to FEM. |
| `geometry/horn_<id>.step` | STEP for each shortlisted geometry. |
| `<id>_results.csv` | Solver output per shortlisted geometry (columns above). `<id>_rejected.json` instead when meshing, CAD or band boundary checks rejected it. |
| `ranked_results.json` | Ranking after the FEM shortlist, before refinement. |
| `refinement/ranked_results.json` | Final ranking after refinement. |
| `refinement/refine_<profile>_<n>.step`, `refinement/*_results.csv` | Geometry and solver output for each refinement evaluation. |
| `refinement/search_audit.json` | Refinement seed and its source, every evaluation, budget use, bounds reached and widened. |
| `report/auto_report.html` | Self-contained report. Start here. |
| `report/auto_summary.txt` | Plain-text summary. |
| `report/auto_ranking.json` | Feasible rows, sorted by score, with near-tie notes. |
| `report/auto_comparison.png` | Coupled response of the top candidates. |
| `report/coupled_<rank>_<driver>_<horn>.csv` | Per ranked candidate: response and electrical columns, below. |
| `report/auto_electrical.png` | Electrical impedance, excursion against Xmax and group delay for the ranked candidates. |
| `report/directivity_<horn>.csv`, `report/auto_beamwidth.png` | Experimental directivity, `modal_baffled` horns only. |
| `report/profile_<rank>_<horn>.csv` | Wall profile for builders: `z_mm`, `r_mm`, from the same profile function as the STEP. |
| `report/auto_rejections.json` | Rejection counts by reason, with margins. |

When nothing can be ranked (no driver passes the screen, the size limits leave no geometry, or no geometry is inside the radiation model's domain), the run still completes and writes `ranked_results.json` with `status: no_feasible_design`, a `reason`, and a report that says why.

### Ranking rows

Each row in `ranked_results.json`, `refinement/ranked_results.json` and `auto_ranking.json` is one driver on one horn. The main fields:

| Field | Meaning |
| --- | --- |
| `driver_id`, `horn_label`, `profile`, `throat_radius`, `mouth_radius`, `length` | The pair. |
| `composite_score`, `score_terms`, `score_gap_from_best` | Score and its four terms. See [SEARCH_AND_SCORING](SEARCH_AND_SCORING.md). |
| `avg_sensitivity_db` | Log-frequency mean observer level over the target band at `voltage_rms` and `observation_distance`. |
| `passband_ripple_db` | Peak-to-peak level over the target band. |
| `target_min_level_db`, `target_max_level_db` | Lowest and highest level in the band. |
| `f3_low_hz`, `f3_high_hz`, `f3_*_is_bound`, `extension_octaves` | f3 points 3 dB below the band mean, found over the full sweep. A `_is_bound` flag means no crossing was found and the value is the sweep limit. |
| `model_feasible`, `rejection_reasons` | Whether every hard limit is met, and which are not. |
| `eligibility_status` | `model_feasible`, `insufficient_evidence` (feasible but with evidence gaps) or `infeasible`. |
| `evidence_gaps` | Missing driver or interface data, for example `driver_interface_unverified` or `front_chamber_volume_estimated`. |
| `front_chamber_estimated`, `estimated_chamber_reasons`, `without_front_chamber` | Front chamber policy; see [DRIVERS](DRIVERS.md#front-and-rear-chambers). |
| `max_displacement_peak_m`, `max_input_power_w` | Peak excursion and input power over the band at the drive voltage. |
| `test_power_*` | With `test_power_w`: the test voltage, excursion and Xmax margin at that power over the band and the whole sweep, and the largest in-band input power with its margin to the driver rating. |
| `radiation_model`, `output_metric`, `loss_model`, `element_degree` | Model used. `output_metric` is `modal_baffled_on_axis` or `uniform_baffled_piston_on_axis`. |
| `kpi` | Older peak-lobe metrics, as in single-mode `kpis.json`. |
| `validation_status` | Always `experimental_prediction`. |

### Coupled CSV columns

| Column | Meaning |
| --- | --- |
| `frequency`, `spl` | Hz, and the on-axis observer level at `observation_distance` with the driver at `voltage_rms`. |
| `ze_real_ohm`, `ze_imag_ohm`, `ze_mag_ohm` | Electrical input impedance. |
| `excursion_peak_mm` | Peak diaphragm excursion at the drive voltage. |
| `input_power_w`, `acoustic_power_w` | Real electrical input power and acoustic power into the horn. |
| `group_delay_ms`, `time_of_flight_removed_ms` | On-axis system group delay, with the straight-line delay (horn length plus observation distance, over c) removed. |
| `drive_voltage_rms` | The drive voltage. |
| `test_voltage_rms`, `test_excursion_peak_mm`, `test_xmax_margin_mm`, `test_input_power_w` | With `test_power_w` only. The model is linear, so these scale with voltage. |

### Directivity CSV

`frequency`, `theta_deg` (0 to 90 degrees in 5 degree steps), `spl_db` relative to the axis, `model` (`modal_baffled_far_field`) and `status` (`experimental, infinite baffle`). It assumes an infinite rigid baffle and the far field. It has not been compared with a measured polar. `auto_beamwidth.png` plots the -6 dB beamwidth from it.

### Rejections file

`auto_rejections.json` has the stage (`FEM` or `analytical prescreen`), the number of rejected pairs and one entry per reason, most frequent first, with a description, the count and `closest_excess` and `largest_excess` in `unit` where a margin is defined. For example a ripple rejection gives dB above `max_ripple_db`; an excursion rejection gives mm above Xmax.

Rejection reasons: `target_band_ripple_exceeded`, `minimum_output_not_met`, `maximum_length_exceeded`, `maximum_mouth_radius_exceeded`, `compression_ratio_exceeded`, `driver_larger_than_mouth`, `outside_driver_usable_band`, `driver_excursion_exceeded`, `driver_nominal_power_exceeded`, `test_power_excursion_exceeded`, `test_power_nominal_power_exceeded`, `radiation_model_out_of_domain`, `profile_geometry_infeasible`, `loft_accuracy_unmet`, `cad_geometry_rejected`, the solver's `mesh_node_limit_exceeded` and `throat_area_guard_failed`, and the band merge's `band_boundary_mismatch`. The last means two neighbouring solver bands, solved on separate meshes, disagree at their shared frequency by more than 0.5 dB in SPL or more than 5 % (plus 1 % of a reference scale) in a complex transfer such as mouth pressure. The candidate's `*_rejected.json` gives the worst quantity, frequency, jump and how many times the allowed difference it is.

## STEP files

Every generated STEP is the horn's air volume: the profile revolved about the z axis, inlet at z = 0, mouth at z = length. Files are written at true size, in millimetre coordinates under a millimetre unit declaration, so they open at the right size in CAD. They have no wall thickness, mounting flange or driver adapter. Treat them as the inside surface of a horn, not as a part to print.

## Reading f3 and ripple

Ripple and feasibility use every solver sample inside the target band plus the two band edges, with no smoothing, so a narrow notch cannot hide.

In rankings, f3 is 3 dB below the target-band mean. Each f3 is found from its band edge: outward over the full sweep if the level at the edge is above the threshold, inward otherwise. Crossings are interpolated in log frequency. The auto sweep runs one octave below `target_f_low` and half an octave above `target_f_high`. When no crossing is found, the value is the sweep limit, the `_is_bound` flag is set and the report prints it with `<=` or `>=`. For example `f3 high >= 2263 Hz` does not mean the horn cuts off at 2263 Hz.

The `kpi` block and single-mode `kpis.json` use a different reference: 3 dB below the largest sampled peak, on the contiguous lobe around that peak. Their ripple and mean level describe that lobe and can differ from the target-band values used for ranking.
