#!/usr/bin/env nextflow

// ========================================================================
// Mode: "single" (default) runs one horn profile; "auto" explores a grid
// of profiles and geometry, ranking driver-horn combinations.
// "fullauto" is an alias for "auto" with all geometry derived.
// ========================================================================
params.mode = "single"

// Horn Geometry — null means "derive from frequency" in auto mode,
// or use sensible defaults in single mode
params.throat_radius = null
params.mouth_radius = null
params.length = null
params.profile = "conical" // Horn flare profile: conical, exponential, hyperbolic, tractrix, os, scaled_tractrix, cd (lecleach: deprecated alias of scaled_tractrix)
params.num_sections = 20   // Minimum profile points in the revolved horn wall

// Simulation Settings
params.min_freq = 500      // Minimum frequency for the sweep in Hz
params.max_freq = 8000     // Maximum frequency for the sweep in Hz
params.num_intervals = 100 // Number of frequency steps in the sweep
params.min_points_per_octave = 24 // Lower limit on log-spaced solver points per octave (0 disables)
params.mesh_size = 0.01    // Largest element size in metres; finer of this and lambda/10 (P1) or lambda/6 (P2)
params.max_mesh_nodes = 100000 // Per-band node limit; P1 coarsens towards lambda/6 to fit, else the candidate is rejected

// Radiation model at the horn mouth. null chooses modal_baffled for generated
// (axisymmetric) geometry with lossless P1 elements, and flanged_piston for an
// imported STEP or other settings. modal_baffled passes the independent MMM
// infinite-baffle reference; flanged_piston misses it by 1.05 dB / 5.7 deg in
// throat impedance (docs/EXTERIOR_RADIATION_VALIDATION.md).
params.radiation_model = null
params.loss_model = "lossless"  // opt-in boundary_layer until horn-domain validation
params.flange_width = 0.0
params.minimum_wall_scale = null  // required for boundary_layer: minimum curvature/gap scale in metres
params.element_degree = 1

// Execution Settings
params.num_bands = 8       // Number of parallel jobs for the solver
params.outdir = null // direct invocation must choose an output directory

// Auto mode settings
params.target_f_low = 500
params.target_f_high = 4000
params.drivers_db = "data/drivers"
params.top_n = 10
params.profiles = null    // Optional comma list restricting auto profiles, e.g. "os,conical"

// Single mode optional driver coupling: when set, single mode also produces
// a driver+horn coupled SPL plot/CSV/KPI alongside the bare-horn outputs.
params.driver_id = null

// Single mode optional pre-built STEP file: when set, skip generate_geometry
// and use this user-supplied STEP (e.g. a horn-with-phase-plug compound).
params.step_file = null

// Single mode optional pre-rendered 3D geometry PNG: when set, skip
// render_horn_3d and use this image (so the report shows the actual STEP
// geometry, not a parametric reconstruction of the bare horn).
params.horn_3d_png = null
params.num_mouth_radii = 3 // Mouth radius grid points (when mouth_radius is null)
params.num_lengths = 3     // Length grid points (when length is null)
params.refinement_budget = 6
params.voltage_rms = 2.83
params.test_power_w = null   // Optional stress power (W into nominal impedance) for full-sweep excursion and Xmax margin
params.observation_distance = 1.0
params.max_ripple_db = 6.0
params.max_compression_ratio = 10.0
params.min_output_db = null  // Optional: minimum observer level anywhere in the target band (dB)
params.score_weights = null  // Optional: e.g. "sensitivity=0.625,ripple=0.25,extension=0.125,size=0" (normalised)
params.lem_top_n = 10       // Number of top candidates to pass from LEM prescreen to FEM
params.min_diameter = null  // Optional: minimum driver nominal diameter (inches)
params.max_diameter = null  // Optional: maximum driver nominal diameter (inches)
params.throat_fractions = null  // Optional CSV (e.g. "0.1,0.3,0.65,1.0") to override prescreen throat sampling
params.max_length = null    // Optional cap/extension of derived length range upper bound (m)
params.min_length = null    // Optional override of derived length range lower bound (m)
params.max_mouth_radius = null // Optional cap/extension of derived mouth radius upper bound (m)
params.min_mouth_radius = null // Optional override of derived mouth radius lower bound (m)

// Effective radiation model; see params.radiation_model.
def radiationModel() {
    if (params.radiation_model != null) return params.radiation_model
    def modalEligible = !params.step_file && params.loss_model == 'lossless' && params.element_degree == 1
    return modalEligible ? 'modal_baffled' : 'flanged_piston'
}

// Solver band edges and point counts come from horn_core.frequency_grid
// (equal log-width bands), so every process uses the same grid.
def bandGridShell(min_f, max_f, band_index) {
    return "band_grid=\$(python3 -m horn_core.frequency_grid --min-freq ${min_f} --max-freq ${max_f} --num-bands ${params.num_bands} --num-intervals ${params.num_intervals} --min-points-per-octave ${params.min_points_per_octave} --band ${band_index})\n" +
        "    read -r band_low band_high band_points <<< \"\$band_grid\""
}

def gridFlags() {
    return "--num-intervals ${params.num_intervals} --min-points-per-octave ${params.min_points_per_octave}"
}

// Score weights and output floor for every TargetSpec built in Python.
def scoringTargetArgs() {
    def weights = params.score_weights == null ? 'None' : "dict(${params.score_weights.toString().replaceAll(/\s/, '')})"
    def minimum = params.min_output_db == null ? 'None' : params.min_output_db
    return "min_output_db=${minimum}, score_weights=${weights}"
}

def scoringCliFlags() {
    def flags = params.min_output_db == null ? '' : "--min-output ${params.min_output_db}"
    if (params.score_weights != null) flags += " --score-weights ${params.score_weights.toString().replaceAll(/\s/, '')}"
    return flags
}

def userBoundsArgs() {
    return 'dict(' + ['throat_radius', 'mouth_radius', 'length', 'min_length', 'max_length', 'min_mouth_radius', 'max_mouth_radius'].collect { key ->
        "${key}=${params[key] == null ? 'None' : params[key]}"
    }.join(', ') + ')'
}

// ========================================================================
// Shared processes
// ========================================================================

process generate_geometry {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val throat_radius
    val mouth_radius
    val length
    val profile
    val num_sections

    output:
    path "horn.step"

    script:
    """
    python3 -m horn_geometry.generator \
        --throat-radius ${throat_radius} \
        --mouth-radius ${mouth_radius} \
        --length ${length} \
        --profile ${profile} \
        --num-sections ${num_sections} \
        --output-file horn.step
    """
}

process run_simulation {
    input:
    tuple path(horn_step, stageAs: "horn_input.step"), val(band_index)

    output:
    path "results_${band_index}.csv"

    script:
    def sim_length = params.length ?: 0.5
    """
    ${bandGridShell(params.min_freq, params.max_freq, band_index)}
    echo "Running band ${band_index}: \${band_low} Hz to \${band_high} Hz"
    python3 -m horn_solver.solver \
        --step-file ${horn_step} \
        --output-file results_${band_index}.csv \
        --min-freq \${band_low} \
        --max-freq \${band_high} \
        --num-intervals \${band_points} \
        --length ${sim_length} \
        --mesh-size ${params.mesh_size} --max-mesh-nodes ${params.max_mesh_nodes} \
        --radiation-model ${radiationModel()} \
        --loss-model ${params.loss_model} --flange-width ${params.flange_width} \
        --element-degree ${params.element_degree} ${params.minimum_wall_scale != null ? "--minimum-wall-scale " + params.minimum_wall_scale : ""}
    """
}

process merge_results {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(csv_files)

    output:
    path "final_results.csv"

    script:
    """
    python3 -m horn_analysis.merge --files results_*.csv \
        --num-bands ${params.num_bands} --min-freq ${params.min_freq} --max-freq ${params.max_freq} \
        ${gridFlags()} \
        --output final_results.csv
    """
}

process extract_kpis {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path final_csv

    output:
    path "kpis.json"

    script:
    """
    python3 -m horn_analysis.kpi ${final_csv} --output kpis.json
    """
}

process generate_plots {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path final_csv

    output:
    path "frequency_response.png"

    script:
    """
    python3 -m horn_analysis.plotter         ${final_csv}         frequency_response.png
    """
}

process generate_impedance_plot {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path final_csv

    output:
    path "impedance.png"

    script:
    """
    python3 -m horn_analysis.impedance_plot ${final_csv} impedance.png
    """
}

process generate_phase_plot {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path final_csv

    output:
    path "phase_response.png"

    script:
    """
    python3 -m horn_analysis.phase_plot ${final_csv} phase_response.png --group-delay
    """
}

process generate_dashboard {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path final_csv

    output:
    path "dashboard.png"

    script:
    """
    python3 -m horn_analysis.dashboard ${final_csv} dashboard.png
    """
}

process couple_with_driver {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path final_csv
    val throat_radius
    val profile
    path drivers_db, stageAs: 'driver_database'

    output:
    path "coupled_spl.csv"
    path "coupled_spl.png"
    path "driver_horn_kpis.json"

    script:
    def throat_flag = !params.step_file && throat_radius != null ? "--throat-radius ${throat_radius}" : ""
    def driver_id_arg = "'" + params.driver_id.toString().replace("'", "'\"'\"'") + "'"
    """
    python3 -m horn_analysis.couple_single ${params.step_file ? "--imported-geometry" : ""} \
        --solver-csv ${final_csv} \
        --drivers-db ${drivers_db} \
        --driver-id ${driver_id_arg} \
        --voltage ${params.voltage_rms} \
        ${throat_flag} \
        --profile ${profile} \
        --output-csv coupled_spl.csv \
        --output-png coupled_spl.png \
        --output-kpis driver_horn_kpis.json
    """
}

process render_horn_3d {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val throat_radius
    val mouth_radius
    val length
    val profile

    output:
    path "horn_3d.png"

    script:
    """
    python3 -m horn_analysis.horn_render ${params.step_file ? "--imported-geometry" : ""} \
        horn_3d.png \
        --throat-radius ${throat_radius} \
        --mouth-radius ${mouth_radius} \
        --length ${length} \
        --profile ${profile}
    """
}

process generate_single_report {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path kpis_json
    path final_csv
    path spl_png
    path impedance_png
    path phase_png
    path dashboard_png
    path horn_3d_png, stageAs: 'horn_render_input.png'

    output:
    path "single_report.html"

    script:
    def throat_r = params.throat_radius ?: 0.05
    def mouth_r = params.mouth_radius ?: 0.2
    def horn_len = params.length ?: 0.5
    """
    horn-single-report ${params.step_file ? "--imported-geometry" : ""} \
        --kpis ${kpis_json} \
        --final-csv ${final_csv} \
        --throat-radius ${throat_r} \
        --mouth-radius ${mouth_r} \
        --length ${horn_len} \
        --profile ${params.profile} \
        --spl-png ${spl_png} \
        --impedance-png ${impedance_png} \
        --phase-png ${phase_png} \
        --dashboard-png ${dashboard_png} \
        --horn-3d-png ${horn_3d_png} \
        --output single_report.html
    """
}

process generate_single_report_with_driver {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path kpis_json
    path final_csv
    path spl_png
    path impedance_png
    path phase_png
    path dashboard_png
    path horn_3d_png, stageAs: 'horn_render_input.png'
    path coupled_png
    path coupled_kpis_json

    output:
    path "single_report.html"

    script:
    def throat_r = params.throat_radius ?: 0.05
    def mouth_r = params.mouth_radius ?: 0.2
    def horn_len = params.length ?: 0.5
    """
    horn-single-report ${params.step_file ? "--imported-geometry" : ""} \
        --kpis ${kpis_json} \
        --final-csv ${final_csv} \
        --throat-radius ${throat_r} \
        --mouth-radius ${mouth_r} \
        --length ${horn_len} \
        --profile ${params.profile} \
        --spl-png ${spl_png} \
        --impedance-png ${impedance_png} \
        --phase-png ${phase_png} \
        --dashboard-png ${dashboard_png} \
        --horn-3d-png ${horn_3d_png} \
        --coupled-png ${coupled_png} \
        --coupled-kpis ${coupled_kpis_json} \
        --output single_report.html
    """
}

// ========================================================================
// Unified auto mode processes
// ========================================================================

process prescreen_drivers {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    val target_f_low
    val target_f_high
    path drivers_db, stageAs: 'driver_database'

    output:
    path "prescreen_result.json"

    script:
    def mouth_flag = params.mouth_radius != null ? "--mouth-radius ${params.mouth_radius}" : ""
    def length_flag = params.length != null ? "--length ${params.length}" : ""
    def min_dia_flag = params.min_diameter != null ? "--min-diameter ${params.min_diameter}" : ""
    def max_dia_flag = params.max_diameter != null ? "--max-diameter ${params.max_diameter}" : ""
    def throat_fractions_flag = params.throat_fractions != null ? "--throat-fractions ${params.throat_fractions}" : ""
    """
    python3 -m horn_analysis.prescreen \
        --drivers-db ${drivers_db} \
        --target-f-low ${target_f_low} \
        --target-f-high ${target_f_high} \
        ${mouth_flag} \
        ${length_flag} \
        ${min_dia_flag} \
        ${max_dia_flag} \
        ${throat_fractions_flag} \
        --output prescreen_result.json
    """
}

process report_no_drivers {
    publishDir "${params.outdir}/auto", mode: 'copy'
    input:
    tuple val(empty_reason), path(empty_metadata)
    output:
    path "report/*"
    path "ranked_results.json"
    script:
    """
    python3 -c '
import json
from pathlib import Path
from horn_analysis.auto_report import generate_auto_report
from horn_analysis.scoring import TargetSpec, describe_scoring
target = TargetSpec(${params.target_f_low}, ${params.target_f_high}, voltage_rms=${params.voltage_rms}, observation_distance_m=${params.observation_distance}, max_ripple_db=${params.max_ripple_db}, max_compression_ratio=${params.max_compression_ratio}, ${scoringTargetArgs()}, test_power_w=${params.test_power_w ?: 'None'})
reasons = {"no_drivers_passed_prescreen": "No drivers passed the initial screening for this request.", "size_constraints_exclude_search_range": "No geometry remains within the requested size limits and current search range.", "no_mouth_larger_than_throat": "The searched mouth sizes are no larger than the throat.", "profile_geometry_infeasible": "Every searched geometry is impossible for its profile, for example a tractrix shorter than its natural length.", "no_geometry_in_radiation_domain": "No searched geometry lies within the selected radiation model domain."}
metadata = json.loads(Path("${empty_metadata}").read_text())
rejected = metadata.get("rankings", [])
count = metadata.get("total_evaluated", 0)
scored = metadata.get("total_pairs", 0)
generate_auto_report([], {}, {}, None, target, "report", total_candidates=0, total_scored=0, lem_results=metadata if "rankings" in metadata else None, no_feasible_reason=reasons["${empty_reason}"])
Path("ranked_results.json").write_text(json.dumps({"status":"no_feasible_design", "reason":"${empty_reason}", "results":[], "rejected":rejected, "total_candidates":0, "total_scored":0, "analytical_candidates_evaluated":count, "analytical_pairs_evaluated":scored, "validation_status":"independent_validation_pending", "scoring":describe_scoring(target, None, ${userBoundsArgs()})}, indent=2))
'
    """
}

process derive_auto_geometry {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    path prescreen_json

    output:
    tuple path("candidates.csv"), path("design.json")

    script:
    def mouth_flag = params.mouth_radius != null ? "--mouth-radius ${params.mouth_radius}" : ""
    def length_flag = params.length != null ? "--length ${params.length}" : ""
    def throat_flag = params.throat_radius != null ? "--throat-radius ${params.throat_radius}" : ""
    def max_length_flag = params.max_length != null ? "--max-length ${params.max_length}" : ""
    def min_length_flag = params.min_length != null ? "--min-length ${params.min_length}" : ""
    def max_mr_flag = params.max_mouth_radius != null ? "--max-mouth-radius ${params.max_mouth_radius}" : ""
    def min_mr_flag = params.min_mouth_radius != null ? "--min-mouth-radius ${params.min_mouth_radius}" : ""
    def profiles_flag = params.profiles != null ? "--profiles ${params.profiles.toString().replaceAll(/\s/, '')}" : ""
    """
    python3 -m horn_core.geometry_designer \
        --target-f-low ${params.target_f_low} \
        --target-f-high ${params.target_f_high} \
        --prescreen-json ${prescreen_json} \
        --num-mouth-radii ${params.num_mouth_radii} \
        --num-lengths ${params.num_lengths} \
        ${throat_flag} \
        ${mouth_flag} \
        ${length_flag} \
        ${max_length_flag} \
        ${min_length_flag} \
        ${max_mr_flag} \
        ${min_mr_flag} \
        ${profiles_flag} \
        --output candidates.csv \
        --design-json design.json
    """
}

process lem_prescreen {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    path candidates_csv
    path prescreen_json
    path drivers_db, stageAs: 'driver_database'
    path design_json

    output:
    tuple path("lem_results.json"), path("lem_filtered_candidates.csv")

    script:
    """
    python3 -m horn_analysis.lem_prescreen \
        --candidates-csv ${candidates_csv} \
        --prescreen-json ${prescreen_json} \
        --drivers-db ${drivers_db} \
        --design-json ${design_json} \
        --target-f-low ${params.target_f_low} \
        --target-f-high ${params.target_f_high} \
        --radiation-model ${radiationModel()} --loss-model ${params.loss_model} --flange-width ${params.flange_width} \
        --voltage ${params.voltage_rms} --distance ${params.observation_distance} \
        --max-ripple ${params.max_ripple_db} --max-compression ${params.max_compression_ratio} ${scoringCliFlags()} \
        ${params.test_power_w ? "--test-power " + params.test_power_w : ""} \
        --top-n ${params.lem_top_n} \
        --output lem_results.json \
        --filtered-csv lem_filtered_candidates.csv
    """
}

process generate_candidate_geometry {
    publishDir "${params.outdir}/auto/geometry", mode: 'copy'
    input:
    tuple val(candidate_id), val(profile), val(throat_radius), val(mouth_radius), val(length)

    output:
    tuple val(candidate_id), val(profile), val(mouth_radius), val(length), path("horn_${candidate_id}.step")

    script:
    """
    python3 -m horn_geometry.generator \
        --throat-radius ${throat_radius} \
        --mouth-radius ${mouth_radius} \
        --length ${length} \
        --profile ${profile} \
        --num-sections ${params.num_sections} \
        --output-file horn_${candidate_id}.step
    """
}

process run_candidate_simulation {
    // Real solver faults stop the run. A per-candidate mesh or CAD rejection
    // (node limit, throat-area guard, non-axisymmetric CAD) exits 0 with a
    // rejection record instead, and the candidate is listed as rejected.
    errorStrategy 'terminate'

    input:
    tuple val(candidate_id), val(profile), val(mouth_radius), val(length), path(horn_step, stageAs: "horn_input.step"), val(band_index), val(sim_min_freq), val(sim_max_freq)

    output:
    tuple val(candidate_id), path("results_${candidate_id}_${band_index}.{csv,rejected.json}")

    script:
    """
    ${bandGridShell(sim_min_freq, sim_max_freq, band_index)}
    echo "Running ${candidate_id} band ${band_index}: \${band_low} Hz to \${band_high} Hz"
    python3 -m horn_solver.solver \
        --step-file ${horn_step} \
        --output-file results_${candidate_id}_${band_index}.csv \
        --rejection-file results_${candidate_id}_${band_index}.rejected.json \
        --min-freq \${band_low} \
        --max-freq \${band_high} \
        --num-intervals \${band_points} \
        --length ${length} \
        --mesh-size ${params.mesh_size} --max-mesh-nodes ${params.max_mesh_nodes} \
        --radiation-model ${radiationModel()} \
        --loss-model ${params.loss_model} --flange-width ${params.flange_width} \
        --element-degree ${params.element_degree} ${params.minimum_wall_scale != null ? "--minimum-wall-scale " + params.minimum_wall_scale : ""}
    """
}

process merge_candidate_results {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    tuple val(candidate_id), path(csv_files)

    output:
    tuple val(candidate_id), path("${candidate_id}_{results.csv,rejected.json}")

    script:
    """
    if compgen -G 'results_${candidate_id}_*.rejected.json' > /dev/null; then
        python3 -m horn_analysis.merge --rejections results_${candidate_id}_*.rejected.json \
            --candidate-id ${candidate_id} --output ${candidate_id}_rejected.json
    else
        python3 -m horn_analysis.merge --files results_${candidate_id}_*.csv \
            --num-bands ${params.num_bands} \
            --min-freq ${params.target_f_low / Math.sqrt(2)} \
            --max-freq ${params.target_f_high * Math.sqrt(2)} \
            ${gridFlags()} \
            --output ${candidate_id}_results.csv
    fi
    """
}

process score_and_rank {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    path solver_csvs
    path prescreen_json
    path drivers_db, stageAs: 'driver_database'
    path candidates_csv
    path design_json

    output:
    path "ranked_results.json"

    script:
    """
    python3 -c "
import json, glob, csv
from pathlib import Path
from horn_drivers.loader import load_drivers
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec, describe_scoring

prescreen = json.loads(Path('${prescreen_json}').read_text())
throat_radius = prescreen['throat_radius_m']
target = TargetSpec(f_low_hz=${params.target_f_low}, f_high_hz=${params.target_f_high}, voltage_rms=${params.voltage_rms}, observation_distance_m=${params.observation_distance}, max_ripple_db=${params.max_ripple_db}, max_compression_ratio=${params.max_compression_ratio}, ${scoringTargetArgs()}, test_power_w=${params.test_power_w ?: 'None'})

# Load only pre-screened drivers
all_drivers = load_drivers('${drivers_db}')
driver_ids = set(prescreen['drivers'])
drivers = [d for d in all_drivers if d.driver_id in driver_ids]

# Build candidate lookup for geometry annotation
candidates_lookup = {}
with open('${candidates_csv}') as f:
    for row in csv.DictReader(f):
        candidates_lookup[row['candidate_id']] = row

csv_files = sorted(glob.glob('*_results.csv'))
total_candidates = len(csv_files)

all_results = []
for csv_path in csv_files:
    candidate_id = Path(csv_path).stem.replace('_results', '')
    cand = candidates_lookup.get(candidate_id, {})
    cand_throat = float(cand.get('throat_radius', throat_radius))
    results = rank_horn_drivers(
        solver_csv=csv_path,
        horn_label=candidate_id,
        throat_radius=cand_throat,
        drivers=drivers,
        target=target,
        top_n=max(1, len(drivers)),
        length=float(cand['length']) if cand.get('length') else None,
    )
    # Annotate each result with geometry info
    for r in results:
        r['mouth_radius'] = float(cand.get('mouth_radius', 0))
        r['length'] = float(cand.get('length', 0))
        r['throat_radius'] = float(cand.get('throat_radius', 0))
        r['profile'] = cand.get('profile', '')
    all_results.extend(results)

# Candidates rejected by the solver before producing a response
from horn_analysis.rank_pipeline import rejection_rows
solver_rejected = []
for json_path in sorted(glob.glob('*_rejected.json')):
    record = json.loads(Path(json_path).read_text())
    solver_rejected.extend(rejection_rows(record, candidates_lookup.get(record['candidate_id'], {}), driver_ids))

# Sort all by composite score and take overall top N
total_scored = len(all_results)
rejected = [r for r in all_results if not r.get('model_feasible', False)] + solver_rejected
all_results = [r for r in all_results if r.get('model_feasible', False)]
all_results.sort(key=lambda r: r['composite_score'], reverse=True)
from horn_analysis.search import annotate_comparable_candidates
all_results = annotate_comparable_candidates(all_results[:${params.top_n}])

output = {
    'total_candidates': total_candidates,
    'solver_rejected_candidates': sorted({r['horn_label'] for r in solver_rejected}),
    'total_scored': total_scored,
    'results': all_results,
    'rejected': rejected,
    'status': 'experimental_candidates' if all_results else 'no_feasible_design',
    'validation_status': 'independent_validation_pending',
    'scoring': describe_scoring(target, json.loads(Path('${design_json}').read_text()), ${userBoundsArgs()}),
}
Path('ranked_results.json').write_text(json.dumps(output, indent=2))
print(f'Ranked {total_scored} driver-horn combinations ({total_candidates} geometries)')
"
    """
}

process refine_ranked_candidates {
    publishDir "${params.outdir}/auto", mode: 'copy'
    input:
    path ranked_json
    path solver_csvs
    path drivers_db, stageAs: 'driver_database'
    path prescreen_json
    path design_json
    output:
    path "refinement/ranked_results.json"
    path "refinement/*_results.csv"
    path "refinement/search_audit.json"
    path "refinement/*.step", optional: true
    script:
    """
    python3 -m horn_analysis.refine \
        --ranked-json ${ranked_json} --solver-csvs ${solver_csvs instanceof List ? solver_csvs.join(' ') : solver_csvs} \
        --drivers-db ${drivers_db} --prescreen-json ${prescreen_json} --design-json ${design_json} \
        --budget ${params.refinement_budget} --top-n ${params.top_n} --num-frequencies ${params.num_intervals} \
        --num-bands ${params.num_bands} --min-points-per-octave ${params.min_points_per_octave} --num-sections ${params.num_sections} \
        --mesh-size ${params.mesh_size} --max-mesh-nodes ${params.max_mesh_nodes} --radiation-model ${radiationModel()} \
        --loss-model ${params.loss_model} --flange-width ${params.flange_width} \
        --element-degree ${params.element_degree} ${params.minimum_wall_scale != null ? "--minimum-wall-scale " + params.minimum_wall_scale : ""} \
        --voltage ${params.voltage_rms} --distance ${params.observation_distance} \
        --max-ripple ${params.max_ripple_db} --max-compression ${params.max_compression_ratio} ${scoringCliFlags()} \
        ${params.test_power_w ? "--test-power " + params.test_power_w : ""} \
        --output-dir refinement
    """
}

process generate_auto_report {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    path ranked_json
    path solver_csvs
    path drivers_db, stageAs: 'driver_database'
    path prescreen_json
    path design_json
    path lem_results_json

    output:
    path "report/auto_ranking.json"
    path "report/auto_comparison.png"
    path "report/auto_summary.txt"
    path "report/auto_report.html"
    path "report/coupled_*.csv", optional: true
    path "report/auto_electrical.png", optional: true
    path "report/directivity_*.csv", optional: true
    path "report/auto_beamwidth.png", optional: true

    script:
    """
    python3 -c "
import json, glob
from pathlib import Path
from horn_drivers.loader import load_drivers
from horn_analysis.scoring import TargetSpec
from horn_analysis.auto_report import generate_auto_report

prescreen = json.loads(Path('${prescreen_json}').read_text())
throat_radius = prescreen['throat_radius_m']
design = json.loads(Path('${design_json}').read_text())

ranked_data = json.loads(Path('${ranked_json}').read_text())
all_ranked = ranked_data['results']
total_candidates = ranked_data.get('total_candidates', len(all_ranked))
total_scored = ranked_data.get('total_scored', len(all_ranked))

# Load LEM prescreen results
lem_results = json.loads(Path('${lem_results_json}').read_text())

solver_csvs = {}
for csv_path in sorted(glob.glob('*_results.csv')):
    candidate_id = Path(csv_path).stem.replace('_results', '')
    solver_csvs[candidate_id] = csv_path

driver_list = load_drivers('${drivers_db}')
driver_ids = set(prescreen['drivers'])
drivers = {d.driver_id: d for d in driver_list if d.driver_id in driver_ids}

target = TargetSpec(f_low_hz=${params.target_f_low}, f_high_hz=${params.target_f_high}, voltage_rms=${params.voltage_rms}, observation_distance_m=${params.observation_distance}, max_ripple_db=${params.max_ripple_db}, max_compression_ratio=${params.max_compression_ratio}, ${scoringTargetArgs()}, test_power_w=${params.test_power_w ?: 'None'})

generate_auto_report(
    all_ranked=all_ranked,
    solver_csvs=solver_csvs,
    drivers=drivers,
    throat_radius=throat_radius,
    target=target,
    output_dir='report',
    top_n=5,
    derived_geometry=design,
    total_candidates=total_candidates,
    total_scored=total_scored,
    lem_results=lem_results,
    scoring=ranked_data.get('scoring'),
)
"
    """
}

// ========================================================================
// Workflow definitions
// ========================================================================

workflow single {
    if (params.step_file && (params.length == null || !Double.isFinite(params.length as double) || (params.length as double) <= 0)) {
        error "Imported STEP geometry requires an explicit positive --length matching its outlet coordinate"
    }
    // Apply defaults for single mode when params are null
    def throat_r = params.throat_radius ?: 0.05
    def mouth_r = params.mouth_radius ?: 0.2
    def horn_len = params.length ?: 0.5

    // 1. Generate geometry once — or use a user-supplied STEP file if provided
    if (params.step_file) {
        ch_step_file = Channel.value(file(params.step_file))
    } else {
        ch_step_file = generate_geometry(
            throat_r,
            mouth_r,
            horn_len,
            params.profile,
            params.num_sections
        )
    }

    // 2. Create a channel of band indices
    ch_band_indices = Channel.from(0..<params.num_bands)

    // 3. Combine the geometry file with the band indices
    ch_sim_inputs = ch_step_file.combine(ch_band_indices)

    // 4. Run simulations in parallel
    ch_band_results = run_simulation(ch_sim_inputs)

    // 5. Collect and merge results
    ch_merged_results = merge_results(ch_band_results.collect())

    // 6. Extract KPIs
    extract_kpis(ch_merged_results)

    // 7. Plot final results
    generate_plots(ch_merged_results)

    // 8. Impedance and phase plots
    generate_impedance_plot(ch_merged_results)
    generate_phase_plot(ch_merged_results)

    // 9. Combined dashboard
    generate_dashboard(ch_merged_results)

    // 9b. Optional driver coupling — produces coupled_spl.csv/.png + KPI JSON
    ch_coupled_png = null
    ch_coupled_kpis = null
    if (params.driver_id) {
        couple_outputs = couple_with_driver(
            ch_merged_results,
            throat_r,
            params.step_file ? "imported_STEP" : params.profile,
            file(params.drivers_db),
        )
        ch_coupled_png  = couple_outputs[1]
        ch_coupled_kpis = couple_outputs[2]
    }

    // 10. 3D horn geometry render — use a pre-rendered PNG if supplied,
    //     otherwise generate one parametrically from the profile.
    if (params.horn_3d_png) {
        ch_horn_3d_png = Channel.value(file(params.horn_3d_png))
    } else {
        ch_horn_3d_png = render_horn_3d(
            throat_r,
            mouth_r,
            horn_len,
            params.profile
        )
    }

    if (params.driver_id) {
        // 12. HTML report with driver-coupled section
        generate_single_report_with_driver(
            extract_kpis.out, ch_merged_results,
            generate_plots.out, generate_impedance_plot.out,
            generate_phase_plot.out, generate_dashboard.out, ch_horn_3d_png,
            ch_coupled_png, ch_coupled_kpis,
        )
    } else {
        // 12. HTML report without driver coupling
        generate_single_report(
            extract_kpis.out, ch_merged_results,
            generate_plots.out, generate_impedance_plot.out,
            generate_phase_plot.out, generate_dashboard.out, ch_horn_3d_png,
        )
    }
}

workflow auto {
    // 1. Pre-screen drivers
    ch_drivers_db = Channel.value(file(params.drivers_db))
    ch_prescreen = prescreen_drivers(
        params.target_f_low,
        params.target_f_high,
        ch_drivers_db,
    )

    ch_driver_branches = ch_prescreen.branch { path ->
        eligible: new groovy.json.JsonSlurper().parse(path).count > 0
        empty: true
    }
    ch_prescreen = ch_driver_branches.eligible

    // 2. Derive geometry grid from frequency band + prescreen throat radii
    //    Fixed params (mouth_radius, length) are passed via CLI flags in the process
    ch_geom_derived = derive_auto_geometry(ch_prescreen)
    ch_geometry_branches = ch_geom_derived.branch { csv, json ->
        eligible: new groovy.json.JsonSlurper().parse(json).candidate_count > 0
        empty: true
    }
    ch_empty_drivers = ch_driver_branches.empty.map { path -> tuple('no_drivers_passed_prescreen', path) }
    ch_empty_geometry = ch_geometry_branches.empty.map { csv, json ->
        def data = new groovy.json.JsonSlurper().parse(json)
        tuple(data.reason, json)
    }
    ch_candidates_csv = ch_geometry_branches.eligible.map { csv, json -> csv }
    ch_design_json = ch_geometry_branches.eligible.map { csv, json -> json }

    // 3. LEM/Webster prescreening — score all candidates analytically,
    //    pass only the top N to expensive STEP + FEM stages
    ch_lem = lem_prescreen(
        ch_candidates_csv,
        ch_prescreen,
        ch_drivers_db,
        ch_design_json,
    )
    ch_lem_branches = ch_lem.branch { json, csv ->
        eligible: new groovy.json.JsonSlurper().parse(json).filtered_candidate_ids.size() > 0
        empty: true
    }
    ch_empty_lem = ch_lem_branches.empty.map { json, csv -> tuple('no_geometry_in_radiation_domain', json) }
    report_no_drivers(ch_empty_drivers.mix(ch_empty_geometry).mix(ch_empty_lem))
    ch_lem_results = ch_lem_branches.eligible.map { json, csv -> json }
    ch_filtered_csv = ch_lem_branches.eligible.map { json, csv -> csv }

    // 4. Parse filtered candidates CSV into channel of tuples
    ch_candidates = ch_filtered_csv
        .splitCsv(header: true)
        .map { row ->
            tuple(row.candidate_id, row.profile, row.throat_radius as double,
                  row.mouth_radius as double, row.length as double)
        }

    // 5. Generate STEP geometry for filtered candidates only
    ch_geometries = generate_candidate_geometry(ch_candidates)

    // 6. Read sim freq range from design.json and combine with band indices
    ch_sim_range = ch_design_json.map { json_file ->
        def data = new groovy.json.JsonSlurper().parse(json_file)
        tuple(data.sim_freq_range[0] as double, data.sim_freq_range[1] as double)
    }

    ch_band_indices = Channel.from(0..<params.num_bands)
    ch_sim_inputs = ch_geometries
        .combine(ch_band_indices)
        .combine(ch_sim_range)

    // 7. Run FEM simulations (filtered candidates x bands)
    ch_band_results = run_candidate_simulation(ch_sim_inputs)

    // 8. Group by candidate_id and merge bands
    ch_grouped = ch_band_results.groupTuple()
    ch_merged = merge_candidate_results(ch_grouped)

    // 9. Score and rank all driver-horn combinations
    ch_all_csvs = ch_merged.map { candidate_id, csv -> csv }.collect()
    ch_ranked = score_and_rank(
        ch_all_csvs,
        ch_prescreen,
        ch_drivers_db,
        ch_filtered_csv,
        ch_design_json,
    )

    // 10. Refine promising geometries inside the solver container.
    ch_final_ranked = ch_ranked
    ch_final_csvs = ch_all_csvs
    if (params.refinement_budget > 0) {
        ch_refined = refine_ranked_candidates(ch_ranked, ch_all_csvs, ch_drivers_db, ch_prescreen, ch_design_json)
        ch_final_ranked = ch_refined[0]
        ch_final_csvs = ch_refined[1]
    }
    // 11. Generate report with design summary + LEM stats
    generate_auto_report(
        ch_final_ranked,
        ch_final_csvs,
        ch_drivers_db,
        ch_prescreen,
        ch_design_json,
        ch_lem_results,
    )
}

workflow {
    if (params.outdir == null || params.outdir instanceof Boolean || !(params.outdir as String).trim()) error "Choose --outdir for direct Nextflow runs, or use scripts/run_pipeline.py for isolated outputs and checked resume"
    if (!(params.loss_model in ['lossless', 'boundary_layer'])) error "Unknown loss model"
    if (!(params.element_degree in [1, 2])) error "Element degree must be 1 or 2"
    if (!Double.isFinite(params.flange_width as double) || params.flange_width < 0) error "Invalid flange width"
    if (params.loss_model == 'boundary_layer' && (params.minimum_wall_scale == null || !Double.isFinite(params.minimum_wall_scale as double) || params.minimum_wall_scale <= 0)) error "Boundary-layer walls require minimum_wall_scale in metres"
    if (params.mode != 'single' && params.radiation_model == 'closed') error "Closed termination is a single-case impedance validation boundary, not a radiating horn"
    if (params.radiation_model != null && !(params.radiation_model in ['plane_wave', 'flanged_piston', 'modal_baffled', 'unflanged_piston', 'finite_flange', 'closed'])) error "Unknown radiation model: choose plane_wave, flanged_piston, modal_baffled, unflanged_piston, finite_flange or closed"
    if (params.profiles != null && !(params.profiles.toString().replaceAll(/\s/, '') ==~ /(conical|exponential|hyperbolic|tractrix|os|scaled_tractrix|cd)(,(conical|exponential|hyperbolic|tractrix|os|scaled_tractrix|cd))*/)) error "profiles must be a comma list of conical, exponential, hyperbolic, tractrix, os, scaled_tractrix, cd"
    if (params.radiation_model == 'modal_baffled' && (params.loss_model != 'lossless' || params.element_degree != 1)) error "Experimental modal aperture requires lossless P1"
    if (!(params.mode in ['single', 'auto', 'fullauto'])) error "Unknown mode: ${params.mode}"
    def low = params.mode == 'single' ? params.min_freq : params.target_f_low
    def high = params.mode == 'single' ? params.max_freq : params.target_f_high
    if (!(low > 0 && high > low && Double.isFinite(high as double))) error "Require 0 < low frequency < high frequency"
    for (key in ['voltage_rms', 'observation_distance', 'max_ripple_db', 'max_compression_ratio', 'mesh_size', 'max_mesh_nodes', 'num_bands', 'num_intervals', 'num_sections', 'num_mouth_radii', 'num_lengths', 'lem_top_n', 'top_n']) {
        if (!(params[key] > 0 && Double.isFinite(params[key] as double))) error "${key} must be finite and positive"
    }
    for (key in ['max_mesh_nodes', 'num_bands', 'num_intervals', 'num_sections', 'num_mouth_radii', 'num_lengths', 'lem_top_n', 'top_n', 'refinement_budget']) {
        if (!Double.isFinite(params[key] as double) || (params[key] as double) != Math.floor(params[key] as double)) error "${key} must be an integer"
    }
    if (params.num_intervals < 2 || params.num_sections < 2) error "At least two frequency points and geometry sections are required"
    if (params.refinement_budget < 0) error "Refinement budget must be nonnegative"
    if (!(params.min_points_per_octave >= 0 && Double.isFinite(params.min_points_per_octave as double))) error "min_points_per_octave must be finite and nonnegative"
    if (params.min_output_db != null && !Double.isFinite(params.min_output_db as double)) error "min_output_db must be finite"
    if (params.score_weights != null && !(params.score_weights.toString().replaceAll(/\s/, '') ==~ /(sensitivity|ripple|extension|size)=[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?(,(sensitivity|ripple|extension|size)=[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)*/)) error "score_weights must look like sensitivity=0.625,ripple=0.25,extension=0.125,size=0"
    for (key in ['throat_radius', 'mouth_radius', 'length', 'min_length', 'max_length', 'min_mouth_radius', 'max_mouth_radius', 'test_power_w']) {
        if (params[key] != null && !(params[key] > 0 && Double.isFinite(params[key] as double))) error "${key} must be finite and positive"
    }
    def outputRoot = new File(params.outdir as String)
    outputRoot.mkdirs()
    def resolvedParameters = [:]
    params.each { key, value -> resolvedParameters[key] = value }
    resolvedParameters.radiation_model = radiationModel()
    if (params.mode == 'single' && !params.step_file) {
        [throat_radius: 0.05, mouth_radius: 0.2, length: 0.5].each { key, value ->
            if (resolvedParameters[key] == null) resolvedParameters[key] = value
        }
    }
    new File(outputRoot, 'resolved_specification.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([specification_version: 2, parameters: resolvedParameters]))
    if (params.mode == "fullauto" || params.mode == "auto") {
        auto()
    } else {
        single()
    }
}
