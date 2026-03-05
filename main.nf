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
params.profile = "conical" // Horn flare profile: conical, exponential, hyperbolic, tractrix, os, lecleach, cd
params.num_sections = 20   // Number of cross-sections for lofting

// Simulation Settings
params.min_freq = 500      // Minimum frequency for the sweep in Hz
params.max_freq = 8000     // Maximum frequency for the sweep in Hz
params.num_intervals = 100 // Number of frequency steps in the sweep
params.mesh_size = 0.01    // Target mesh element size in meters

// Radiation impedance model at the horn mouth
params.radiation_model = "plane_wave"  // plane_wave, flanged_piston, unflanged_piston, bem

// Execution Settings
params.num_bands = 8       // Number of parallel jobs for the solver
params.outdir = "./results"

// Directivity (opt-in, single mode only, requires BEM)
params.directivity = false

// Auto mode settings
params.target_f_low = 500
params.target_f_high = 4000
params.drivers_db = "data/drivers"
params.top_n = 10
params.num_mouth_radii = 3 // Mouth radius grid points (when mouth_radius is null)
params.num_lengths = 3     // Length grid points (when length is null)
params.lem_top_n = 3       // Number of top candidates to pass from LEM prescreen to FEM
params.min_diameter = null  // Optional: minimum driver nominal diameter (inches)
params.max_diameter = null  // Optional: maximum driver nominal diameter (inches)

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
    tuple path(horn_step), val(band_index)

    output:
    path "results_${band_index}.csv"

    script:
    def sim_length = params.length ?: 0.5
    def band_width = (params.max_freq - params.min_freq) / (params.num_bands as double)
    def min_f = params.min_freq + band_width * band_index
    def max_f = params.min_freq + band_width * (band_index + 1)
    def num_intervals_per_band = Math.ceil(params.num_intervals / (params.num_bands as double)) as int
    """
    echo "Running band ${band_index}: ${min_f} Hz to ${max_f} Hz"
    python3 -m horn_solver.solver \
        --step-file ${horn_step} \
        --output-file results_${band_index}.csv \
        --min-freq ${min_f} \
        --max-freq ${max_f} \
        --num-intervals ${num_intervals_per_band} \
        --length ${sim_length} \
        --mesh-size ${params.mesh_size} \
        --radiation-model ${params.radiation_model}
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
    python3 -c "import pandas as pd; import glob; all_files = glob.glob('results_*.csv'); df = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True); df.sort_values(by='frequency').to_csv('final_results.csv', index=False)"
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
    python3 -m horn_analysis.horn_render \
        horn_3d.png \
        --throat-radius ${throat_radius} \
        --mouth-radius ${mouth_radius} \
        --length ${length} \
        --profile ${profile}
    """
}

process run_simulation_directivity {
    input:
    tuple path(horn_step), val(band_index)

    output:
    path "directivity_${band_index}.csv"

    script:
    def sim_length = params.length ?: 0.5
    def band_width = (params.max_freq - params.min_freq) / (params.num_bands as double)
    def min_f = params.min_freq + band_width * band_index
    def max_f = params.min_freq + band_width * (band_index + 1)
    def num_intervals_per_band = Math.ceil(params.num_intervals / (params.num_bands as double)) as int
    """
    echo "Running directivity band ${band_index}: ${min_f} Hz to ${max_f} Hz"
    python3 -m horn_solver.solver \
        --step-file ${horn_step} \
        --output-file solver_directivity_${band_index}.csv \
        --min-freq ${min_f} \
        --max-freq ${max_f} \
        --num-intervals ${num_intervals_per_band} \
        --length ${sim_length} \
        --mesh-size ${params.mesh_size} \
        --radiation-model bem \
        --compute-directivity \
        --directivity-file directivity_${band_index}.csv
    """
}

process merge_directivity_results {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(csv_files)

    output:
    path "directivity.csv"

    script:
    """
    python3 -c "import pandas as pd; import glob; all_files = glob.glob('directivity_*.csv'); df = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True); df = df.drop_duplicates(subset=['frequency', 'theta_deg'], keep='first'); df.sort_values(by=['frequency', 'theta_deg']).to_csv('directivity.csv', index=False)"
    """
}

process generate_directivity_plots {
    publishDir "${params.outdir}/directivity", mode: 'copy'

    input:
    path directivity_csv

    output:
    path "polar_directivity.png"
    path "directivity_contour.png"
    path "beamwidth.png"
    path "directivity_index.png"

    script:
    """
    python3 -m horn_analysis.directivity_plot ${directivity_csv} --output-dir .
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
    path horn_3d_png

    output:
    path "single_report.html"

    script:
    def throat_r = params.throat_radius ?: 0.05
    def mouth_r = params.mouth_radius ?: 0.2
    def horn_len = params.length ?: 0.5
    """
    horn-single-report \
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

process generate_single_report_with_directivity {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path kpis_json
    path final_csv
    path spl_png
    path impedance_png
    path phase_png
    path dashboard_png
    path horn_3d_png
    path polar_png
    path contour_png
    path beamwidth_png
    path di_png

    output:
    path "single_report.html"

    script:
    def throat_r = params.throat_radius ?: 0.05
    def mouth_r = params.mouth_radius ?: 0.2
    def horn_len = params.length ?: 0.5
    """
    horn-single-report \
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
        --polar-png ${polar_png} \
        --contour-png ${contour_png} \
        --beamwidth-png ${beamwidth_png} \
        --di-png ${di_png} \
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
    path drivers_db

    output:
    path "prescreen_result.json"

    script:
    def mouth_flag = params.mouth_radius != null ? "--mouth-radius ${params.mouth_radius}" : ""
    def length_flag = params.length != null ? "--length ${params.length}" : ""
    def min_dia_flag = params.min_diameter != null ? "--min-diameter ${params.min_diameter}" : ""
    def max_dia_flag = params.max_diameter != null ? "--max-diameter ${params.max_diameter}" : ""
    """
    python3 -m horn_analysis.prescreen \
        --drivers-db ${drivers_db} \
        --target-f-low ${target_f_low} \
        --target-f-high ${target_f_high} \
        ${mouth_flag} \
        ${length_flag} \
        ${min_dia_flag} \
        ${max_dia_flag} \
        --output prescreen_result.json
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
    """
    python3 -m horn_core.geometry_designer \
        --target-f-low ${params.target_f_low} \
        --target-f-high ${params.target_f_high} \
        --prescreen-json ${prescreen_json} \
        --num-mouth-radii ${params.num_mouth_radii} \
        --num-lengths ${params.num_lengths} \
        ${mouth_flag} \
        ${length_flag} \
        --output candidates.csv \
        --design-json design.json
    """
}

process lem_prescreen {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    path candidates_csv
    path prescreen_json
    path drivers_db
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
        --top-n ${params.lem_top_n} \
        --output lem_results.json \
        --filtered-csv lem_filtered_candidates.csv
    """
}

process generate_candidate_geometry {
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
    errorStrategy 'ignore'

    input:
    tuple val(candidate_id), val(profile), val(mouth_radius), val(length), path(horn_step), val(band_index), val(sim_min_freq), val(sim_max_freq)

    output:
    tuple val(candidate_id), path("results_${candidate_id}_${band_index}.csv")

    script:
    def band_width = (sim_max_freq - sim_min_freq) / (params.num_bands as double)
    def min_f = sim_min_freq + band_width * band_index
    def max_f = sim_min_freq + band_width * (band_index + 1)
    def num_intervals_per_band = Math.ceil(params.num_intervals / (params.num_bands as double)) as int
    """
    echo "Running ${candidate_id} band ${band_index}: ${min_f} Hz to ${max_f} Hz"
    python3 -m horn_solver.solver \
        --step-file ${horn_step} \
        --output-file results_${candidate_id}_${band_index}.csv \
        --min-freq ${min_f} \
        --max-freq ${max_f} \
        --num-intervals ${num_intervals_per_band} \
        --length ${length} \
        --mesh-size ${params.mesh_size}
    """
}

process merge_candidate_results {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    tuple val(candidate_id), path(csv_files)

    output:
    tuple val(candidate_id), path("${candidate_id}_results.csv")

    script:
    """
    python3 -c "
import pandas as pd, glob
files = glob.glob('results_${candidate_id}_*.csv')
df = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
df.sort_values(by='frequency').to_csv('${candidate_id}_results.csv', index=False)
"
    """
}

process score_and_rank {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    path solver_csvs
    path prescreen_json
    path drivers_db
    path candidates_csv

    output:
    path "ranked_results.json"

    script:
    """
    python3 -c "
import json, glob, csv
from pathlib import Path
from horn_drivers.loader import load_drivers
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec

prescreen = json.loads(Path('${prescreen_json}').read_text())
throat_radius = prescreen['throat_radius_m']
target = TargetSpec(f_low_hz=${params.target_f_low}, f_high_hz=${params.target_f_high})

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
        top_n=${params.top_n},
    )
    # Annotate each result with geometry info
    for r in results:
        r['mouth_radius'] = float(cand.get('mouth_radius', 0))
        r['length'] = float(cand.get('length', 0))
        r['throat_radius'] = float(cand.get('throat_radius', 0))
        r['profile'] = cand.get('profile', '')
    all_results.extend(results)

# Sort all by composite score and take overall top N
total_scored = len(all_results)
all_results.sort(key=lambda r: r['composite_score'], reverse=True)
all_results = all_results[:${params.top_n}]

output = {
    'total_candidates': total_candidates,
    'total_scored': total_scored,
    'results': all_results,
}
Path('ranked_results.json').write_text(json.dumps(output, indent=2))
print(f'Ranked {total_scored} driver-horn combinations ({total_candidates} geometries)')
"
    """
}

process generate_auto_report {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    path ranked_json
    path solver_csvs
    path drivers_db
    path prescreen_json
    path design_json
    path lem_results_json

    output:
    path "report/auto_ranking.json"
    path "report/auto_comparison.png"
    path "report/auto_summary.txt"
    path "report/auto_report.html"

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
drivers = {d.driver_id: d for d in driver_list}

target = TargetSpec(f_low_hz=${params.target_f_low}, f_high_hz=${params.target_f_high})

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
)
"
    """
}

// ========================================================================
// Workflow definitions
// ========================================================================

workflow single {
    // Apply defaults for single mode when params are null
    def throat_r = params.throat_radius ?: 0.05
    def mouth_r = params.mouth_radius ?: 0.2
    def horn_len = params.length ?: 0.5

    // 1. Generate geometry once
    ch_step_file = generate_geometry(
        throat_r,
        mouth_r,
        horn_len,
        params.profile,
        params.num_sections
    )

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

    // 10. 3D horn geometry render (runs in parallel with simulation)
    render_horn_3d(
        throat_r,
        mouth_r,
        horn_len,
        params.profile
    )

    // 11. Directivity (opt-in, requires BEM) -- parallelized across frequency bands
    if (params.directivity) {
        ch_dir_band_indices = Channel.from(0..<params.num_bands)
        ch_dir_inputs = ch_step_file.combine(ch_dir_band_indices)
        ch_dir_band_results = run_simulation_directivity(ch_dir_inputs)
        ch_directivity_csv = merge_directivity_results(ch_dir_band_results.collect())
        generate_directivity_plots(ch_directivity_csv)

        // 12. HTML report with directivity
        generate_single_report_with_directivity(
            extract_kpis.out, ch_merged_results,
            generate_plots.out, generate_impedance_plot.out,
            generate_phase_plot.out, generate_dashboard.out, render_horn_3d.out,
            generate_directivity_plots.out[0],
            generate_directivity_plots.out[1],
            generate_directivity_plots.out[2],
            generate_directivity_plots.out[3],
        )
    } else {
        // 12. HTML report without directivity
        generate_single_report(
            extract_kpis.out, ch_merged_results,
            generate_plots.out, generate_impedance_plot.out,
            generate_phase_plot.out, generate_dashboard.out, render_horn_3d.out,
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

    // 2. Derive geometry grid from frequency band + prescreen throat radii
    //    Fixed params (mouth_radius, length) are passed via CLI flags in the process
    ch_geom_derived = derive_auto_geometry(ch_prescreen)
    ch_candidates_csv = ch_geom_derived.map { csv, json -> csv }
    ch_design_json = ch_geom_derived.map { csv, json -> json }

    // 3. LEM/Webster prescreening — score all candidates analytically,
    //    pass only the top N to expensive STEP + FEM stages
    ch_lem = lem_prescreen(
        ch_candidates_csv,
        ch_prescreen,
        ch_drivers_db,
        ch_design_json,
    )
    ch_lem_results = ch_lem.map { json, csv -> json }
    ch_filtered_csv = ch_lem.map { json, csv -> csv }

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
    )

    // 10. Generate report with design summary + LEM stats
    generate_auto_report(
        ch_ranked,
        ch_all_csvs,
        ch_drivers_db,
        ch_prescreen,
        ch_design_json,
        ch_lem_results,
    )
}

workflow {
    if (params.mode == "fullauto" || params.mode == "auto") {
        auto()
    } else {
        single()
    }
}
