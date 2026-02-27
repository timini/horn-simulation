#!/usr/bin/env nextflow

// ========================================================================
// Mode: "single" (default) runs one horn profile; "auto" runs all 3
// profiles and ranks driver-horn combinations; "fullauto" derives
// geometry from a target frequency band and explores a grid.
// ========================================================================
params.mode = "single"

// Horn Geometry
params.throat_radius = 0.05 // Radius of the horn's throat in meters
params.mouth_radius = 0.2  // Radius of the horn's mouth in meters
params.length = 0.5        // Length of the horn in meters
params.profile = "conical" // Horn flare profile: conical, exponential, hyperbolic
params.num_sections = 20   // Number of cross-sections for lofting

// Simulation Settings
params.min_freq = 500      // Minimum frequency for the sweep in Hz
params.max_freq = 8000     // Maximum frequency for the sweep in Hz
params.num_intervals = 100 // Number of frequency steps in the sweep
params.mesh_size = 0.01    // Target mesh element size in meters

// Execution Settings
params.num_bands = 8       // Number of parallel jobs for the solver
params.outdir = "./results"

// Auto mode settings
params.target_f_low = 500
params.target_f_high = 4000
params.drivers_db = "data/drivers.json"
params.top_n = 10

// Fullauto mode settings
params.num_mouth_radii = 3 // Mouth radius grid points for fullauto
params.num_lengths = 3     // Length grid points for fullauto

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
        --length ${params.length} \
        --mesh-size ${params.mesh_size}
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

// ========================================================================
// Auto mode processes
// ========================================================================

process prescreen_drivers {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    val target_f_low
    val target_f_high
    val mouth_radius
    val length
    path drivers_db

    output:
    path "prescreen_result.json"

    script:
    """
    python3 -m horn_analysis.prescreen \
        --drivers-db ${drivers_db} \
        --target-f-low ${target_f_low} \
        --target-f-high ${target_f_high} \
        --mouth-radius ${mouth_radius} \
        --length ${length} \
        --output prescreen_result.json
    """
}

process generate_auto_geometry {
    input:
    tuple val(profile), val(throat_radius)

    output:
    tuple val(profile), path("horn_${profile}.step")

    script:
    """
    python3 -m horn_geometry.generator \
        --throat-radius ${throat_radius} \
        --mouth-radius ${params.mouth_radius} \
        --length ${params.length} \
        --profile ${profile} \
        --num-sections ${params.num_sections} \
        --output-file horn_${profile}.step
    """
}

process run_auto_simulation {
    input:
    tuple val(profile), path(horn_step), val(band_index)

    output:
    tuple val(profile), path("results_${profile}_${band_index}.csv")

    script:
    def band_width = (params.max_freq - params.min_freq) / (params.num_bands as double)
    def min_f = params.min_freq + band_width * band_index
    def max_f = params.min_freq + band_width * (band_index + 1)
    def num_intervals_per_band = Math.ceil(params.num_intervals / (params.num_bands as double)) as int
    """
    echo "Running ${profile} band ${band_index}: ${min_f} Hz to ${max_f} Hz"
    python3 -m horn_solver.solver \
        --step-file ${horn_step} \
        --output-file results_${profile}_${band_index}.csv \
        --min-freq ${min_f} \
        --max-freq ${max_f} \
        --num-intervals ${num_intervals_per_band} \
        --length ${params.length} \
        --mesh-size ${params.mesh_size}
    """
}

process merge_auto_results {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    tuple val(profile), path(csv_files)

    output:
    tuple val(profile), path("${profile}_results.csv")

    script:
    """
    python3 -c "
import pandas as pd, glob
files = glob.glob('results_${profile}_*.csv')
df = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
df.sort_values(by='frequency').to_csv('${profile}_results.csv', index=False)
"
    """
}

process score_and_rank {
    publishDir "${params.outdir}/auto", mode: 'copy'

    input:
    path solver_csvs
    path prescreen_json
    path drivers_db

    output:
    path "ranked_results.json"

    script:
    """
    python3 -c "
import json, glob
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

all_results = []
for csv_path in sorted(glob.glob('*_results.csv')):
    profile = Path(csv_path).stem.replace('_results', '')
    results = rank_horn_drivers(
        solver_csv=csv_path,
        horn_label=profile,
        throat_radius=throat_radius,
        drivers=drivers,
        target=target,
        top_n=${params.top_n},
    )
    all_results.extend(results)

# Sort all by composite score and take overall top N
all_results.sort(key=lambda r: r['composite_score'], reverse=True)
all_results = all_results[:${params.top_n}]

Path('ranked_results.json').write_text(json.dumps(all_results, indent=2))
print(f'Ranked {len(all_results)} driver-horn combinations')
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

    output:
    path "report/auto_ranking.json"
    path "report/auto_comparison.png"
    path "report/auto_summary.txt"
    path "report/auto_report.html"

    script:
    """
    python3 -c "
import json
from pathlib import Path
from horn_drivers.loader import load_drivers
from horn_analysis.scoring import TargetSpec
from horn_analysis.auto_report import generate_auto_report

prescreen = json.loads(Path('${prescreen_json}').read_text())
throat_radius = prescreen['throat_radius_m']

all_ranked = json.loads(Path('${ranked_json}').read_text())

import glob
solver_csvs = {}
for csv_path in sorted(glob.glob('*_results.csv')):
    profile = Path(csv_path).stem.replace('_results', '')
    solver_csvs[profile] = csv_path

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
    mouth_radius=${params.mouth_radius},
    horn_length=${params.length},
)
"
    """
}

// ========================================================================
// Fullauto mode processes
// ========================================================================

process derive_fullauto_geometry {
    publishDir "${params.outdir}/fullauto", mode: 'copy'

    input:
    path prescreen_json

    output:
    path "candidates.csv"
    path "design.json"

    script:
    """
    python3 -m horn_core.geometry_designer \
        --target-f-low ${params.target_f_low} \
        --target-f-high ${params.target_f_high} \
        --prescreen-json ${prescreen_json} \
        --num-mouth-radii ${params.num_mouth_radii} \
        --num-lengths ${params.num_lengths} \
        --output candidates.csv \
        --design-json design.json
    """
}

process generate_fullauto_geometry {
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

process run_fullauto_simulation {
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

process merge_fullauto_results {
    publishDir "${params.outdir}/fullauto", mode: 'copy'

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

process score_and_rank_fullauto {
    publishDir "${params.outdir}/fullauto", mode: 'copy'

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

all_results = []
for csv_path in sorted(glob.glob('*_results.csv')):
    candidate_id = Path(csv_path).stem.replace('_results', '')
    cand = candidates_lookup.get(candidate_id, {})
    results = rank_horn_drivers(
        solver_csv=csv_path,
        horn_label=candidate_id,
        throat_radius=throat_radius,
        drivers=drivers,
        target=target,
        top_n=${params.top_n},
    )
    # Annotate each result with geometry info
    for r in results:
        r['mouth_radius'] = float(cand.get('mouth_radius', 0))
        r['length'] = float(cand.get('length', 0))
        r['profile'] = cand.get('profile', '')
    all_results.extend(results)

# Sort all by composite score and take overall top N
all_results.sort(key=lambda r: r['composite_score'], reverse=True)
all_results = all_results[:${params.top_n}]

Path('ranked_results.json').write_text(json.dumps(all_results, indent=2))
print(f'Ranked {len(all_results)} driver-horn combinations')
"
    """
}

process generate_fullauto_report {
    publishDir "${params.outdir}/fullauto", mode: 'copy'

    input:
    path ranked_json
    path solver_csvs
    path drivers_db
    path prescreen_json
    path design_json

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

all_ranked = json.loads(Path('${ranked_json}').read_text())

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
)
"
    """
}

// ========================================================================
// Workflow definitions
// ========================================================================

workflow single {
    // 1. Generate geometry once
    ch_step_file = generate_geometry(
        params.throat_radius,
        params.mouth_radius,
        params.length,
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
}

workflow auto {
    // 1. Pre-screen drivers
    ch_drivers_db = Channel.fromPath(params.drivers_db)
    ch_prescreen = prescreen_drivers(
        params.target_f_low,
        params.target_f_high,
        params.mouth_radius,
        params.length,
        ch_drivers_db,
    )

    // 2. Generate geometry for each profile using throat radius from prescreen
    ch_profiles = Channel.from("conical", "exponential", "hyperbolic")

    // Extract throat radius from prescreen result
    ch_throat_radius = ch_prescreen.map { json_file ->
        def data = new groovy.json.JsonSlurper().parse(json_file)
        data.throat_radius_m
    }

    // Combine profiles with throat radius
    ch_geom_inputs = ch_profiles.combine(ch_throat_radius)

    // 3. Generate geometries (3 parallel jobs)
    ch_geometries = generate_auto_geometry(ch_geom_inputs)

    // 4. Create band indices and combine with geometries
    ch_band_indices = Channel.from(0..<params.num_bands)
    ch_sim_inputs = ch_geometries.combine(ch_band_indices)

    // 5. Run simulations (3 profiles x num_bands parallel jobs)
    ch_band_results = run_auto_simulation(ch_sim_inputs)

    // 6. Group by profile and merge
    ch_grouped = ch_band_results.groupTuple()
    ch_merged = merge_auto_results(ch_grouped)

    // 7. Score and rank all combinations
    ch_all_csvs = ch_merged.map { profile, csv -> csv }.collect()
    ch_ranked = score_and_rank(
        ch_all_csvs,
        ch_prescreen,
        ch_drivers_db,
    )

    // 8. Generate report
    generate_auto_report(
        ch_ranked,
        ch_all_csvs,
        ch_drivers_db,
        ch_prescreen,
    )
}

workflow fullauto {
    // 1. Pre-screen drivers (no mouth/length needed)
    ch_drivers_db = Channel.fromPath(params.drivers_db)
    ch_prescreen = prescreen_drivers(
        params.target_f_low,
        params.target_f_high,
        0,  // placeholder — not used by prescreen filtering
        0,  // placeholder — not used by prescreen filtering
        ch_drivers_db,
    )

    // 2. Derive geometry from frequency band + prescreen throat radius
    ch_geom_derived = derive_fullauto_geometry(ch_prescreen)
    ch_candidates_csv = ch_geom_derived.map { csv, json -> csv }
    ch_design_json = ch_geom_derived.map { csv, json -> json }

    // 3. Parse candidates CSV into channel of tuples
    ch_candidates = ch_candidates_csv
        .splitCsv(header: true)
        .map { row ->
            tuple(row.candidate_id, row.profile, row.throat_radius as double,
                  row.mouth_radius as double, row.length as double)
        }

    // 4. Generate STEP geometry for each candidate
    ch_geometries = generate_fullauto_geometry(ch_candidates)

    // 5. Read sim freq range from design.json and combine with band indices
    ch_sim_range = ch_design_json.map { json_file ->
        def data = new groovy.json.JsonSlurper().parse(json_file)
        tuple(data.sim_freq_range[0] as double, data.sim_freq_range[1] as double)
    }

    ch_band_indices = Channel.from(0..<params.num_bands)
    ch_sim_inputs = ch_geometries
        .combine(ch_band_indices)
        .combine(ch_sim_range)

    // 6. Run FEM simulations (candidates x bands)
    ch_band_results = run_fullauto_simulation(ch_sim_inputs)

    // 7. Group by candidate_id and merge bands
    ch_grouped = ch_band_results.groupTuple()
    ch_merged = merge_fullauto_results(ch_grouped)

    // 8. Score and rank all driver-horn combinations
    ch_all_csvs = ch_merged.map { candidate_id, csv -> csv }.collect()
    ch_ranked = score_and_rank_fullauto(
        ch_all_csvs,
        ch_prescreen,
        ch_drivers_db,
        ch_candidates_csv,
    )

    // 9. Generate report with design summary
    generate_fullauto_report(
        ch_ranked,
        ch_all_csvs,
        ch_drivers_db,
        ch_prescreen,
        ch_design_json,
    )
}

workflow {
    if (params.mode == "fullauto") {
        fullauto()
    } else if (params.mode == "auto") {
        auto()
    } else {
        single()
    }
}
