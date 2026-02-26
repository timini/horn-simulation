#!/usr/bin/env nextflow

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
params.test_outdir = null  // for testing purposes

// Auto-select mode settings
params.mode = "single"              // "single" or "auto_select"
params.target_f_low = 500           // Target low frequency (Hz)
params.target_f_high = 4000         // Target high frequency (Hz)
params.max_length = 0.5             // Max horn length (m)
params.max_mouth_radius = 0.2       // Max mouth radius (m)
params.drivers_db = "data/drivers.json"
params.top_n = 5                    // Number of top candidates for Phase B

// ─────────────────────────────────────────────────────────────────────────────
// SHARED PROCESSES (used by both single and auto_select modes)
// ─────────────────────────────────────────────────────────────────────────────

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

// ─────────────────────────────────────────────────────────────────────────────
// AUTO-SELECT PROCESSES
// ─────────────────────────────────────────────────────────────────────────────

process generate_candidates {
    publishDir "${params.outdir}/auto_select", mode: 'copy'

    output:
    path "candidates.csv"

    script:
    def profiles_arg = ""
    """
    python3 -m horn_core.candidates \
        --target-f-low ${params.target_f_low} \
        --target-f-high ${params.target_f_high} \
        --max-length ${params.max_length} \
        --max-mouth-radius ${params.max_mouth_radius} \
        --output candidates.csv
    """
}

process generate_candidate_geometry {
    input:
    tuple val(candidate_id), val(profile), val(throat_radius), val(mouth_radius), val(length)

    output:
    tuple val(candidate_id), val(throat_radius), val(length), path("${candidate_id}.step")

    script:
    """
    python3 -m horn_geometry.generator \
        --throat-radius ${throat_radius} \
        --mouth-radius ${mouth_radius} \
        --length ${length} \
        --profile ${profile} \
        --num-sections 20 \
        --output-file ${candidate_id}.step
    """
}

process run_candidate_simulation {
    input:
    tuple val(candidate_id), val(throat_radius), val(length), path(step_file)

    output:
    tuple val(candidate_id), val(throat_radius), path("${candidate_id}_results.csv")

    script:
    """
    python3 -m horn_solver.solver \
        --step-file ${step_file} \
        --output-file ${candidate_id}_results.csv \
        --min-freq ${params.target_f_low} \
        --max-freq ${params.target_f_high} \
        --num-intervals ${params.num_intervals} \
        --length ${length} \
        --mesh-size ${params.mesh_size}
    """
}

process screen_drivers {
    publishDir "${params.outdir}/auto_select/screening", mode: 'copy'

    input:
    tuple val(candidate_id), val(throat_radius), path(solver_csv)

    output:
    tuple val(candidate_id), path("${candidate_id}_screening.csv")

    script:
    """
    python3 -m horn_analysis.transfer_function \
        --solver-csv ${solver_csv} \
        --drivers-db ${params.drivers_db} \
        --throat-radius ${throat_radius} \
        --output-dir . && \
    mv screening_summary.csv ${candidate_id}_screening.csv
    """
}

process rank_candidates {
    publishDir "${params.outdir}/auto_select", mode: 'copy'

    input:
    path(screening_csvs)
    path(solver_csvs)

    output:
    path "ranked_candidates.json"
    path "top_candidates.csv"

    script:
    """
    python3 << 'PYEOF'
import json, glob
import pandas as pd
from pathlib import Path
from horn_analysis.kpi import extract_kpis
from horn_analysis.scoring import TargetSpec, compute_selection_score, rank_candidates

target = TargetSpec(f_low_hz=${params.target_f_low}, f_high_hz=${params.target_f_high})

# Collect all screening results
screening_files = sorted(glob.glob("*_screening.csv"))
solver_files = sorted(glob.glob("*_results.csv"))

all_scores = []
for sf in solver_files:
    candidate_id = Path(sf).stem.replace("_results", "")
    kpi = extract_kpis(sf)

    # Find matching screening file
    screen_file = f"{candidate_id}_screening.csv"
    if Path(screen_file).exists():
        screen_df = pd.read_csv(screen_file)
        for _, row in screen_df.iterrows():
            score = compute_selection_score(
                kpi, target,
                driver_id=row["driver_id"],
                horn_label=candidate_id,
            )
            all_scores.append(score)

top = rank_candidates(all_scores, top_n=${params.top_n})

# Write ranked JSON
ranked_data = [s.to_dict() for s in top]
Path("ranked_candidates.json").write_text(json.dumps(ranked_data, indent=2))

# Write CSV for Phase B
rows = [{"candidate_id": s.horn_label, "driver_id": s.driver_id,
         "composite_score": s.composite_score} for s in top]
pd.DataFrame(rows).to_csv("top_candidates.csv", index=False)

print(f"Ranked {len(all_scores)} combinations, selected top {len(top)}")
for s in top:
    print(f"  {s.horn_label} + {s.driver_id}: score={s.composite_score:.3f}")
PYEOF
    """
}

process run_neumann_simulation {
    publishDir "${params.outdir}/auto_select/phase_b", mode: 'copy'

    input:
    tuple val(candidate_id), val(driver_id), val(throat_radius), val(length), path(step_file), path(phase_a_csv)

    output:
    tuple val(candidate_id), val(driver_id), path("${candidate_id}_${driver_id}_neumann.csv")

    script:
    """
    python3 -m horn_solver.solver \
        --step-file ${step_file} \
        --output-file ${candidate_id}_${driver_id}_neumann.csv \
        --min-freq ${params.target_f_low} \
        --max-freq ${params.target_f_high} \
        --num-intervals ${params.num_intervals} \
        --length ${length} \
        --mesh-size ${params.mesh_size} \
        --bc-mode neumann \
        --driver-json ${params.drivers_db} \
        --driver-id ${driver_id} \
        --throat-area \$(python3 -c "import math; print(math.pi * ${throat_radius}**2)") \
        --phase-a-csv ${phase_a_csv}
    """
}

process generate_selection_report {
    publishDir "${params.outdir}/auto_select", mode: 'copy'

    input:
    path(neumann_csvs)
    path(ranked_json)

    output:
    path "selection_report.png"
    path "selection_summary.json"

    script:
    """
    python3 << 'PYEOF'
import json, glob
from pathlib import Path
from horn_analysis.compare import plot_multi_comparison
from horn_analysis.kpi import extract_kpis

# Collect Phase B results
neumann_files = sorted(glob.glob("*_neumann.csv"))
file_label_pairs = []
summary = []

for f in neumann_files:
    name = Path(f).stem.replace("_neumann", "")
    file_label_pairs.append((f, name))
    kpi = extract_kpis(f)
    summary.append({"label": name, "kpis": kpi.to_dict()})

if file_label_pairs:
    plot_multi_comparison(file_label_pairs, "selection_report.png", kpi_table=True)

# Combine with ranking data
ranked = json.loads(Path("ranked_candidates.json").read_text())
report = {"ranked_candidates": ranked, "phase_b_results": summary}
Path("selection_summary.json").write_text(json.dumps(report, indent=2))

print(f"Selection report generated with {len(neumann_files)} Phase B results")
PYEOF
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW DEFINITIONS
// ─────────────────────────────────────────────────────────────────────────────

workflow single_horn {
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

workflow auto_select {
    // ── Phase A: Generate candidates and run Dirichlet simulations ──

    // 1. Generate candidate geometries
    ch_candidates_csv = generate_candidates()

    // 2. Parse CSV into per-candidate channels
    ch_candidates = ch_candidates_csv
        .splitCsv(header: true)
        .map { row -> tuple(
            row.candidate_id, row.profile,
            row.throat_radius, row.mouth_radius, row.length
        )}

    // 3. Generate STEP geometry for each candidate
    ch_candidate_steps = generate_candidate_geometry(ch_candidates)

    // 4. Run Phase A (Dirichlet) simulation for each candidate
    ch_candidate_results = run_candidate_simulation(ch_candidate_steps)

    // 5. Screen all drivers against each candidate geometry
    ch_screening = screen_drivers(ch_candidate_results)

    // 6. Rank all driver-horn combinations
    ch_ranked = rank_candidates(
        ch_screening.map { it[1] }.collect(),
        ch_candidate_results.map { it[2] }.collect()
    )

    // ── Phase B: Re-simulate top candidates with Neumann BC ──

    // 7. Parse top candidates and join with geometry/solver data
    ch_top = ch_ranked[1]
        .splitCsv(header: true)
        .map { row -> tuple(row.candidate_id, row.driver_id) }

    // Join top candidates with their step files and Phase A CSV
    // ch_candidate_steps emits: (candidate_id, throat_radius, length, step_file)
    // ch_candidate_results emits: (candidate_id, throat_radius, results_csv)
    ch_step_lookup = ch_candidate_steps
        .map { tuple(it[0], it[1], it[2], it[3]) }  // (id, throat_r, length, step)

    ch_csv_lookup = ch_candidate_results
        .map { tuple(it[0], it[2]) }  // (id, csv)

    ch_neumann_inputs = ch_top
        .combine(ch_step_lookup, by: 0)  // join on candidate_id
        .combine(ch_csv_lookup, by: 0)   // join on candidate_id
        .map { tuple(it[0], it[1], it[2], it[3], it[4], it[5]) }
        // (candidate_id, driver_id, throat_radius, length, step_file, phase_a_csv)

    // 8. Run Neumann simulations for top candidates
    ch_neumann_results = run_neumann_simulation(ch_neumann_inputs)

    // 9. Generate final selection report
    generate_selection_report(
        ch_neumann_results.map { it[2] }.collect(),
        ch_ranked[0]
    )
}

workflow {
    if (params.mode == "auto_select") {
        auto_select()
    } else {
        single_horn()
    }
}
