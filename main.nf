#!/usr/bin/env nextflow

// Horn Geometry
params.throat_radius = 0.05 // Radius of the horn's throat in meters
params.mouth_radius = 0.2  // Radius of the horn's mouth in meters
params.length = 0.5        // Length of the horn in meters

// Simulation Settings
params.min_freq = 500      // Minimum frequency for the sweep in Hz
params.max_freq = 8000     // Maximum frequency for the sweep in Hz
params.num_intervals = 100 // Number of frequency steps in the sweep
params.mesh_size = 0.01    // Target mesh element size in meters

// Execution Settings
params.num_bands = 8       // Number of parallel jobs for the solver
params.outdir = "./results"
params.test_outdir = null  // for testing purposes

process generate_geometry {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val throat_radius
    val mouth_radius
    val length

    output:
    path "horn.step"

    script:
    """
    python3 -m horn_geometry.generator \
        --throat-radius ${throat_radius} \
        --mouth-radius ${mouth_radius} \
        --length ${length} \
        --output-file horn.step
    """
}

process run_simulation {
    input:
    tuple path(horn_step), val(band_index)

    output:
    path "results_${band_index}.csv"

    script:
    def band_width = (params.max_freq - params.min_freq) / params.num_bands
    def min_f = params.min_freq + band_width * band_index
    def max_f = params.min_freq + band_width * (band_index + 1)
    def num_intervals_per_band = params.num_intervals / params.num_bands
    """
    echo "Running band ${band_index}: ${min_f} Hz to ${max_f} Hz"
    python3 -m horn_solver.solver \
        --step-file ${horn_step} \
        --output-file results_${band_index}.csv \
        --min-freq ${min_f} \
        --max-freq ${max_f} \
        --num-intervals ${num_intervals_per_band as int} \
        --length ${params.length}
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

workflow {
    // 1. Generate geometry once
    ch_step_file = generate_geometry(
        params.throat_radius,
        params.mouth_radius,
        params.length
    )

    // 2. Create a channel of band indices
    ch_band_indices = Channel.from(0..<params.num_bands)

    // 3. Combine the geometry file with the band indices
    ch_sim_inputs = ch_step_file.combine(ch_band_indices)

    // 4. Run simulations in parallel
    ch_band_results = run_simulation(ch_sim_inputs)

    // 5. Collect and merge results
    ch_merged_results = merge_results(ch_band_results.collect())

    // 6. Plot final results
    generate_plots(ch_merged_results)
} 