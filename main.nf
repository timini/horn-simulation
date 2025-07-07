#!/usr/bin/env nextflow

params.throat_radius = 0.05
params.mouth_radius = 0.2
params.length = 0.5
params.min_freq = 100
params.max_freq = 2000
params.mesh_size = 0.01

process generate_geometry {
    output:
    path "horn.step"

    """
    python3 -m horn_geometry.generator \
        --throat-radius ${params.throat_radius} \
        --mouth-radius ${params.mouth_radius} \
        --length ${params.length} \
        --output-file horn.step
    """
}

process run_simulation {
    input:
    path step_file

    output:
    path "results.csv"

    """
    python3 -m horn_solver.solver \
        --step-file ${step_file} \
        --output-file results.csv \
        --min-freq ${params.min_freq} \
        --max-freq ${params.max_freq} \
        --mesh-size ${params.mesh_size}
    """
}

process generate_plots {
    input:
    path csv_file
    
    output:
    path "frequency_response.png"

    """
    python3 -m horn_analysis.plotter \
        --input ${csv_file} \
        --output frequency_response.png
    """
}

workflow {
    geom_ch = generate_geometry()
    sim_ch = run_simulation(geom_ch)
    generate_plots(sim_ch)
} 