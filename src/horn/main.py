import os
import tempfile
from horn.data_ingestion import driver_db
from horn.geometry import horn_generator
from horn.geometry.parameters import HornParameters
from horn.simulation import meshing, solver
from horn.analysis import analyzer
from typing import Tuple, Dict, Any

def run_pipeline(driver_id: str, horn_params: HornParameters, freq_range_hz: Tuple[float, float]) -> Dict[str, Any]:
    """
    Executes the full horn simulation and analysis pipeline.

    Args:
        driver_id: The ID of the driver to use for the simulation.
        horn_params: The geometric parameters of the horn to simulate.
        freq_range_hz: A tuple containing the start and end frequencies in Hz.

    Returns:
        A dictionary containing the final analysis report, including
        metrics, scores, and paths to generated plots.
    """
    start_freq, end_freq = freq_range_hz
    
    # Create a temporary directory for this pipeline run
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Created temporary directory for pipeline run: {tmpdir}")

        # Stage 1: Data Ingestion
        print(f"\n--- Stage 1: Data Ingestion ---")
        driver_params = driver_db.get_driver_parameters(driver_id)
        print(f"Successfully fetched driver data: {driver_params['manufacturer']} {driver_params['model_name']}")

        # Stage 2: Parametric Geometry Generation
        print(f"\n--- Stage 2: Geometry Generation ---")
        step_file_path = horn_generator.create_horn(horn_params, output_dir=tmpdir)
        print(f"Successfully generated geometry file: {step_file_path}")

        # Stage 3.1: Meshing
        print(f"\n--- Stage 3.1: Meshing ---")
        # The mesh must be fine enough for the highest frequency
        mesh_file_path = meshing.create_mesh(step_file_path, end_freq)
        print(f"Successfully generated mesh file: {mesh_file_path}")
        
        # Stage 3.2: Solving
        print(f"\n--- Stage 3.2: Simulation Solver ---")
        results_path = solver.run_simulation(mesh_file_path, driver_params, freq_range_hz)
        print(f"Successfully generated results file: {results_path}")

        # Stage 4: Analysis
        print(f"\n--- Stage 4: Analysis & Scoring ---")
        analysis_report = analyzer.analyze_results(results_path)
        print(f"Successfully generated analysis report.")

        print(f"\n--- Pipeline Complete ---")
        # Note: The temp directory and its contents are deleted upon exiting the 'with' block.
        # In a real application, we would copy the final report to a persistent location.
        return analysis_report
