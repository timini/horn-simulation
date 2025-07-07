import argparse
import json
from pathlib import Path
from horn_solver import solver

def main():
    """
    Command-line entry point for the solver module.
    
    This script is executed inside the solver Docker container.
    """
    print("--- Running Solver Script ---")
    parser = argparse.ArgumentParser(description="Run the acoustic simulation.")
    parser.add_argument("--step-file", required=True, help="Path to the input STEP file.")
    parser.add_argument("--output-file", required=True, help="Path to save the results CSV file.")
    parser.add_argument("--freq-min", type=float, required=True, help="Minimum frequency.")
    parser.add_argument("--freq-max", type=float, required=True, help="Maximum frequency.")
    parser.add_argument("--driver-params-json", required=True, help="JSON string of driver parameters.")
    args = parser.parse_args()

    driver_params = json.loads(args.driver_params_json)
    freq_range = (args.freq_min, args.freq_max)

    results_file = solver.run_simulation_from_step(
        step_file=args.step_file,
        freq_range=freq_range,
        driver_params=driver_params,
        output_file=args.output_file,
        max_freq_mesh=args.freq_max, # Use max sim freq for meshing for now
    )
    print(f"Successfully wrote results to {results_file}")
    print("--- Solver Script Complete ---")

if __name__ == "__main__":
    main() 