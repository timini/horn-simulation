import argparse
from pathlib import Path
import json
from horn.simulation import solver

def main():
    """
    Command-line entry point for the solver module.
    
    This script is executed inside the solver Docker container.
    """
    parser = argparse.ArgumentParser(description="Run the FEM simulation.")
    parser.add_argument("--mesh-file", type=str, required=True, help="Path to the input mesh file.")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to save the results file.")
    parser.add_argument("--freq-min", type=float, required=True, help="Minimum frequency.")
    parser.add_argument("--freq-max", type=float, required=True, help="Maximum frequency.")
    parser.add_argument("--driver-params-json", type=str, required=True, help="JSON string of driver params.")
    args = parser.parse_args()

    driver_params = json.loads(args.driver_params_json)
    freq_range = (args.freq_min, args.freq_max)

    print("--- Running Solver Script ---")
    results_file = solver.run_simulation(
        mesh_file=args.mesh_file,
        freq_range=freq_range,
        driver_params=driver_params,
        output_dir=Path(args.output_dir)
    )
    print(f"Successfully created results file: {results_file}")
    print("--- Solver Script Complete ---")

if __name__ == "__main__":
    main() 