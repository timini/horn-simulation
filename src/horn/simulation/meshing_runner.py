import argparse
from pathlib import Path
from horn.simulation import meshing

def main():
    """
    Command-line entry point for the meshing module.
    
    This script is executed inside the solver Docker container.
    """
    parser = argparse.ArgumentParser(description="Generate a mesh from a STEP file.")
    parser.add_argument("--step-file", type=str, required=True, help="Path to the input STEP file.")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to save the mesh file.")
    parser.add_argument("--max-freq", type=float, required=True, help="Maximum frequency for mesh sizing.")
    args = parser.parse_args()

    print("--- Running Meshing Script ---")
    mesh_file = meshing.create_mesh(
        step_file=args.step_file,
        max_freq=args.max_freq,
        output_dir=args.output_dir
    )
    print(f"Successfully created mesh file: {mesh_file}")
    print("--- Meshing Script Complete ---")

if __name__ == "__main__":
    main() 