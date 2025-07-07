import gmsh
from pathlib import Path
import argparse

def create_conical_horn(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
) -> Path:
    """
    Generates a simple conical horn geometry using gmsh's OpenCASCADE kernel
    and saves it as a STEP file.

    Args:
        throat_radius: The radius of the horn's throat (inlet).
        mouth_radius: The radius of the horn's mouth (outlet).
        length: The length of the horn along the z-axis.
        output_file: The exact path to save the output STEP file to.

    Returns:
        The path to the generated STEP file.
    """
    gmsh.initialize()
    gmsh.model.add("conical_horn")
    gmsh.option.setNumber("General.Terminal", 1)

    # Create the throat and mouth circles (wires)
    throat_wire = gmsh.model.occ.addCircle(0, 0, 0, throat_radius)
    mouth_wire = gmsh.model.occ.addCircle(0, 0, length, mouth_radius)

    # Create curve loops from the wires
    throat_loop = gmsh.model.occ.addCurveLoop([throat_wire])
    mouth_loop = gmsh.model.occ.addCurveLoop([mouth_wire])

    # Create plane surfaces from the loops
    gmsh.model.occ.addPlaneSurface([throat_loop])
    gmsh.model.occ.addPlaneSurface([mouth_loop])

    # Create the lofted solid (conical section)
    gmsh.model.occ.addThruSections([throat_loop, mouth_loop])

    gmsh.model.occ.synchronize()

    # Export the model to a STEP file
    gmsh.write(str(output_file))

    gmsh.finalize()

    # Print the final path to stdout so it can be captured by the calling script.
    print(str(output_file))
    return output_file

def main():
    """A simple command-line interface for the generator."""
    parser = argparse.ArgumentParser(description="Generate a conical horn STEP file.")
    parser.add_argument("--throat-radius", type=float, required=True)
    parser.add_argument("--mouth-radius", type=float, required=True)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--output-file", type=str, required=True)
    args = parser.parse_args()

    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    create_conical_horn(
        throat_radius=args.throat_radius,
        mouth_radius=args.mouth_radius,
        length=args.length,
        output_file=output_path,
    )

if __name__ == "__main__":
    main()
