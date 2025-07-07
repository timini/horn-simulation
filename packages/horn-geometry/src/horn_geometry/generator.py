import argparse
import FreeCAD
import Part
from pathlib import Path

from horn_core.parameters import HornParameters, FlareProfile

def create_conical_horn(params: HornParameters):
    """
    Creates a conical horn solid geometry using FreeCAD.
    """
    cone = Part.makeCone(params.throat_radius, params.mouth_radius, params.length)
    return cone

def main():
    parser = argparse.ArgumentParser(description="Generate a horn geometry STEP file.")
    parser.add_argument("--profile", required=True, choices=["conical", "exponential", "tractrix"])
    parser.add_argument("--throat", type=float, required=True, help="Throat radius in meters.")
    parser.add_argument("--mouth", type=float, required=True, help="Mouth radius in meters.")
    parser.add_argument("--length", type=float, required=True, help="Length in meters.")
    parser.add_argument("--output-file", type=str, required=True, help="Path to save the STEP file.")
    args = parser.parse_args()

    params = HornParameters(
        flare_profile=FlareProfile(args.profile),
        throat_radius=args.throat,
        mouth_radius=args.mouth,
        length=args.length,
    )

    if params.flare_profile == FlareProfile.CONICAL:
        horn_solid = create_conical_horn(params)
    else:
        raise NotImplementedError(f"Profile '{params.flare_profile.value}' is not yet implemented.")

    output_path = Path(args.output_file)
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        horn_solid.exportStep(str(output_path))
        print(f"Successfully exported horn geometry to: {output_path}")
    except Exception as e:
        print(f"Error exporting STEP file: {e}")
        print(f"Is horn solid valid? {horn_solid.isValid()}")
        raise

if __name__ == "__main__":
    main()
