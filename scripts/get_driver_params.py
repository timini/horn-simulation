import argparse
import json
from pathlib import Path

def get_driver_parameters(driver_id: str) -> dict:
    """
    Reads driver parameters from the JSON file.
    """
    # Assumes this script is in /scripts and data is in /data
    data_path = Path(__file__).parent.parent / "data/drivers.json"
    with open(data_path, 'r') as f:
        drivers = json.load(f)
    
    if driver_id not in drivers:
        raise ValueError(f"Driver ID '{driver_id}' not found in drivers.json")
    
    return drivers[driver_id]

def main():
    """
    Entry point to print driver parameters as a JSON string.
    """
    parser = argparse.ArgumentParser(description="Fetch driver parameters and print as JSON.")
    parser.add_argument("--driver-id", type=str, required=True, help="ID of the driver to fetch.")
    args = parser.parse_args()
    
    params = get_driver_parameters(args.driver_id)
    # Print the JSON string to stdout so it can be captured by the Makefile
    print(json.dumps(params))

if __name__ == "__main__":
    main() 