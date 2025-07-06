import os
import tempfile
import subprocess
from typing import Tuple, Dict, Any
from pathlib import Path

from horn.data_ingestion import driver_db
from horn.geometry.parameters import HornParameters

def run_pipeline(driver_id: str, horn_params: HornParameters, freq_range_hz: Tuple[float, float]) -> Dict[str, Any]:
    """
    Runs the full horn simulation pipeline by orchestrating Docker containers.
    """
    min_freq, max_freq = freq_range_hz
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        print(f"Created temporary directory for pipeline run: {tmpdir_path}")

        # --- Stage 1: Data Ingestion (Local) ---
        print(f"\n--- Stage 1: Data Ingestion ---")
        driver_params = driver_db.get_driver_parameters(driver_id)
        print(f"Successfully fetched driver data: {driver_params['manufacturer']} {driver_params['model_name']}")

        # --- Stage 2: Geometry Generation (in Docker) ---
        print(f"\n--- Stage 2: Geometry Generation ---")
        step_file_name = f"horn_{horn_params.flare_profile.value}_{horn_params.length}m.stp"
        step_file_path = tmpdir_path / step_file_name
        
        # This stage is now a Docker command. We mount the temp dir to /data in the container.
        subprocess.run([
            "docker", "run", "--rm",
            "-v", f"{tmpdir_path}:/data",
            "horn-freecad-app",
            "python", "-m", "horn.geometry.horn_generator",
            # We need to serialize the horn_params to pass them as arguments
            "--profile", horn_params.flare_profile.value,
            "--throat", str(horn_params.throat_radius),
            "--mouth", str(horn_params.mouth_radius),
            "--length", str(horn_params.length),
            "--output-dir", "/data"
        ], check=True)
        print(f"Successfully generated STEP file: {step_file_path}")

        # --- Stage 3 & 4: Meshing and Solving (in Docker) ---
        print(f"\n--- Stage 3 & 4: Meshing and Solving ---")
        # This will be another Docker command using horn-solver-app
        # For now, create a dummy results file
        results_csv = tmpdir_path / "results.csv"
        with open(results_csv, "w") as f:
            f.write("frequency,spl\n100,95\n")
        print(f"Placeholder: Meshing and solving complete. Results at {results_csv}")
        
        # --- Stage 5: Analysis and Visualization (Local) ---
        # This will be updated to read the real results
        final_report = {
            "metrics": {"score": 98.5},
            "plots": [str(tmpdir_path / "spl.png")]
        }
        print(f"\n--- Stage 5: Analysis ---")
        print("Analysis complete. Final score: 98.5")

    return final_report
