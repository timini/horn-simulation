import argparse
import pandas as pd
from pathlib import Path

def analyze_results(results_file: str) -> dict:
    """
    Performs a simple analysis of the simulation results.
    """
    results_df = pd.read_csv(results_file)
    
    avg_spl = results_df['spl'].mean()
    
    report = {
        "source_file": str(results_file),
        "metrics": {
            "average_spl_db": round(avg_spl, 2)
        }
    }
    return report

def main():
    """
    Entry point for running the analysis.
    """
    parser = argparse.ArgumentParser(description="Analyze horn simulation results.")
    parser.add_argument("--results-file", required=True, help="Path to the results CSV file.")
    args = parser.parse_args()

    report = analyze_results(args.results_file)
    
    print("\n--- Analysis Report ---")
    print(f"Source: {report['source_file']}")
    print(f"Average SPL: {report['metrics']['average_spl_db']} dB")
    print("-----------------------\n")

if __name__ == "__main__":
    main() 