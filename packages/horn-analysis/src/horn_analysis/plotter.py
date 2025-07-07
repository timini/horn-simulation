import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def plot_results(input_file: str, output_file: str):
    """
    Reads simulation results from a CSV and generates a frequency response plot.
    """
    df = pd.read_csv(input_file)
    plt.figure()
    plt.plot(df["frequency"], df["spl"])
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("SPL (dB)")
    plt.title("Frequency Response")
    plt.grid(True)
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Plot horn simulation results.")
    parser.add_argument("--input", type=str, required=True, help="Input CSV file.")
    parser.add_argument("--output", type=str, required=True, help="Output PNG file.")
    args = parser.parse_args()
    
    # Ensure output directory exists
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    plot_results(args.input, args.output)

if __name__ == "__main__":
    main() 