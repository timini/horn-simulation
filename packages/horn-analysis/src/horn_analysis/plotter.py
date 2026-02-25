import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

def plot_spl_vs_frequency(csv_file: str, output_image_file: str):
    """
    Reads simulation results from a CSV and plots SPL vs. Frequency.

    Args:
        csv_file: Path to the input CSV file.
        output_image_file: Path to save the output plot image.
    """
    # Read the data using pandas
    data = pd.read_csv(csv_file)

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(data['frequency'], data['spl'], marker='o', linestyle='-')
    plt.grid(True)
    plt.title('Sound Pressure Level (SPL) vs. Frequency')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('SPL (dB)')
    plt.xscale('log') # Frequency is often better viewed on a log scale
    
    # Save the plot to a file
    plt.savefig(output_image_file)
    print(f"Plot saved to {output_image_file}")

def main():
    """Command-line interface for the plotting script."""
    parser = argparse.ArgumentParser(description="Plot SPL vs. Frequency from a CSV file.")
    parser.add_argument("csv_file", type=str, help="Path to the input CSV file.")
    parser.add_argument("output_image_file", type=str, help="Path to save the output plot image.")
    args = parser.parse_args()
    
    plot_spl_vs_frequency(args.csv_file, args.output_image_file)

if __name__ == "__main__":
    main()