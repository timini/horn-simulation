import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_comparison(file_a, label_a, file_b, label_b, output_file):
    """
    Reads two CSV files containing frequency response data and plots them on the same graph.

    Args:
        file_a (str): Path to the first CSV file.
        label_a (str): Label for the first horn.
        file_b (str): Path to the second CSV file.
        label_b (str): Label for the second horn.
        output_file (str): Path to save the output plot image.
    """
    # Read the data
    df_a = pd.read_csv(file_a)
    df_b = pd.read_csv(file_b)

    # Create the plot
    plt.figure(figsize=(12, 8))
    plt.plot(df_a['frequency'], df_a['spl'], label=label_a)
    plt.plot(df_b['frequency'], df_b['spl'], label=label_b)

    # Formatting
    plt.xscale('log')
    plt.title('Horn Frequency Response Comparison')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Sound Pressure Level (dB)')
    plt.grid(True, which="both", ls="--")
    plt.legend()

    # Save the plot
    plt.savefig(output_file)
    print(f"Comparison plot saved to {output_file}")

if __name__ == "__main__":
    # Expects: python compare_horns.py <file_a> <label_a> <file_b> <label_b> <output_file>
    if len(sys.argv) != 6:
        print("Usage: python compare_horns.py <file_a> <label_a> <file_b> <label_b> <output_file>")
        sys.exit(1)
    
    plot_comparison(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]) 