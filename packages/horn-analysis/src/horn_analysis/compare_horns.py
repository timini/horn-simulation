import pandas as pd
import matplotlib
matplotlib.use("Agg")
import sys

from horn_analysis import plot_theme


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

    fig, ax = plot_theme.create_figure(figsize=(12, 8))
    ax.plot(df_a['frequency'], df_a['spl'], label=label_a,
            color=plot_theme.MULTI_COLORS[0], linewidth=1.4)
    ax.plot(df_b['frequency'], df_b['spl'], label=label_b,
            color=plot_theme.MULTI_COLORS[1], linewidth=1.4)

    all_freq = list(df_a['frequency'].values) + list(df_b['frequency'].values)
    plot_theme.setup_freq_axis(ax, min(all_freq), max(all_freq))
    plot_theme.setup_grid(ax)

    ax.set_title('Horn Frequency Response Comparison')
    ax.set_ylabel('Sound Pressure Level (dB)')
    ax.legend()

    plot_theme.save_figure(fig, output_file)
    print(f"Comparison plot saved to {output_file}")

if __name__ == "__main__":
    # Expects: python compare_horns.py <file_a> <label_a> <file_b> <label_b> <output_file>
    if len(sys.argv) != 6:
        print("Usage: python compare_horns.py <file_a> <label_a> <file_b> <label_b> <output_file>")
        sys.exit(1)

    plot_comparison(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
