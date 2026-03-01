import pandas as pd
import matplotlib
matplotlib.use('Agg')
import argparse

from horn_analysis import plot_theme


def plot_spl_vs_frequency(csv_file: str, output_image_file: str):
    """
    Reads simulation results from a CSV and plots SPL vs. Frequency.

    Args:
        csv_file: Path to the input CSV file.
        output_image_file: Path to save the output plot image.
    """
    # Read the data using pandas
    data = pd.read_csv(csv_file)

    fig, ax = plot_theme.create_figure(figsize=(10, 6))
    ax.plot(data['frequency'], data['spl'], color=plot_theme.COLORS["primary"], linewidth=1.4)
    ax.set_title('Sound Pressure Level (SPL) vs. Frequency')

    plot_theme.setup_freq_axis(ax, data['frequency'].min(), data['frequency'].max())
    plot_theme.setup_spl_axis(ax, data['spl'].values)
    plot_theme.setup_grid(ax)

    plot_theme.save_figure(fig, output_image_file)
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
