import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def plot_frequency_response(input_file: Path, output_file: Path):
    """
    Reads simulation results from a CSV and plots the frequency response.
    """
    df = pd.read_csv(input_file)

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots()

    ax.plot(df['frequency'], df['spl'], marker='o', linestyle='-')

    ax.set_title('Horn Frequency Response')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Sound Pressure Level (dB)')
    ax.grid(True)
    
    # Ensure the y-axis is sensible, handling potential -inf values
    valid_spl = df.loc[df['spl'] > -float('inf'), 'spl']
    if not valid_spl.empty:
        ax.set_ylim(max(0, valid_spl.min() - 10), valid_spl.max() + 10)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Successfully generated plot at: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Plot frequency response from a CSV file.")
    parser.add_argument("--input", type=str, required=True, help="Input CSV file path.")
    parser.add_argument("--output", type=str, required=True, help="Output plot image path.")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    plot_frequency_response(input_path, output_path)

if __name__ == "__main__":
    main() 