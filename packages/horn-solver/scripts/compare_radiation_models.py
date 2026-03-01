#!/usr/bin/env python3
"""Compare SPL curves from different radiation impedance models.

Runs the horn solver with plane_wave, flanged_piston, and bem radiation
models on the same STEP geometry, then plots an overlay comparison.

Usage:
    python compare_radiation_models.py \
        --step-file horn.step \
        --length 0.5 \
        --min-freq 200 \
        --max-freq 4000 \
        --num-intervals 50 \
        --output-dir comparison_results
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


RADIATION_MODELS = ["plane_wave", "flanged_piston", "bem"]


def run_model(
    step_file: str,
    length: float,
    freq_range: tuple,
    num_intervals: int,
    mesh_size: float,
    radiation_model: str,
    output_dir: Path,
) -> pd.DataFrame:
    """Run the solver for a single radiation model and return results."""
    from horn_solver.solver import run_simulation_from_step

    output_file = output_dir / f"results_{radiation_model}.csv"
    driver_params = {"length": length}

    try:
        run_simulation_from_step(
            step_file=step_file,
            freq_range=freq_range,
            num_intervals=num_intervals,
            driver_params=driver_params,
            output_file=str(output_file),
            max_freq_mesh=freq_range[1],
            mesh_size=mesh_size,
            radiation_model=radiation_model,
        )
        return pd.read_csv(output_file)
    except (ImportError, RuntimeError) as exc:
        print(f"WARNING: {radiation_model} failed: {exc}")
        return None


def plot_comparison(all_results: dict, output_dir: Path):
    """Plot overlaid SPL curves for all models."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {"plane_wave": "blue", "flanged_piston": "orange", "bem": "red"}
    labels = {
        "plane_wave": "Plane Wave (Z=rho*c)",
        "flanged_piston": "Flanged Piston",
        "bem": "BEM (nonlocal)",
    }

    for model, df in all_results.items():
        ax.semilogx(
            df["frequency"], df["spl"],
            color=colors.get(model, "gray"),
            label=labels.get(model, model),
            linewidth=1.5,
        )

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("SPL (dB)")
    ax.set_title("Radiation Model Comparison")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plot_path = output_dir / "radiation_model_comparison.png"
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved: {plot_path}")


def compute_deviations(all_results: dict, output_dir: Path):
    """Compute pairwise deviations between models and save to CSV."""
    models = list(all_results.keys())
    rows = []

    for i in range(len(models)):
        for j in range(i + 1, len(models)):
            m1, m2 = models[i], models[j]
            df1, df2 = all_results[m1], all_results[m2]

            # Interpolate to common frequencies if needed
            if np.allclose(df1["frequency"].values, df2["frequency"].values):
                diff = np.abs(df1["spl"].values - df2["spl"].values)
            else:
                from scipy.interpolate import interp1d
                f_common = np.union1d(df1["frequency"].values, df2["frequency"].values)
                spl1 = interp1d(df1["frequency"], df1["spl"], fill_value="extrapolate")(f_common)
                spl2 = interp1d(df2["frequency"], df2["spl"], fill_value="extrapolate")(f_common)
                diff = np.abs(spl1 - spl2)

            rows.append({
                "model_1": m1,
                "model_2": m2,
                "max_deviation_dB": float(np.max(diff)),
                "mean_deviation_dB": float(np.mean(diff)),
                "std_deviation_dB": float(np.std(diff)),
            })
            print(f"  {m1} vs {m2}: max={np.max(diff):.2f} dB, mean={np.mean(diff):.2f} dB")

    dev_df = pd.DataFrame(rows)
    dev_path = output_dir / "model_deviations.csv"
    dev_df.to_csv(dev_path, index=False)
    print(f"Deviations saved: {dev_path}")


def main():
    parser = argparse.ArgumentParser(description="Compare radiation impedance models.")
    parser.add_argument("--step-file", type=str, required=True)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--min-freq", type=float, default=200.0)
    parser.add_argument("--max-freq", type=float, default=4000.0)
    parser.add_argument("--num-intervals", type=int, default=50)
    parser.add_argument("--mesh-size", type=float, default=0.01)
    parser.add_argument("--output-dir", type=str, default="comparison_results")
    parser.add_argument(
        "--models", type=str, nargs="+", default=RADIATION_MODELS,
        help="Radiation models to compare (default: all)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    freq_range = (args.min_freq, args.max_freq)
    all_results = {}

    for model in args.models:
        print(f"\n{'='*60}")
        print(f"Running: {model}")
        print(f"{'='*60}")

        df = run_model(
            step_file=args.step_file,
            length=args.length,
            freq_range=freq_range,
            num_intervals=args.num_intervals,
            mesh_size=args.mesh_size,
            radiation_model=model,
            output_dir=output_dir,
        )
        if df is not None:
            all_results[model] = df

    if len(all_results) < 2:
        print("ERROR: Need at least 2 successful models to compare")
        return 1

    print(f"\n{'='*60}")
    print("Computing deviations...")
    print(f"{'='*60}")
    compute_deviations(all_results, output_dir)

    print(f"\n{'='*60}")
    print("Generating comparison plot...")
    print(f"{'='*60}")
    plot_comparison(all_results, output_dir)

    print(f"\nDone. Results in: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
