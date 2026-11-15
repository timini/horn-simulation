"""Electrical and excursion outputs for ranked driver-horn candidates.

Per candidate this gives what an amplifier, crossover and builder need from
the linear model: electrical impedance in ohms, peak cone excursion and
input power at the drive voltage, and on-axis system group delay with the
acoustic time of flight removed. An optional stress level (``test_power_w``)
is evaluated over the full simulated sweep, separately from the drive
voltage used for SPL and scoring.

The stress voltage is sqrt(P * Z), with Z the driver's nominal impedance, or
Re when no nominal impedance is recorded (the basis is reported). The model is
linear, so excursion scales with voltage; it is not a thermal or large-signal
prediction.
"""
import math

import numpy as np
import pandas as pd

from horn_core.acoustics import C0


def stress_voltage(driver, test_power_w):
    """RMS voltage for test_power_w into the nominal impedance (or Re)."""
    if driver.nominal_impedance_ohm:
        return math.sqrt(test_power_w * driver.nominal_impedance_ohm), "nominal_impedance"
    return math.sqrt(test_power_w * driver.re_ohm), "re_ohm"


def _band_grid(frequencies, low, high):
    f = np.asarray(frequencies, dtype=float)
    return np.log(np.r_[low, f[(f > low) & (f < high)], high])


def stress_assessment(frequencies, operating_point, driver, target):
    """Excursion at the test power over the full sweep and in the target band.

    Returns (metrics, rejection_reasons). Exceeding Xmax inside the target band
    rejects the design; exceeding it only below the band is reported as the
    frequency a high-pass crossover has to sit above.
    """
    f = np.asarray(frequencies, dtype=float)
    voltage, basis = stress_voltage(driver, target.test_power_w)
    scale = voltage / operating_point["drive_voltage_rms"]
    excursion = np.asarray(operating_point["displacement_peak_m"]) * scale
    grid = _band_grid(f, target.f_low_hz, target.f_high_hz)
    band_max = float(np.max(np.interp(grid, np.log(f), excursion)))
    peak_index = int(np.argmax(excursion))
    metrics = {
        "test_power_w": target.test_power_w,
        "test_voltage_rms": voltage,
        "test_voltage_basis": basis,
        "test_power_band_max_excursion_mm": band_max * 1e3,
        "test_power_sweep_max_excursion_mm": float(excursion[peak_index]) * 1e3,
        "test_power_sweep_max_excursion_hz": float(f[peak_index]),
        "test_power_xmax_margin_mm": None,
        "test_power_sweep_xmax_margin_mm": None,
        "test_power_xmax_exceeded_below_hz": None,
        "test_power_scope": "Linear model over the full simulated sweep; no crossover high-pass applied",
    }
    reasons = []
    if driver.xmax_m is not None:
        metrics["test_power_xmax_margin_mm"] = (driver.xmax_m - band_max) * 1e3
        metrics["test_power_sweep_xmax_margin_mm"] = (driver.xmax_m - float(excursion[peak_index])) * 1e3
        over = excursion > driver.xmax_m
        if over.any():
            metrics["test_power_xmax_exceeded_below_hz"] = float(f[np.flatnonzero(over)[-1]])
        if band_max > driver.xmax_m:
            reasons.append("test_power_excursion_exceeded")
    return metrics, reasons


def time_of_flight_s(distance_m, horn_length_m=None):
    """Straight-line acoustic delay through the horn axis and on to the observer."""
    return ((horn_length_m or 0.0) + distance_m) / C0


def system_group_delay(frequencies, pressure, time_of_flight):
    """Group delay -d(phase)/d(omega) in seconds, with time_of_flight removed.

    pressure uses the exp(+i w t) convention, where a pure delay t has phase
    -w t. The delay is removed before unwrapping so the remaining phase changes
    slowly enough between samples to unwrap reliably.
    """
    omega = 2 * np.pi * np.asarray(frequencies, dtype=float)
    phase = np.unwrap(np.angle(np.asarray(pressure) * np.exp(1j * omega * time_of_flight)))
    return -np.gradient(phase, omega)


def electrical_table(frame, driver, target, legacy_radius=None, horn_length_m=None):
    """Per-frequency SPL, impedance, excursion, power and group delay for one candidate."""
    from horn_analysis.evaluation import coupled_observer
    from horn_analysis.transfer_function import compute_driver_operating_point

    observed = coupled_observer(frame, driver, target, legacy_radius)
    f = frame["frequency"].to_numpy(dtype=float)
    point = compute_driver_operating_point(driver, f, frame["z_real"].to_numpy(), frame["z_imag"].to_numpy(),
                                           observed["area"], target.voltage_rms)
    ze = np.asarray(point["electrical_impedance"])
    tof = time_of_flight_s(target.observation_distance_m, horn_length_m)
    if observed["pressure"] is None:
        delay = np.full(f.shape, np.nan)  # legacy mouth-plane frames carry no observer phase
    else:
        delay = system_group_delay(f, observed["pressure"], tof)
    table = {
        "frequency": f,
        "spl": observed["levels"],
        "ze_real_ohm": ze.real,
        "ze_imag_ohm": ze.imag,
        "ze_mag_ohm": np.abs(ze),
        "excursion_peak_mm": np.asarray(point["displacement_peak_m"]) * 1e3,
        "input_power_w": np.asarray(point["input_power_w"]),
        "acoustic_power_w": np.asarray(point["horn_power_w"]),
        "group_delay_ms": delay * 1e3,
        "time_of_flight_removed_ms": np.full(f.shape, tof * 1e3),
        "drive_voltage_rms": np.full(f.shape, target.voltage_rms),
    }
    if target.test_power_w is not None:
        voltage, _ = stress_voltage(driver, target.test_power_w)
        scale = voltage / target.voltage_rms
        excursion = table["excursion_peak_mm"] * scale
        table.update({
            "test_voltage_rms": np.full(f.shape, voltage),
            "test_excursion_peak_mm": excursion,
            "test_xmax_margin_mm": (driver.xmax_m * 1e3 - excursion) if driver.xmax_m is not None else np.full(f.shape, np.nan),
            "test_input_power_w": table["input_power_w"] * scale**2,
        })
    return pd.DataFrame(table)


def plot_electrical(tables, output_png):
    """|Ze|, excursion (with Xmax) and group delay for the ranked candidates."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    for label, table, xmax_m in tables:
        line, = axes[0].semilogx(table.frequency, table.ze_mag_ohm, label=label)
        column = "test_excursion_peak_mm" if "test_excursion_peak_mm" in table else "excursion_peak_mm"
        axes[1].semilogx(table.frequency, table[column], color=line.get_color())
        if xmax_m is not None:
            axes[1].axhline(xmax_m * 1e3, color=line.get_color(), linestyle=":", linewidth=1)
        axes[2].semilogx(table.frequency, table.group_delay_ms, color=line.get_color())
    stress = any("test_excursion_peak_mm" in t for _, t, _ in tables)
    axes[0].set_ylabel("|Ze| (ohm)")
    axes[1].set_ylabel("Peak excursion (mm)" + (" at test power" if stress else " at drive voltage") + "; dotted: Xmax")
    axes[2].set_ylabel("Group delay (ms), time of flight removed")
    axes[2].set_xlabel("Frequency (Hz)")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_png, dpi=110)
    plt.close(fig)
