"""Shared labels and tables for the auto report, its summary and its CSV outputs.

Every plot legend in the auto report is built from ``candidate_label`` or
``horn_legend``, so it names a rank and an ID that appear in the ranking table.
"""
from collections import OrderedDict

import numpy as np
import pandas as pd


def candidate_id(row):
    """Driver-horn pair ID shown in the ranking table: '<horn_label> + <driver_id>'."""
    return f"{row['horn_label']} + {row['driver_id']}"


def candidate_label(rank, row):
    """Legend label for a driver-horn pair: '#<rank> <horn_label> + <driver_id>'."""
    return f"#{rank} {candidate_id(row)}"


def horn_legend(rank, row):
    """Legend label for a horn-only curve, using the best rank that uses this horn."""
    return f"#{rank} {row['horn_label']}"


def ranked_horns(rows):
    """(rank, row) for the first ranked row of each distinct horn, in rank order."""
    seen = OrderedDict()
    for rank, row in enumerate(rows, 1):
        seen.setdefault(row["horn_label"], (rank, row))
    return list(seen.values())


def has_parametric_geometry(row):
    return all(row.get(k) for k in ("profile", "throat_radius", "mouth_radius", "length")) \
        and row["mouth_radius"] > row["throat_radius"]


def profile_table(row, points=101):
    """Wall radius against axial position, both in mm, from horn_core.profiles.

    z = 0 is the throat. This is the acoustic air passage, not a wall with
    thickness.
    """
    from horn_core.profiles import radius_profile
    length = float(row["length"])
    z = np.linspace(0.0, length, points)
    r = radius_profile(z, float(row["throat_radius"]), float(row["mouth_radius"]), length, row["profile"])
    return pd.DataFrame({"z_mm": z * 1e3, "r_mm": r * 1e3})


REASON_DESCRIPTIONS = {
    "target_band_ripple_exceeded": "Target-band ripple above the limit",
    "minimum_output_not_met": "Target-band level below the minimum output",
    "maximum_length_exceeded": "Horn longer than the length limit",
    "maximum_mouth_radius_exceeded": "Mouth larger than the mouth limit",
    "compression_ratio_exceeded": "Sd / throat area above the compression limit",
    "driver_larger_than_mouth": "Diaphragm larger than the mouth",
    "outside_driver_usable_band": "Target band outside the driver's usable band",
    "driver_excursion_exceeded": "Excursion above Xmax at the drive voltage",
    "driver_nominal_power_exceeded": "Input power above the driver rating at the drive voltage",
    "test_power_excursion_exceeded": "In-band excursion above Xmax at the test power",
    "test_power_nominal_power_exceeded": "In-band input power above the driver rating at the test power",
    "radiation_model_out_of_domain": "Outside the radiation model's valid range",
    "profile_geometry_infeasible": "Profile cannot reach the mouth in this length",
    "loft_accuracy_unmet": "STEP loft could not follow the profile closely enough",
    "cad_geometry_rejected": "STEP geometry rejected",
    "band_boundary_mismatch": ("Solver bands disagree where they join; mesh not resolved. "
                               "Not a property of the design: rerun with a smaller --mesh_size"),
}


def _excess(reason, row, target, driver):
    """How far a row is past the limit behind ``reason``, as (value, unit), or (None, unit)."""
    interface = row.get("driver_interface") or {}
    def diff(a, b, scale=1.0):
        return None if a is None or b is None else (a - b) * scale
    if reason == "target_band_ripple_exceeded":
        return diff(row.get("passband_ripple_db"), target.max_ripple_db), "dB"
    if reason == "minimum_output_not_met":
        return diff(target.min_output_db, row.get("target_min_level_db")), "dB"
    if reason == "maximum_length_exceeded":
        return diff(row.get("length"), target.max_length_m, 1e3), "mm"
    if reason == "maximum_mouth_radius_exceeded":
        return diff(row.get("mouth_radius"), target.max_mouth_radius_m, 1e3), "mm"
    if reason == "compression_ratio_exceeded":
        ratio = interface.get("compression_ratio")
        if ratio is None and driver is not None and row.get("inlet_area_m2"):
            ratio = driver.sd_m2 / row["inlet_area_m2"]
        return diff(ratio, target.max_compression_ratio), "Sd/St"
    if reason == "driver_excursion_exceeded":
        return diff(row.get("max_displacement_peak_m"), getattr(driver, "xmax_m", None), 1e3), "mm"
    if reason == "driver_nominal_power_exceeded":
        return diff(row.get("max_input_power_w"), getattr(driver, "power_w", None)), "W"
    if reason == "test_power_excursion_exceeded":
        margin = row.get("test_power_xmax_margin_mm")
        return (None if margin is None else -margin), "mm"
    if reason == "test_power_nominal_power_exceeded":
        margin = row.get("test_power_rating_margin_w")
        return (None if margin is None else -margin), "W"
    if reason == "band_boundary_mismatch":
        return (row.get("band_boundary_mismatch") or {}).get("excess_ratio"), "times allowed"
    return None, ""


def rejection_summary(rejected, target, drivers=None):
    """Count rejection reasons and how far past each limit the rejected designs were.

    Returns one dict per reason, most frequent first, with ``closest_excess``
    (the nearest miss) and ``largest_excess`` in ``unit`` where a margin is
    defined, else None.
    """
    drivers = drivers or {}
    table = {}
    for row in rejected:
        driver = drivers.get(row.get("driver_id"))
        for reason in row.get("rejection_reasons") or []:
            entry = table.setdefault(reason, {"reason": reason, "description": REASON_DESCRIPTIONS.get(reason, reason.replace("_", " ")),
                                              "count": 0, "excess": [], "unit": ""})
            entry["count"] += 1
            value, unit = _excess(reason, row, target, driver)
            entry["unit"] = unit
            if value is not None and np.isfinite(value):
                entry["excess"].append(float(value))
    out = []
    for entry in sorted(table.values(), key=lambda e: (-e["count"], e["reason"])):
        excess = entry.pop("excess")
        entry["closest_excess"] = min(excess) if excess else None
        entry["largest_excess"] = max(excess) if excess else None
        out.append(entry)
    return out


def mesh_rejection_note(rejected):
    """Plain warning when the solver mesh, not the brief, rejected candidates; None otherwise.

    A band boundary mismatch means the mesh was not resolved where two solver
    bands join. Without this note a run whose every candidate failed that way
    reads as "no design fits the brief".
    """
    horns, frequencies = set(), set()
    for row in rejected:
        if "band_boundary_mismatch" in (row.get("rejection_reasons") or []):
            horns.add(row.get("horn_label"))
            frequency = (row.get("band_boundary_mismatch") or {}).get("frequency_hz")
            if frequency is not None:
                frequencies.add(round(float(frequency), 1))
    if not horns:
        return None
    where = ""
    if frequencies:
        where = " at " + ", ".join(f"{f:.0f}" for f in sorted(frequencies)) + " Hz"
    return (f"{len(horns)} simulated horn{'s were' if len(horns) != 1 else ' was'} rejected because the solver mesh "
            f"was not resolved where two frequency bands join{where}. This says nothing about whether a design "
            "fits the brief. Rerun with a smaller --mesh_size to evaluate them.")


def chamber_assumption(row):
    """One-line description of the front chamber used for a ranked row."""
    interface = row.get("driver_interface") or {}
    source = interface.get("front_chamber_source")
    if source is None:
        return "not recorded"
    volume = interface.get("front_chamber_volume_cm3") or 0.0
    text = {"specified": f"supplied {volume:.1f} cm3",
            "estimated": f"estimated {volume:.1f} cm3 (cone air, flat adapter plate)",
            "unknown": "none modelled (ideal transformer)"}.get(source, source)
    if row.get("estimated_chamber_reasons"):
        text += "; over limit only because of the estimate: " + ", ".join(row["estimated_chamber_reasons"])
    return text
