#!/usr/bin/env python3
"""Hornresp input sheets for the cross-software horns, and our matching output.

Hornresp is a Windows program and is not run here. This script writes, for
each horn, the values to type into the Hornresp input window (areas in cm2,
lengths in cm, flare cut-offs in Hz) and our own TMM throat impedance for the
same geometry, so an exported Hornresp result can be compared with
compare_hornresp.py.

Hornresp flares: Con (conical), Exp (exponential), Hyp with T (Salmon hypex,
the same family as horn_core's hypex) and Tra (tractrix). OS and the scaled
tractrix have no Hornresp flare, so each is approximated by four Con or Exp
segments; the sheet states the largest radius error, and our reference output
is computed for that same segmented geometry.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from horn_core.duct import AirProperties
from horn_core.profiles import (CD_THROAT_FRACTION, PROFILE_NAMES, cd_flare_rate, get_radius_func,
                                hypex_flare_rate, natural_tractrix_length)
from horn_core.webster import compute_horn_transfer_tmm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "validation/hornresp"
THROAT, MOUTH, LENGTH = 0.025, 0.150, 0.300
# Hornresp's usual air at 20 C; check Tools > Options in your copy and pass
# the same values to compare_hornresp.py if they differ.
HORNRESP_AIR = AirProperties(c=344.0, rho=1.205)
FREQS = np.geomspace(20.0, 20000.0, 533)
SEGMENTS = 4

EXAMPLE = dict(id="example_800_1600", profile="hyperbolic", throat_radius=0.032498,
               mouth_radius=0.08870898640584517, length=0.11790625)
DRIVER_JSON = ROOT / "data/drivers-curated/18Sound/18sound-6nmb420.json"


def area_cm2(r):
    return math.pi * r * r * 1e4


def exp_cutoff(r0, r1, length, c):
    """Hornresp Exp segment: S = S0 exp(m x), cut-off m c / (4 pi)."""
    m = 2 * math.log(r1 / r0) / length
    return m * c / (4 * math.pi)


def segment_radius(kind, r0, r1, length):
    if kind == "Con":
        return lambda x: r0 + (r1 - r0) * x / length
    return lambda x: r0 * (r1 / r0) ** (x / length)


def fit_segments(radius, length, count=SEGMENTS, grid=61):
    """Choose breakpoints and a Con or Exp flare per segment.

    Breakpoints lie on a uniform grid; a small dynamic programme minimises the
    largest relative radius error over the whole horn.
    """
    knots = np.linspace(0, length, grid)
    r_knots = np.array([radius(z) for z in knots])
    z_fine = np.linspace(0, length, 3001)
    r_fine = np.array([radius(z) for z in z_fine])
    cost = np.full((grid, grid), np.inf)
    kind = {}
    for i in range(grid):
        for j in range(i + 1, grid):
            a, b = knots[i], knots[j]
            mask = (z_fine >= a) & (z_fine <= b)
            x = z_fine[mask] - a
            con = r_knots[i] + (r_knots[j] - r_knots[i]) * x / (b - a)
            exp = r_knots[i] * (r_knots[j] / r_knots[i]) ** (x / (b - a))
            e_con = np.max(np.abs(con - r_fine[mask]) / r_fine[mask])
            e_exp = np.max(np.abs(exp - r_fine[mask]) / r_fine[mask])
            cost[i, j], kind[i, j] = (e_con, "Con") if e_con <= e_exp else (e_exp, "Exp")
    best = np.full((count + 1, grid), np.inf)
    back = np.zeros((count + 1, grid), int)
    best[0, 0] = 0.0
    for n in range(1, count + 1):
        for j in range(1, grid):
            options = np.maximum(best[n - 1, :j], cost[:j, j])
            back[n, j] = int(np.argmin(options))
            best[n, j] = options[back[n, j]]
    cuts, j = [grid - 1], grid - 1
    for n in range(count, 0, -1):
        j = back[n, j]
        cuts.append(j)
    cuts = cuts[::-1]
    segs = [(knots[i], knots[j], r_knots[i], r_knots[j], kind[i, j]) for i, j in zip(cuts[:-1], cuts[1:])]
    return segs, float(best[count, grid - 1])


def horn_sheet(name, profile, r_t, r_m, length, c):
    """Hornresp segment list and the matching radius function."""
    radius = get_radius_func(profile, r_t, r_m, length)
    note, max_err = "exact", 0.0
    if profile == "conical":
        segs = [dict(flare="Con", s_in=area_cm2(r_t), s_out=area_cm2(r_m), length_cm=length * 100)]
    elif profile == "exponential":
        segs = [dict(flare="Exp", s_in=area_cm2(r_t), s_out=area_cm2(r_m), length_cm=length * 100,
                     cutoff_hz=exp_cutoff(r_t, r_m, length, c))]
    elif profile == "hyperbolic":
        u = hypex_flare_rate(r_t, r_m, length, 0.0)
        segs = [dict(flare="Hyp", s_in=area_cm2(r_t), s_out=area_cm2(r_m), length_cm=length * 100,
                     cutoff_hz=c * u / (2 * math.pi * length), T=0.0)]
        note = "exact: Hornresp Hyp with T = 0 is S = S1 cosh^2(x / x0), F12 = c / (2 pi x0)"
    elif profile == "tractrix":
        segs = [dict(flare="Tra", s_in=area_cm2(r_t), s_out=area_cm2(r_m), length_cm=length * 100,
                     cutoff_hz=c / (2 * math.pi * r_m))]
        note = "exact if Hornresp's tractrix length for S1, S2, F12 equals the length given here"
    elif profile == "cd":
        m = cd_flare_rate(r_t, r_m, length)
        zj = CD_THROAT_FRACTION * length
        rj = r_t * math.exp(m * zj)
        segs = [dict(flare="Exp", s_in=area_cm2(r_t), s_out=area_cm2(rj), length_cm=zj * 100,
                     cutoff_hz=exp_cutoff(r_t, rj, zj, c)),
                dict(flare="Con", s_in=area_cm2(rj), s_out=area_cm2(r_m), length_cm=(length - zj) * 100)]
    else:
        fitted, max_err = fit_segments(radius, length)
        segs = []
        for a, b, ra, rb, kind in fitted:
            seg = dict(flare=kind, s_in=area_cm2(ra), s_out=area_cm2(rb), length_cm=(b - a) * 100)
            if kind == "Exp":
                seg["cutoff_hz"] = exp_cutoff(ra, rb, b - a, c)
            segs.append(seg)
        note = f"approximated by {len(segs)} segments, max radius error {max_err * 100:.2f} %"
    return {"id": name, "profile": profile, "throat_radius_m": r_t, "mouth_radius_m": r_m, "length_m": length,
            "geometry_note": note, "max_radius_error": max_err, "segments": segs}


def segmented_radius(sheet):
    """Radius function of the geometry Hornresp is given (segments as entered)."""
    if sheet["profile"] in ("hyperbolic", "tractrix"):
        return get_radius_func(sheet["profile"], sheet["throat_radius_m"], sheet["mouth_radius_m"], sheet["length_m"])
    pieces, z0 = [], 0.0
    for seg in sheet["segments"]:
        L = seg["length_cm"] / 100
        r0, r1 = math.sqrt(seg["s_in"] / 1e4 / math.pi), math.sqrt(seg["s_out"] / 1e4 / math.pi)
        pieces.append((z0, z0 + L, segment_radius(seg["flare"], r0, r1, L)))
        z0 += L

    def radius(z):
        for a, b, f in pieces:
            if z <= b:
                return f(z - a)
        return pieces[-1][2](pieces[-1][1] - pieces[-1][0])
    return radius


def our_output(sheet, air=HORNRESP_AIR, freqs=FREQS):
    out = compute_horn_transfer_tmm(freqs, segmented_radius(sheet), sheet["length_m"], sheet["throat_radius_m"],
                                    sheet["mouth_radius_m"], n_segments=2000, radiation_model="flanged_piston",
                                    air=air)
    z = (out["z_real"] + 1j * out["z_imag"]) / (air.rho * air.c)
    return freqs, z


def halfspace_spl(sheet, air=HORNRESP_AIR, freqs=None, voltage=2.83, distance=1.0):
    """On-axis SPL built the way a 1D horn program builds it.

    Our TMM horn (baffled piston mouth load) driven through an ideal area
    transformer by the sheet's driver; the mouth volume velocity then radiates
    as a point source into half space, p = rho w |U| / (2 pi r), with no
    piston directivity. This is the like-for-like reference for a Hornresp
    SPL chart at Ang = 2 pi. Returns (freqs, spl_db).
    """
    freqs = np.geomspace(500.0, 3000.0, 241) if freqs is None else np.asarray(freqs, float)
    d = sheet["driver"]
    sd, s_t = d["Sd_cm2"] / 1e4, math.pi * sheet["throat_radius_m"] ** 2
    out = compute_horn_transfer_tmm(freqs, segmented_radius(sheet), sheet["length_m"], sheet["throat_radius_m"],
                                    sheet["mouth_radius_m"], n_segments=2000, radiation_model="flanged_piston",
                                    air=air)
    z_horn = out["z_real"] + 1j * out["z_imag"]          # specific, Pa s/m
    w = 2 * np.pi * freqs
    z_e = d["Re_ohm"] + 1j * w * d["Le_mH"] / 1e3
    z_mech = (d["Rms_Ns_per_m"] + 1j * w * d["Mmd_g"] / 1e3 + 1 / (1j * w * d["Cms_m_per_N"])
              + z_horn * sd ** 2 / s_t)
    velocity = d["Bl_Tm"] * voltage / (z_e * z_mech + d["Bl_Tm"] ** 2)
    p_throat = z_horn * velocity * sd / s_t
    u_mouth = p_throat * out["mouth_volume_velocity_transfer"]
    p_far = air.rho * w * np.abs(u_mouth) / (2 * np.pi * distance)
    return freqs, 20 * np.log10(p_far / 20e-6)


def driver_values():
    p = json.loads(DRIVER_JSON.read_text())["parameters"]
    w = 2 * math.pi * p["fs_hz"]
    cms = 1 / (w * w * p["mms_kg"])
    rms = w * p["mms_kg"] / p["qms"]
    return {"Sd_cm2": p["sd_m2"] * 1e4, "Bl_Tm": p["bl_tm"], "Cms_m_per_N": cms, "Rms_Ns_per_m": rms,
            "Mmd_g": p["mms_kg"] * 1e3, "Le_mH": p["le_h"] * 1e3, "Re_ohm": p["re_ohm"],
            "note": ("Cms and Rms derived from fs, Mms and Qms as in this project. Enter Mms as Mmd: our model uses "
                     "the published Mms as the whole moving mass with no rear chamber (Vrc = 0, open back) and an "
                     "ideal area transformer (Vtc = 0).")}


def sheet_text(sheet, c, rho, driver=None):
    lines = [f"Hornresp input sheet: {sheet['id']}  ({sheet['profile']})",
             f"Geometry: throat radius {sheet['throat_radius_m']*1e3:.3f} mm, mouth radius "
             f"{sheet['mouth_radius_m']*1e3:.3f} mm, length {sheet['length_m']*1e3:.3f} mm",
             f"Geometry match: {sheet['geometry_note']}",
             "", "Radiation and source",
             "  Ang = 2.0 x Pi   (half space: horn mouth in an infinite baffle)",
             "  Eg  = 2.83 V     Rg = 0", "  Cir = 0",
             f"  Air (Tools > Options): c = {c:g} m/s, rho = {rho:g} kg/m3", "", "Horn"]
    s_keys = ["S1", "S2", "S3", "S4", "S5"]
    for i, seg in enumerate(sheet["segments"]):
        j = i + 1
        lines.append(f"  {s_keys[i]} = {seg['s_in']:.2f} cm2")
        entry = f"  L{j}{j+1} = {seg['length_cm']:.2f} cm, flare {seg['flare']}"
        if "cutoff_hz" in seg:
            entry += f", F{j}{j+1} = {seg['cutoff_hz']:.2f} Hz"
        if "T" in seg:
            entry += f", T = {seg['T']:.2f}"
        lines.append(entry)
    lines.append(f"  {s_keys[len(sheet['segments'])]} = {sheet['segments'][-1]['s_out']:.2f} cm2")
    if driver:
        lines += ["", "Driver (18Sound 6NMB420, from data/drivers-curated)",
                  f"  Sd = {driver['Sd_cm2']:.2f} cm2   Bl = {driver['Bl_Tm']:.2f} Tm",
                  f"  Cms = {driver['Cms_m_per_N']:.4E} m/N   Rms = {driver['Rms_Ns_per_m']:.3f} Ns/m",
                  f"  Mmd = {driver['Mmd_g']:.2f} g   Le = {driver['Le_mH']:.3f} mH   Re = {driver['Re_ohm']:.2f} ohm",
                  "  Vrc = 0, Lrc = 0, Vtc = 0, Atc = S1", f"  Note: {driver['note']}"]
    else:
        lines += ["", "Driver: none needed. Read the throat acoustic impedance (Ra, Xa, normalised to rho c / S1)."]
    lines += ["", "Export: Tools > Export (chart data, text) over the default 10 Hz to 20 kHz range,",
              "then run validation/hornresp/compare_hornresp.py with the exported file."]
    return "\n".join(lines) + "\n"


def import_record_text(sheet, driver=None):
    """Best-effort Hornresp import record. UNVERIFIED.

    Written from the section and key layout of Hornresp text records as
    published in user forums (ID, Comment, "|...PARAMETER VALUES:" sections,
    one key=value per line, each horn segment given as its flare keyword
    holding the segment length in cm, followed by F<ij> for the flare
    cut-off). It has not been checked against Hornresp itself. If the import
    is rejected or a value lands in the wrong field, enter the sheet in
    records/<id>.txt by hand instead.
    """
    segs = sheet["segments"]
    lines = ["ID=0.00", f"Comment=horn-simulation cross-software check: {sheet['id']} ({sheet['profile']})",
             "", "|RADIATION, SOURCE AND MOUTH PARAMETER VALUES:", "", "Ang=2.0 x Pi", "Eg=2.83", "Rg=0.00",
             "Cir=0.00", "", "|HORN PARAMETER VALUES:", ""]
    for i, seg in enumerate(segs):
        j = i + 1
        lines.append(f"S{j}={seg['s_in']:.2f}")
        lines.append(f"S{j+1}={seg['s_out']:.2f}")
        lines.append(f"{seg['flare']}={seg['length_cm']:.2f}")
        lines.append(f"F{j}{j+1}={seg.get('cutoff_hz', 0.0):.2f}")
        if "T" in seg:
            lines.append(f"T={seg['T']:.2f}")
    if driver:
        lines += ["", "|TRADITIONAL DRIVER PARAMETER VALUES:", "",
                  f"Sd={driver['Sd_cm2']:.2f}", f"Bl={driver['Bl_Tm']:.2f}", f"Cms={driver['Cms_m_per_N']:.2E}",
                  f"Rms={driver['Rms_Ns_per_m']:.2f}", f"Mmd={driver['Mmd_g']:.2f}", f"Le={driver['Le_mH']:.2f}",
                  f"Re={driver['Re_ohm']:.2f}", "Nd=1", "", "|CHAMBER PARAMETER VALUES:", "",
                  "Vrc=0.00", "Lrc=0.00", "Vtc=0.00", f"Atc={segs[0]['s_in']:.2f}"]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    c, rho = HORNRESP_AIR.c, HORNRESP_AIR.rho
    sheets = []
    for profile in PROFILE_NAMES:
        L = natural_tractrix_length(THROAT, MOUTH) if profile == "tractrix" else LENGTH
        sheets.append(horn_sheet(profile, profile, THROAT, MOUTH, L, c))
    ex = horn_sheet(EXAMPLE["id"], EXAMPLE["profile"], EXAMPLE["throat_radius"], EXAMPLE["mouth_radius"],
                    EXAMPLE["length"], c)
    ex["driver"] = driver_values()
    ex["driver_response_reference"] = "examples/800-1600-modal/response.csv (modal_baffled FEM, ideal transformer)"
    sheets.append(ex)
    (args.out / "records").mkdir(parents=True, exist_ok=True)
    (args.out / "reference").mkdir(parents=True, exist_ok=True)
    (args.out / "records" / "import_unverified").mkdir(parents=True, exist_ok=True)
    for sheet in sheets:
        (args.out / "records" / f"{sheet['id']}.txt").write_text(sheet_text(sheet, c, rho, sheet.get("driver")))
        (args.out / "records" / "import_unverified" / f"{sheet['id']}.txt").write_text(
            import_record_text(sheet, sheet.get("driver")))
        f, z = our_output(sheet)
        rows = ["frequency_hz,ra_norm,xa_norm"] + [f"{a:.6f},{b.real:.8g},{b.imag:.8g}" for a, b in zip(f, z)]
        (args.out / "reference" / f"{sheet['id']}_tmm.csv").write_text("\n".join(rows) + "\n")
    f, spl = halfspace_spl(ex)
    rows = ["frequency_hz,spl_db"] + [f"{a:.6f},{b:.6f}" for a, b in zip(f, spl)]
    (args.out / "reference" / f"{ex['id']}_spl_halfspace.csv").write_text("\n".join(rows) + "\n")
    (args.out / "records" / "horns.json").write_text(json.dumps(
        {"air": {"c": c, "rho": rho}, "horns": sheets}, indent=1) + "\n")
    for s in sheets:
        print(f"{s['id']:18s} {s['geometry_note']}")


if __name__ == "__main__":
    main()
