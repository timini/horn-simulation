#!/usr/bin/env python3
"""Exterior radiation checks with bempp-cl (MIT) boundary elements.

Runs in the separate environment from validation/cross_software/requirements.txt
(numba backend, no OpenCL) with packages/horn-core/src on PYTHONPATH.

Cases
  sphere   pulsating sphere against its closed form (checks the BEM set-up).
  piston   rigid piston in a large finite circular baffle (a thick disc)
           against the analytic infinite-baffle piston and the project's
           modal far field (horn_core.modal_radiation) with uniform velocity.
  horn     the archived worked-example hyperbolic horn cut into the front of
           the same kind of baffle, driven by a uniform throat velocity.
           Throat impedance is compared with the archived modal_baffled FEM
           solve; the far field with the project's experimental directivity
           computed from the archived FEM mouth velocities.

bempp uses exp(-i w t); results are conjugated into the project's exp(+i w t)
before comparison. A finite baffle is not an infinite one: the expected
differences are written next to each tolerance.
"""
import argparse
import csv
import io
import json
import tarfile
import time
from pathlib import Path

import numpy as np
from scipy.special import j1

from horn_core.modal_radiation import modal_baffled_far_field
from horn_core.profiles import get_radius_func
from horn_core.webster import compute_horn_transfer_tmm, piston_radiation_impedance

ROOT = Path(__file__).resolve().parents[2]
RHO, C = 1.225, 343.0
ANGLES = np.arange(0.0, 91.0, 5.0)
FAR_R = 50.0

# Archived worked example (docs/MODAL_APERTURE.md): refine_hyperbolic_0001.
HORN = dict(profile="hyperbolic", throat_radius=0.032498, mouth_radius=0.08870898640584517, length=0.11790625)
HORN_ARCHIVE = ROOT / "data/validation/worked_example_modal_800_1600.tar.gz"
HORN_CSV = "worked-example-modal/outputs/auto/refinement/refine_hyperbolic_0001_results.csv"
HORN_FREQS = [800.0, 1000.0, 1200.0, 1400.0, 1600.0]

PISTON_A = 0.1
PISTON_KA = [0.5, 1.0, 2.0, 4.0]
BAFFLE_RATIO = 5.0      # baffle radius / radiator radius
BAFFLE_THICKNESS = 0.05

# Frozen before the first BEM result was inspected. Do not loosen.
TOLERANCES = {
    "sphere": {"level_db": 0.1, "phase_deg": 1.0},
    # Finite baffle of radius 5a: edge diffraction is expected, so impedance
    # and directivity are checked only for ka >= 2 and theta <= 60 deg, and
    # directivity only where the analytic level is above -20 dB (off nulls).
    "piston_impedance": {"level_db": 0.5, "phase_deg": 5.0, "ka_min": 2.0},
    "piston_directivity": {"level_db": 1.0, "theta_max_deg": 60.0, "floor_db": -20.0, "ka_min": 2.0},
    "modal_uniform_vs_analytic": {"level_db": 1e-6},
    "horn_throat_impedance": {"level_db": 0.5, "phase_deg": 5.0},
    "horn_directivity": {"level_db": 1.0, "theta_max_deg": 60.0, "floor_db": -20.0},
    "horn_on_axis_level": {"level_db": 1.5},
    "mesh_refinement": {"level_db": 0.1, "phase_deg": 1.0},
}


def bempp():
    import bempp_cl.api as api
    api.DEFAULT_DEVICE_INTERFACE = "numba"
    api.enable_console_logging("warning") if hasattr(api, "enable_console_logging") else None
    return api


def oriented_surface(volume_tag_pairs, driven_surfaces):
    """Mesh the current gmsh model and return outward-oriented boundary triangles.

    Orientation comes from the tetrahedron behind each boundary face, so the
    normal always points out of the body into the air.
    """
    import gmsh
    gmsh.model.mesh.generate(3)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    index = {int(t): i for i, t in enumerate(node_tags)}
    xyz = coords.reshape(-1, 3)
    _, _, tet_nodes = gmsh.model.mesh.getElements(3)
    tets = np.vectorize(index.get)(np.concatenate(tet_nodes).astype(int)).reshape(-1, 4)
    opposite = {}
    for t in tets:
        for i in range(4):
            face = tuple(sorted(np.delete(t, i)))
            opposite[face] = t[i]
    triangles, domains = [], []
    for dim, tag in gmsh.model.getBoundary(volume_tag_pairs, oriented=False):
        _, _, tri_nodes = gmsh.model.mesh.getElements(2, abs(tag))
        tri = np.vectorize(index.get)(np.concatenate(tri_nodes).astype(int)).reshape(-1, 3)
        for t in tri:
            inner = opposite.get(tuple(sorted(t)))
            if inner is None:
                raise RuntimeError("Boundary face without a tetrahedron")
            n = np.cross(xyz[t[1]] - xyz[t[0]], xyz[t[2]] - xyz[t[0]])
            if np.dot(n, xyz[inner] - xyz[t[0]]) > 0:
                t = t[[0, 2, 1]]
            triangles.append(t)
            domains.append(1 if abs(tag) in driven_surfaces else 2)
    triangles = np.array(triangles)
    used = np.unique(triangles)
    remap = -np.ones(len(xyz), int)
    remap[used] = np.arange(len(used))
    return xyz[used].T.copy(), remap[triangles].T.copy(), np.array(domains, dtype=np.uint32)


def baffle_with(radiator, h_near, h_far, horn=None):
    """Thick circular baffle with a flush piston, or with a horn cut into it."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("baffle")
    occ = gmsh.model.occ
    b = BAFFLE_RATIO * radiator
    if horn is None:
        body = occ.addCylinder(0, 0, -BAFFLE_THICKNESS, 0, 0, BAFFLE_THICKNESS, b)
        disk = occ.addDisk(0, 0, 0, radiator, radiator)
        out, _ = occ.fragment([(3, body)], [(2, disk)])
        occ.synchronize()
        vols = [p for p in out if p[0] == 3]
    else:
        L = horn["length"]
        depth = L + BAFFLE_THICKNESS
        body = occ.addCylinder(0, 0, -depth, 0, 0, depth, b)
        radius = get_radius_func(horn["profile"], horn["throat_radius"], horn["mouth_radius"], L)
        zs = np.linspace(0.0, L, 60)
        pts = [occ.addPoint(float(radius(z)), 0, float(z - L)) for z in zs]
        spline = occ.addSpline(pts)
        top = occ.addPoint(0, 0, 0)
        bottom = occ.addPoint(0, 0, -L)
        wire = occ.addCurveLoop([occ.addLine(bottom, pts[0]), spline, occ.addLine(pts[-1], top), occ.addLine(top, bottom)])
        face = occ.addPlaneSurface([wire])
        cavity = occ.revolve([(2, face)], 0, 0, 0, 0, 0, 1, 2 * np.pi)
        cavity_vol = [p for p in cavity if p[0] == 3]
        out, _ = occ.cut([(3, body)], cavity_vol)
        occ.synchronize()
        vols = [p for p in out if p[0] == 3]
    driven = set()
    for dim, tag in gmsh.model.getBoundary(vols, oriented=False):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, abs(tag))
        eps = 1e-5   # OCC bounding boxes carry a small tolerance
        flat = abs(zmax - zmin) < eps
        if horn is None and flat and abs(zmax) < eps and xmax < radiator * 1.001:
            driven.add(abs(tag))
        if horn is not None and flat and abs(zmax + horn["length"]) < eps and xmax < horn["throat_radius"] * 1.001:
            driven.add(abs(tag))
    if len(driven) != 1:
        raise RuntimeError(f"Expected one driven surface, found {driven}")
    field = gmsh.model.mesh.field
    dist = field.add("Distance")
    # Refine near the radiator: the piston disc, or the horn wall and throat.
    near = [abs(t) for d, t in gmsh.model.getBoundary(vols, oriented=False)
            if gmsh.model.getBoundingBox(2, abs(t))[3] < radiator * 1.001]
    field.setNumbers(dist, "SurfacesList", near)
    thr = field.add("Threshold")
    field.setNumber(thr, "InField", dist)
    field.setNumber(thr, "SizeMin", h_near)
    field.setNumber(thr, "SizeMax", h_far)
    field.setNumber(thr, "DistMin", 0.5 * radiator)
    field.setNumber(thr, "DistMax", 3 * radiator)
    field.setAsBackgroundMesh(thr)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    try:
        return oriented_surface(vols, driven)
    finally:
        gmsh.finalize()


def solve_exterior(api, vertices, elements, domains, k, velocity_fn, points):
    """Exterior Neumann problem. velocity_fn(centroids) gives the normal velocity on driven faces."""
    grid = api.Grid(vertices, elements, domains)
    p1 = api.function_space(grid, "P", 1)
    dp0 = api.function_space(grid, "DP", 0)
    centroids = grid.centroids
    v = np.zeros(grid.number_of_elements, dtype=complex)
    driven = grid.domain_indices == 1
    v[driven] = velocity_fn(centroids[driven])
    neumann = api.GridFunction(dp0, coefficients=1j * k * C * RHO * v)   # dp/dn = i w rho v_n
    ident = api.operators.boundary.sparse.identity(p1, p1, p1)
    dlp = api.operators.boundary.helmholtz.double_layer(p1, p1, p1, k)
    slp = api.operators.boundary.helmholtz.single_layer(dp0, p1, p1, k)
    from bempp_cl.api.linalg import gmres
    trace, info = gmres(dlp - 0.5 * ident, slp * neumann, tol=1e-8)
    if info != 0:
        raise RuntimeError(f"GMRES did not converge (info={info})")
    far = (api.operators.potential.helmholtz.double_layer(p1, points, k) * trace
           - api.operators.potential.helmholtz.single_layer(dp0, points, k) * neumann).ravel()
    surface = trace.evaluate_on_element_centers().ravel()
    area = grid.volumes
    mean_p = np.sum(surface[driven] * area[driven]) / np.sum(area[driven])
    return {"far": far, "mean_driven_pressure": mean_p, "dofs": p1.global_dof_count,
            "elements": int(grid.number_of_elements)}


def far_points(radius=FAR_R):
    t = np.radians(ANGLES)
    return np.vstack([radius * np.sin(t), np.zeros_like(t), radius * np.cos(t)])


def db_err(a, b):
    return float(np.max(np.abs(20 * np.log10(np.abs(a) / np.abs(b)))))


def ph_err(a, b):
    return float(np.max(np.abs(np.degrees(np.angle(a / b)))))


def directivity_db(far):
    return 20 * np.log10(np.abs(far) / np.abs(far[0]))


def directivity_check(ours, ref, tol):
    keep = (ANGLES <= tol["theta_max_deg"]) & (ref > tol["floor_db"])
    err = np.abs(ours[keep] - ref[keep])
    return float(err.max()), [float(x) for x in ANGLES[keep]]


def sphere_case(api):
    a, results = 0.1, []
    for ka in (0.5, 2.0, 5.0):
        k = ka / a
        g = api.shapes.regular_sphere(5)
        pts = far_points()
        sol = solve_exterior(api, g.vertices * a, g.elements, np.ones(g.number_of_elements, np.uint32), k,
                             lambda x: np.ones(x.shape[0]), pts)
        exact = RHO * C * a / FAR_R * (1j * ka) / (1j * ka - 1) * np.exp(1j * k * (FAR_R - a))
        lvl, ph = db_err(sol["far"], exact), ph_err(sol["far"], exact)
        tol = TOLERANCES["sphere"]
        results.append({"ka": ka, "dofs": sol["dofs"], "level_db_max": lvl, "phase_deg_max": ph,
                        "pass": lvl <= tol["level_db"] and ph <= tol["phase_deg"]})
        print(f"sphere ka={ka}: {lvl:.4f} dB {ph:.3f} deg", flush=True)
    return results


def piston_case(api, refine):
    rows = []
    for ka in PISTON_KA:
        k = ka / PISTON_A
        lam = 2 * np.pi / k
        h_near = min(lam / 8, PISTON_A / 6) / refine
        h_far = min(lam / 6, BAFFLE_RATIO * PISTON_A / 4) / refine
        verts, elems, dom = baffle_with(PISTON_A, h_near, h_far)
        t0 = time.time()
        sol = solve_exterior(api, verts, elems, dom, k, lambda x: np.ones(x.shape[0]), far_points())
        far = np.conj(sol["far"] * FAR_R * np.exp(-1j * k * FAR_R))       # R p e^{ikR}, exp(+iwt)
        z_bem = np.conj(sol["mean_driven_pressure"]) / (RHO * C)
        z_piston = piston_radiation_impedance(k, PISTON_A)
        f = np.array([k * C / (2 * np.pi)])
        modal = modal_baffled_far_field(f, np.r_[1.0, np.zeros(15)][None, :], PISTON_A, ANGLES, rho=RHO, c=C)[0]
        x = k * PISTON_A * np.sin(np.radians(ANGLES))
        with np.errstate(invalid="ignore", divide="ignore"):
            pattern = np.where(x > 1e-12, 2 * j1(x) / x, 1.0)
        analytic = 1j * RHO * C * k * PISTON_A**2 / 2 * pattern      # infinite-baffle piston, v = 1
        dir_bem, dir_analytic = directivity_db(far), directivity_db(analytic)
        tol_z, tol_d = TOLERANCES["piston_impedance"], TOLERANCES["piston_directivity"]
        dir_err, dir_angles = directivity_check(dir_bem, dir_analytic, tol_d)
        z_lvl, z_ph = db_err(z_bem, z_piston), ph_err(z_bem, z_piston)
        checked = ka >= tol_z["ka_min"]
        row = {"ka": ka, "frequency_hz": float(f[0]), "dofs": sol["dofs"], "elements": sol["elements"],
               "seconds": round(time.time() - t0, 1),
               "z_bem_over_rhoc": [z_bem.real, z_bem.imag], "z_piston_over_rhoc": [z_piston.real, z_piston.imag],
               "z_level_db": z_lvl, "z_phase_deg": z_ph,
               "on_axis_level_bem_minus_infinite_db": float(20 * np.log10(abs(far[0]) / abs(analytic[0]))),
               "directivity_bem_db": dir_bem.tolist(), "directivity_analytic_db": dir_analytic.tolist(),
               "directivity_err_db": dir_err, "directivity_angles_checked": dir_angles,
               "modal_uniform_vs_analytic_db": db_err(modal, analytic),
               "checked": checked,
               "pass": (not checked) or (z_lvl <= tol_z["level_db"] and z_ph <= tol_z["phase_deg"]
                                         and dir_err <= tol_d["level_db"]),
               "modal_pass": db_err(modal, analytic) <= TOLERANCES["modal_uniform_vs_analytic"]["level_db"]}
        rows.append(row)
        print(f"piston ka={ka} dofs={sol['dofs']}: Z {z_lvl:.3f} dB {z_ph:.2f} deg, dir {dir_err:.3f} dB, "
              f"axis {row['on_axis_level_bem_minus_infinite_db']:+.2f} dB ({row['seconds']} s)", flush=True)
    return rows


def archived_horn():
    with tarfile.open(HORN_ARCHIVE) as tar:
        text = tar.extractfile(HORN_CSV).read().decode()
    rows = list(csv.DictReader(io.StringIO(text)))
    return rows


def horn_case(api, refine):
    rows = archived_horn()
    freqs_fem = np.array([float(r["frequency"]) for r in rows])
    out = []
    radius = HORN["mouth_radius"]
    tmm = compute_horn_transfer_tmm(np.array(HORN_FREQS), get_radius_func(HORN["profile"], HORN["throat_radius"],
                                    HORN["mouth_radius"], HORN["length"]), HORN["length"], HORN["throat_radius"],
                                    HORN["mouth_radius"], n_segments=200, radiation_model="flanged_piston")
    for i, f in enumerate(HORN_FREQS):
        j = int(np.argmin(np.abs(freqs_fem - f)))
        r = rows[j]
        f_fem = float(r["frequency"])
        k = 2 * np.pi * f_fem / C
        lam = 2 * np.pi / k
        h_near = min(lam / 10, HORN["throat_radius"] / 3) / refine
        h_far = (lam / 6) / refine
        verts, elems, dom = baffle_with(radius, h_near, h_far, horn=HORN)
        t0 = time.time()
        sol = solve_exterior(api, verts, elems, dom, k, lambda x: np.ones(x.shape[0]), far_points())
        z_bem = np.conj(sol["mean_driven_pressure"])                 # specific, v = 1
        # inlet_u is a volume velocity for unit inlet pressure; convert to particle velocity.
        u_in = (float(r["inlet_u_real"]) + 1j * float(r["inlet_u_imag"])) / float(r["inlet_area_m2"])
        z_fem = (float(r["z_real"]) + 1j * float(r["z_imag"]))
        z_tmm = tmm["z_real"][i] + 1j * tmm["z_imag"][i]
        far_bem = np.conj(sol["far"] * FAR_R * np.exp(-1j * k * FAR_R))  # per unit throat velocity
        modal_v = np.array([[float(r[f"modal_v_{n}_real"]) + 1j * float(r[f"modal_v_{n}_imag"])
                             for n in range(16)]]) / u_in            # per unit throat velocity
        far_modal = modal_baffled_far_field(np.array([f_fem]), modal_v, radius, ANGLES, rho=RHO, c=C)[0]
        dir_bem, dir_modal = directivity_db(far_bem), directivity_db(far_modal)
        tol_d = TOLERANCES["horn_directivity"]
        dir_err, dir_angles = directivity_check(dir_bem, dir_modal, tol_d)
        axis_db = float(20 * np.log10(abs(far_bem[0]) / abs(far_modal[0])))
        z_lvl, z_ph = db_err(z_bem, z_fem), ph_err(z_bem, z_fem)
        tz = TOLERANCES["horn_throat_impedance"]
        row = {"frequency_hz": f_fem, "ka_mouth": k * radius, "dofs": sol["dofs"], "elements": sol["elements"],
               "seconds": round(time.time() - t0, 1),
               "z_bem_over_rhoc": [z_bem.real / (RHO * C), z_bem.imag / (RHO * C)],
               "z_fem_modal_over_rhoc": [z_fem.real / (RHO * C), z_fem.imag / (RHO * C)],
               "z_tmm_flanged_over_rhoc": [z_tmm.real / (RHO * C), z_tmm.imag / (RHO * C)],
               "z_level_db": z_lvl, "z_phase_deg": z_ph,
               "z_tmm_vs_bem_level_db": db_err(z_tmm, z_bem), "z_tmm_vs_bem_phase_deg": ph_err(z_tmm, z_bem),
               "directivity_bem_db": dir_bem.tolist(), "directivity_modal_db": dir_modal.tolist(),
               "directivity_err_db": dir_err, "directivity_angles_checked": dir_angles,
               "on_axis_level_bem_minus_modal_db": axis_db,
               "pass_impedance": z_lvl <= tz["level_db"] and z_ph <= tz["phase_deg"],
               "pass_directivity": dir_err <= tol_d["level_db"],
               "pass_on_axis": abs(axis_db) <= TOLERANCES["horn_on_axis_level"]["level_db"]}
        out.append(row)
        print(f"horn f={f_fem:.0f} dofs={sol['dofs']}: Z {z_lvl:.3f} dB {z_ph:.2f} deg (TMM {row['z_tmm_vs_bem_level_db']:.3f} dB), "
              f"dir {dir_err:.3f} dB, axis {axis_db:+.2f} dB ({row['seconds']} s)", flush=True)
    return out


def refinement(coarse, fine, key):
    out = []
    for a, b in zip(coarse, fine):
        za, zb = complex(*a[key]), complex(*b[key])
        lvl, ph = db_err(za, zb), ph_err(za, zb)
        d = float(np.max(np.abs(np.array(a["directivity_bem_db"]) - np.array(b["directivity_bem_db"]))[ANGLES <= 60]))
        tol = TOLERANCES["mesh_refinement"]
        out.append({"dofs": [a["dofs"], b["dofs"]], "z_level_db": lvl, "z_phase_deg": ph, "directivity_db_max_le60": d,
                    "pass": lvl <= tol["level_db"] and ph <= tol["phase_deg"] and d <= tol["level_db"]})
    return out


def plot(path, piston, horn):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for row in piston:
        line, = axes[0].plot(ANGLES, row["directivity_bem_db"], "-", lw=1.2, label=f"BEM ka={row['ka']}")
        axes[0].plot(ANGLES, row["directivity_analytic_db"], "--", color=line.get_color(), lw=1)
    axes[0].set_title("Piston, finite baffle (BEM, solid) vs infinite (analytic, dashed)")
    for row in horn:
        line, = axes[1].plot(ANGLES, row["directivity_bem_db"], "-", lw=1.2, label=f"BEM {row['frequency_hz']:.0f} Hz")
        axes[1].plot(ANGLES, row["directivity_modal_db"], "--", color=line.get_color(), lw=1)
    axes[1].set_title("Hyperbolic horn: BEM (solid) vs modal FEM far field (dashed)", fontsize=10)
    for ax in axes:
        ax.set_xlabel("angle from axis (deg)")
        ax.set_ylabel("level re axis (dB)")
        ax.set_ylim(-30, 3)
        ax.grid(alpha=.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=ROOT / "data/validation/cross_software_bempp.json")
    parser.add_argument("--plot", type=Path, default=ROOT / "docs/validation/cross_software/bempp_directivity.png")
    parser.add_argument("--skip-refinement", action="store_true")
    args = parser.parse_args()
    api = bempp()
    t0 = time.time()
    sphere = sphere_case(api)
    piston = piston_case(api, 1.0)
    horn = horn_case(api, 1.0)
    summary = {"tool": "bempp-cl", "tool_version": api.__version__, "tool_license": "MIT",
               "backend": "numba", "phasor_note": "bempp exp(-iwt) results conjugated to exp(+iwt)",
               "air": {"rho": RHO, "c": C}, "far_field_radius_m": FAR_R,
               "baffle": {"radius_over_radiator": BAFFLE_RATIO, "thickness_m": BAFFLE_THICKNESS},
               "horn": HORN, "horn_reference": HORN_CSV, "tolerances": TOLERANCES,
               "sphere": sphere, "piston": piston, "horn_results": horn}
    if not args.skip_refinement:
        piston_fine = piston_case(api, 1.4)
        horn_fine = horn_case(api, 1.4)
        summary["refinement"] = {"piston": refinement(piston, piston_fine, "z_bem_over_rhoc"),
                                 "horn": refinement(horn, horn_fine, "z_bem_over_rhoc")}
        summary["piston_fine"], summary["horn_fine"] = piston_fine, horn_fine
    summary["seconds"] = round(time.time() - t0, 1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    args.plot.parent.mkdir(parents=True, exist_ok=True)
    plot(args.plot, piston, horn)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
