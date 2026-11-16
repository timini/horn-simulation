"""Frozen, unfitted comparison with the Hajicek-Neely conical horn measurements.

Source: Zenodo record 19710026 (CC BY 4.0), imported by horn_analysis.reference_data.
The horn is a 3D-printed conical standard: 7 mm throat tube, 120 mm cone to
28.57 mm, 10 mm mouth cylinder, then either a rigid cap or an open mouth. The
measured quantity is probe load impedance inside the throat tube.

Every setting below is fixed before any measured value is read. Nominal-geometry
predictions and the protocol are written and hashed first. The only quantity
taken from the measurements is the probe plane for the time-windowed check,
estimated by the authors' published time-of-flight procedure (std_hrn_rfl.m,
ported, not executed); it is applied identically to model and measurement and
is never adjusted to improve agreement. Nothing else is fitted.

Prediction and code hashes in protocol.json are reproducible only in the same
environment: CSV bytes differ in the last digits between the horn-solver image
and a macOS host. The committed summary comes from the horn-solver image.

Run from the repository root. TMM needs only horn-core. --with-fem needs the
horn-solver image (dolfinx) with the repository mounted.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from horn_core.duct import AirProperties
from horn_core.webster import compute_horn_transfer_tmm

# Authors' 25 C air (std_hrn_rfl.m), with the 25 C transport properties already
# used for the Ernoult benchmark in validate_duct_physics.py.
T_C = 25.
AIR = AirProperties(c=34723e-2*(1+.00166*(T_C-26.85)), rho=1.1769*(1-.00335*(T_C-26.85)),
                    gamma=1.402, viscosity=1.83183e-5, conductivity=.02613336, heat_capacity=1004.16)
CGS_TO_SI = 1e5  # cgs acoustic ohm to Pa*s/m^3
REPO = Path(__file__).resolve().parents[1]

# Nominal geometry from spec.txt / fun_con.txt, checked against hrn_con.stl.
R_THROAT = .0035
CONE_LENGTH = .12
R_MOUTH = .0035*(1+(np.sqrt(10*.07**2/(np.pi*.035**2))-1)/10*12)  # fun_con.txt at x = 120 mm
MOUTH_CYLINDER = .01
WALL = .002
NOMINAL_PROBE_TO_FLARE = .025  # probe at the throat tube entrance: an assumption
CAP_DISTANCE_AUTHORS = .135     # std_hrn_rfl.m 'corrected distance from flare-start to cap'
SWEET_SPOT = .08                # std_hrn_rfl.m window: 8 cm past flare start
TRUNCATION_HZ = 16000.
BAND = (100., 8000.)
OPEN_BAND = (100., 5500.)       # finite-flange radiation needs ka < 1.5
SEGMENTS = 600
LOSS = 'kirchhoff'

FILES = {'22E27P01.THL': 'open', '22E27P02.THL': 'closed'}
AUTHOR_RESIDUAL_LENGTH_M = {'22E27P01.THL': .0161, '22E27P02.THL': .0132}  # fig1.png, fig2.png
LIMITS = {
    'untruncated_capped': {'median_abs_reflectance_error': .02, 'p95_abs_reflectance_error': .05},
    'time_windowed': {'max_abs_reflectance_error': .05},
    'tmm_vs_analytic_cone_windowed': {'max_abs_reflectance_error': .02},
    'impedance_magnitude_db': {'median': 2., 'p95': 4.},
    'fem_vs_tmm': {'p95_scaled_complex_error': .05},
}
POST_HOC_NOTES = [
    'post_hoc_peak_shift fields (peak frequencies and levels) were added after the first run; they are diagnostics, not gates.',
    'post_hoc_estimated_plane_cap_at_author_distance was added after review: the estimated probe plane is cap echo distance '
    'minus 0.135 m, so a model consistent with it puts the cap 0.135 m past the flare. It is a diagnostic, not a gate; '
    'the frozen nominal impedance gate is unchanged.',
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def radius_function(probe_to_flare):
    def radius(z):
        x = z - probe_to_flare
        if x <= 0:
            return R_THROAT
        if x <= CONE_LENGTH:
            return R_THROAT + (R_MOUTH-R_THROAT)*x/CONE_LENGTH
        return R_MOUTH
    return radius


def tmm_impedance(frequency, probe_to_flare, termination, loss=LOSS, segments=SEGMENTS,
                  flare_to_end=CONE_LENGTH+MOUTH_CYLINDER):
    """Throat impedance normalized by rho*c/S_throat.

    ``flare_to_end`` is the distance from flare start to the cap or mouth; the
    part beyond the cone is a cylinder of the mouth radius.
    """
    length = probe_to_flare + flare_to_end
    kwargs = dict(radiation_model='closed') if termination == 'closed' else \
        dict(radiation_model='finite_flange', flange_width=WALL)
    t = compute_horn_transfer_tmm(frequency, radius_function(probe_to_flare), length, R_THROAT, R_MOUTH,
                                  n_segments=segments, air=AIR, loss_model=loss, **kwargs)
    return (t['z_real']+1j*t['z_imag'])/(AIR.rho*AIR.c)


def analytic_infinite_cone(frequency, probe_to_flare):
    """Authors' lossless ideal: throat tube then an infinite cone, normalized."""
    k = 2*np.pi*frequency/AIR.c
    apex = R_THROAT*CONE_LENGTH/(R_MOUTH-R_THROAT)
    z_flare = 1/(1+1/(1j*k*apex))
    t = np.tan(k*probe_to_flare)
    return (z_flare+1j*t)/(1+1j*z_flare*t)


def reflectance(z_normalized):
    return (z_normalized-1)/(z_normalized+1)


def full_grid(df, bins):
    return np.arange(1, bins)*df


def with_dc(r):
    """Uniform 0..Nyquist spectrum: DC copies the first bin, Nyquist is zero."""
    return np.concatenate([[r[0].real], r, [0.]])


def time_window(r_positive, df, t_end, pad=4):
    """Zero R above TRUNCATION_HZ, keep the reflectance impulse up to t_end."""
    spectrum = with_dc(r_positive)
    n_bins = len(spectrum)
    spectrum[np.arange(n_bins)*df >= TRUNCATION_HZ] = 0
    n = 2*(n_bins-1)*pad
    padded = np.zeros(n//2+1, dtype=complex); padded[:n_bins] = spectrum
    impulse = np.fft.irfft(padded, n)
    rate = n*df
    impulse[int(round(t_end*rate))+1:] = 0
    return np.fft.rfft(impulse)[1:n_bins-1]


def residual_length(r_positive, df, pad=8):
    """Port of std_hrn_rfl.m residual_length: probe-to-flare from cap arrival."""
    spectrum = with_dc(r_positive)
    n = 2*(len(spectrum)-1)*pad
    padded = np.zeros(n//2+1, dtype=complex); padded[:len(spectrum)] = spectrum
    impulse = np.fft.irfft(padded, n)
    start = pad
    peak = start + int(np.argmax(np.abs(impulse[start:n-start])))
    return AIR.c*(peak/(n*df))/2 - CAP_DISTANCE_AUTHORS


def band_mask(frequency, band):
    return (frequency >= band[0]) & (frequency <= band[1])


def reflectance_metrics(predicted, measured):
    error = np.abs(np.abs(predicted)-np.abs(measured))
    return {'median_abs_reflectance_error': float(np.median(error)),
            'p95_abs_reflectance_error': float(np.percentile(error, 95)),
            'max_abs_reflectance_error': float(error.max())}


def impedance_metrics(predicted, measured):
    error = np.abs(20*np.log10(np.abs(predicted)/np.abs(measured)))
    return {'median_db': float(np.median(error)), 'p95_db': float(np.percentile(error, 95)),
            'max_db': float(error.max()),
            'p95_scaled_complex_error': float(np.percentile(np.abs(predicted-measured)/(1+np.abs(measured)), 95))}


def peak_frequencies(frequency, z, band):
    """Impedance maxima with at least 20 dB prominence inside the band."""
    from scipy.signal import find_peaks
    mask = band_mask(frequency, band)
    index, _ = find_peaks(20*np.log10(np.abs(z[mask])), prominence=20)
    return frequency[mask][index]


def peak_shift(frequency, predicted, measured):
    """Post-hoc diagnostic added after the first run; not a gate.

    Peaks are paired in order. Levels are 20 log10 abs(Z/(rho c/S)) at each peak.
    """
    m, q = peak_frequencies(frequency, measured, BAND), peak_frequencies(frequency, predicted, BAND)
    n = min(len(m), len(q))
    m, q = m[:n], q[:n]

    def level(z, f):
        return 20*np.log10(np.abs(z[np.searchsorted(frequency, f)]))
    lm, lq = level(measured, m), level(predicted, q)
    return {'measured_hz': m.tolist(), 'predicted_hz': q.tolist(),
            'relative_shift_percent': (100*(q/m-1)).tolist(),
            'measured_peak_db': lm.tolist(), 'predicted_peak_db': lq.tolist(),
            'peak_level_difference_db': (lq-lm).tolist()}


def gate(metrics, limits):
    return all(metrics[k] <= v for k, v in limits.items())


def write_frame(path, frequency, z):
    pd.DataFrame({'frequency': frequency, 'z_real_normalized': z.real, 'z_imag_normalized': z.imag}).to_csv(path, index=False)


def run_fem(out, frequency_count, mesh_size):
    import gmsh
    from horn_core.cad_units import set_step_units
    from horn_solver.solver import run_simulation_from_step
    length = NOMINAL_PROBE_TO_FLARE + CONE_LENGTH + MOUTH_CYLINDER
    step = out/'capped-horn-air.step'
    gmsh.initialize()
    set_step_units(gmsh)  # write true-size STEP (metre model)
    try:
        gmsh.model.add('hajicek_capped_horn')
        occ = gmsh.model.occ
        parts = [occ.addCylinder(0, 0, 0, 0, 0, NOMINAL_PROBE_TO_FLARE, R_THROAT),
                 occ.addCone(0, 0, NOMINAL_PROBE_TO_FLARE, 0, 0, CONE_LENGTH, R_THROAT, R_MOUTH),
                 occ.addCylinder(0, 0, NOMINAL_PROBE_TO_FLARE+CONE_LENGTH, 0, 0, MOUTH_CYLINDER, R_MOUTH)]
        occ.fuse([(3, parts[0])], [(3, p) for p in parts[1:]])
        occ.synchronize(); gmsh.write(str(step))
    finally:
        gmsh.finalize()
    csv = out/'fem-capped-nominal.csv'
    run_simulation_from_step(str(step), BAND, frequency_count, {'length': length}, str(csv), BAND[1],
                             mesh_size=mesh_size, element_degree=2, loss_model='boundary_layer',
                             minimum_wall_scale=R_THROAT, radiation_model='closed', air=AIR)
    frame = pd.read_csv(csv)
    return frame, (frame.z_real.to_numpy()+1j*frame.z_imag.to_numpy())/(AIR.rho*AIR.c)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-dir', required=True, help='Imported hajicek-neely-horn-reflectance directory')
    p.add_argument('--output-dir', required=True)
    p.add_argument('--catalog', default='data/validation/references.json')
    p.add_argument('--with-fem', action='store_true')
    p.add_argument('--fem-frequencies', type=int, default=161)
    p.add_argument('--fem-mesh', type=float, default=.004)
    a = p.parse_args()
    root, out = Path(a.reference_dir), Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((root/'manifest.json').read_text())
    catalog = next(r for r in json.loads(Path(a.catalog).read_text())['references'] if r['id'] == manifest['reference']['id'])
    if manifest['reference'] != catalog:
        raise ValueError('Imported reference disagrees with the pinned catalog')
    curves = [c for c in manifest['curves'] if c['source_member'] in FILES]
    grids = {(c['df_hz'], c['samples']) for c in curves}
    if len(grids) != 1:
        raise ValueError('Expected one shared measurement grid')
    df, samples = grids.pop()
    frequency = full_grid(df, samples+1)
    band = band_mask(frequency, BAND)

    # 1. Nominal predictions and protocol, before any measured value is read.
    predictions = {
        'tmm-closed-nominal': tmm_impedance(frequency, NOMINAL_PROBE_TO_FLARE, 'closed'),
        'tmm-closed-nominal-1200seg': tmm_impedance(frequency, NOMINAL_PROBE_TO_FLARE, 'closed', segments=2*SEGMENTS),
        'tmm-closed-nominal-boundary-layer': np.where(band, 0, np.nan)+0j,
    }
    predictions['tmm-closed-nominal-boundary-layer'][band] = tmm_impedance(frequency[band], NOMINAL_PROBE_TO_FLARE, 'closed', loss='boundary_layer')
    open_band = band_mask(frequency, OPEN_BAND)
    predictions['tmm-open-nominal'] = np.full(len(frequency), np.nan, dtype=complex)
    predictions['tmm-open-nominal'][open_band] = tmm_impedance(frequency[open_band], NOMINAL_PROBE_TO_FLARE, 'open')
    for name, z in predictions.items():
        write_frame(out/f'{name}.csv', frequency, z)
    protocol = {
        'source': manifest['reference']['source_url'], 'doi': manifest['reference']['doi'],
        'imported_manifest_sha256': sha(root/'manifest.json'),
        'air': asdict(AIR), 'temperature_c': T_C, 'loss_model': LOSS, 'segments': SEGMENTS,
        'geometry_m': {'throat_radius': R_THROAT, 'cone_length': CONE_LENGTH, 'mouth_radius': R_MOUTH,
                       'mouth_cylinder_length': MOUTH_CYLINDER, 'open_mouth_wall_flange': WALL,
                       'nominal_probe_to_flare': NOMINAL_PROBE_TO_FLARE},
        'terminations': FILES, 'channels': ['zl1', 'zl2'],
        'reference_impedance': 'nominal rho*c/(pi*throat_radius^2); not estimated from data',
        'band_hz': BAND, 'open_band_hz': OPEN_BAND, 'truncation_hz': TRUNCATION_HZ,
        'window': 'reflectance impulse kept to t = 2*(0.08 m + probe_to_flare)/c, same for model and measurement',
        'probe_plane_for_window': "authors' time-of-flight estimate, cap/mouth arrival minus 0.135 m",
        'limits': LIMITS, 'fit_parameters': [],
        'prediction_sha256': {f.name: sha(f) for f in sorted(out.glob('tmm-*.csv'))},
        'code_sha256': {f: sha(REPO/f) for f in ['scripts/validate_horn_reflectance.py',
                                                 'packages/horn-core/src/horn_core/duct.py',
                                                 'packages/horn-core/src/horn_core/webster.py']},
    }
    (out/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')

    # 2. Read measurements.
    measured = {}
    for c in curves:
        path = root/c['csv']
        if sha(path) != c['csv_sha256']:
            raise ValueError(f'Changed reference CSV: {path}')
        frame = pd.read_csv(path)
        if not np.allclose(frame.frequency, frequency, rtol=0, atol=1e-6):
            raise ValueError('Measurement grid changed')
        z = (frame.z_real_cgs_acoustic_ohm+1j*frame.z_imag_cgs_acoustic_ohm).to_numpy()*CGS_TO_SI
        measured[(c['source_member'], c['channel'])] = z/(AIR.rho*AIR.c/(np.pi*R_THROAT**2))

    z_nominal = predictions['tmm-closed-nominal']
    rows, gates = [], {}
    probe_planes = {}
    for (member, channel), zm in sorted(measured.items()):
        rm = reflectance(zm)
        xr = residual_length(rm, df)
        probe_planes[f'{member}:{channel}'] = {'estimated_m': xr, 'author_reported_m': AUTHOR_RESIDUAL_LENGTH_M[member],
                                               'reproduced_within_0_5_mm': bool(abs(xr-AUTHOR_RESIDUAL_LENGTH_M[member]) <= 5e-4)}
        row = {'source_member': member, 'channel': channel, 'termination': FILES[member],
               'estimated_probe_to_flare_m': xr}
        # Time-windowed flare reflectance: termination is outside the window.
        t_end = 2*(SWEET_SPOT+xr)/AIR.c
        z_window_model = tmm_impedance(frequency, xr, 'closed')
        r_meas_w = time_window(rm, df, t_end)
        r_model_w = time_window(reflectance(z_window_model), df, t_end)
        r_ideal_w = time_window(reflectance(analytic_infinite_cone(frequency, xr)), df, t_end)
        r_lossless_w = time_window(reflectance(tmm_impedance(frequency, xr, 'closed', loss='lossless')), df, t_end)
        windowed = reflectance_metrics(r_model_w[band], r_meas_w[band])
        row['time_windowed'] = {**windowed, 'window_end_ms': 1e3*t_end,
                                'passed': gate(windowed, LIMITS['time_windowed'])}
        row['time_windowed_nominal_probe_plane'] = reflectance_metrics(
            time_window(reflectance(z_nominal), df, 2*(SWEET_SPOT+NOMINAL_PROBE_TO_FLARE)/AIR.c)[band],
            time_window(rm, df, 2*(SWEET_SPOT+NOMINAL_PROBE_TO_FLARE)/AIR.c)[band])
        row['analytic_ideal_vs_measured_windowed'] = reflectance_metrics(r_ideal_w[band], r_meas_w[band])
        numerical = reflectance_metrics(r_lossless_w[band], r_ideal_w[band])
        row['tmm_lossless_vs_analytic_cone_windowed'] = {**numerical, 'passed': gate(numerical, LIMITS['tmm_vs_analytic_cone_windowed'])}
        if FILES[member] == 'closed':
            full = reflectance_metrics(reflectance(z_nominal[band]), rm[band])
            row['untruncated_capped'] = {**full, 'passed': gate(full, LIMITS['untruncated_capped'])}
            row['untruncated_capped_estimated_probe_plane'] = reflectance_metrics(reflectance(z_window_model[band]), rm[band])
            for label, zp in (('nominal', z_nominal), ('estimated', z_window_model)):
                m = impedance_metrics(zp[band], zm[band])
                row[f'impedance_{label}_probe_plane'] = {**m, 'passed': gate({'median': m['median_db'], 'p95': m['p95_db']}, LIMITS['impedance_magnitude_db']),
                                                         'post_hoc_peak_shift': peak_shift(frequency, zp, zm)}
            # Post-hoc, not a gate: xr was estimated as cap echo distance minus
            # CAP_DISTANCE_AUTHORS, so the model consistent with that estimate
            # puts the cap CAP_DISTANCE_AUTHORS past the flare, not 130 mm.
            z_consistent = tmm_impedance(frequency, xr, 'closed', flare_to_end=CAP_DISTANCE_AUTHORS)
            m = impedance_metrics(z_consistent[band], zm[band])
            row['post_hoc_estimated_plane_cap_at_author_distance'] = {
                'flare_to_cap_m': CAP_DISTANCE_AUTHORS, 'probe_to_cap_m': xr+CAP_DISTANCE_AUTHORS,
                **m, 'would_pass_impedance_limit': gate({'median': m['median_db'], 'p95': m['p95_db']}, LIMITS['impedance_magnitude_db']),
                'untruncated_reflectance': reflectance_metrics(reflectance(z_consistent[band]), rm[band]),
                'peak_shift': peak_shift(frequency, z_consistent, zm)}
        else:
            zo = predictions['tmm-open-nominal'][open_band]
            row['open_mouth_diagnostic'] = {'band_hz': OPEN_BAND,
                                            **reflectance_metrics(reflectance(zo), rm[open_band]),
                                            **impedance_metrics(zo, zm[open_band])}
        rows.append(row)
        pd.DataFrame({'frequency': frequency, 'measured_abs_r': np.abs(rm), 'measured_windowed_abs_r': np.abs(r_meas_w),
                      'tmm_windowed_abs_r': np.abs(r_model_w), 'analytic_cone_windowed_abs_r': np.abs(r_ideal_w)}
                     ).to_csv(out/f'reflectance-{member}-{channel}.csv', index=False)

    seg = impedance_metrics(predictions['tmm-closed-nominal'][band], predictions['tmm-closed-nominal-1200seg'][band])
    bl = impedance_metrics(predictions['tmm-closed-nominal-boundary-layer'][band], predictions['tmm-closed-nominal'][band])
    fem = None
    if a.with_fem:
        frame, zf = run_fem(out, a.fem_frequencies, a.fem_mesh)
        ff = frame.frequency.to_numpy()
        zt = tmm_impedance(ff, NOMINAL_PROBE_TO_FLARE, 'closed', loss='boundary_layer')
        vs_tmm = impedance_metrics(zf, zt)
        fem = {'frequencies': len(ff), 'mesh_size_m': a.fem_mesh, 'element_degree': 2, 'loss_model': 'boundary_layer',
               'mesh_cells': int(frame.mesh_cells.iloc[0]) if 'mesh_cells' in frame else None,
               'vs_same_model_tmm': {**vs_tmm, 'passed': vs_tmm['p95_scaled_complex_error'] <= LIMITS['fem_vs_tmm']['p95_scaled_complex_error']},
               'vs_measurement': {}}
        for (member, channel), zm in sorted(measured.items()):
            if FILES[member] != 'closed':
                continue
            zi = np.interp(ff, frequency, zm.real)+1j*np.interp(ff, frequency, zm.imag)
            r = reflectance_metrics(reflectance(zf), reflectance(zi))
            fem['vs_measurement'][f'{member}:{channel}'] = {'untruncated_capped': {**r, 'passed': gate(r, LIMITS['untruncated_capped'])},
                                                           'impedance': impedance_metrics(zf, zi)}
        write_frame(out/'fem-capped-nominal-normalized.csv', ff, zf)

    result = {
        'protocol': protocol, 'probe_plane_estimates': probe_planes,
        'numerical_checks': {'tmm_600_vs_1200_segments': seg, 'kirchhoff_vs_boundary_layer_tmm': bl},
        'comparisons': rows, 'fem': fem,
        'gate_summary': {
            'time_windowed_passes': f"{sum(r['time_windowed']['passed'] for r in rows)}/{len(rows)}",
            'untruncated_capped_passes': f"{sum(r['untruncated_capped']['passed'] for r in rows if 'untruncated_capped' in r)}/{sum('untruncated_capped' in r for r in rows)}",
            'impedance_nominal_probe_plane_passes': f"{sum(r['impedance_nominal_probe_plane']['passed'] for r in rows if 'impedance_nominal_probe_plane' in r)}/{sum('impedance_nominal_probe_plane' in r for r in rows)}",
            'tmm_vs_analytic_passes': f"{sum(r['tmm_lossless_vs_analytic_cone_windowed']['passed'] for r in rows)}/{len(rows)}"},
        'post_hoc_notes': POST_HOC_NOTES,
        'physical_horn_workflow_validated': False,
        'scope': 'Rigid conical duct input impedance seen by an ear-canal probe. Does not test radiation into free space, drivers or SPL.',
    }
    (out/'validation.json').write_text(json.dumps(result, indent=2)+'\n')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    f_k = frequency/1e3
    capped = measured[('22E27P02.THL', 'zl1')]
    ax = axes[0, 0]
    ax.semilogx(f_k, np.abs(reflectance(capped)), label='Measured, capped (P02 zl1)')
    ax.semilogx(f_k, np.abs(reflectance(z_nominal)), '--', label='TMM, nominal geometry')
    if fem:
        ax.semilogx(ff/1e3, np.abs(reflectance(zf)), 'k.', ms=3, label='Production FEM')
    ax.set_title('Capped horn, untruncated |R|'); ax.set_ylabel('|R|')
    ax = axes[0, 1]
    ax.semilogx(f_k, 20*np.log10(np.abs(capped)), label='Measured, capped (P02 zl1)')
    ax.semilogx(f_k, 20*np.log10(np.abs(z_nominal)), '--', label='TMM, nominal probe plane')
    ax.set_title('Capped horn input impedance'); ax.set_ylabel('20 log10 |Z/(rho c/S)| (dB)')
    for ax, member in zip(axes[1], FILES):
        d = pd.read_csv(out/f'reflectance-{member}-zl1.csv')
        ax.semilogx(d.frequency/1e3, d.measured_windowed_abs_r, label='Measured, windowed')
        ax.semilogx(d.frequency/1e3, d.tmm_windowed_abs_r, '--', label='TMM, windowed')
        ax.semilogx(d.frequency/1e3, d.analytic_cone_windowed_abs_r, ':', label='Analytic infinite cone')
        ax.set_title(f'{member} zl1 ({FILES[member]}), flare reflectance'); ax.set_ylabel('|R|')
    for ax in axes.ravel():
        ax.set_xlim(BAND[0]/1e3, BAND[1]/1e3); ax.grid(alpha=.2); ax.set_xlabel('Frequency (kHz)'); ax.legend(fontsize=7)
    fig.text(.5, .005, 'Hajicek, Harris, Norgaard, Siegel and Neely; Zenodo 19710026 (CC BY 4.0). No fitted parameters.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .02, 1, 1)); fig.savefig(out/'comparison.png', dpi=150); plt.close(fig)
    print(json.dumps(result['gate_summary'], indent=2))


if __name__ == '__main__':
    main()
