"""Diagnose current pipe predictions against measured Ernoult brass cylinders.

Run inside horn-solver with the repository mounted at /workspace. This deliberately
compares the existing lossless/unflanged approximation with real 2 mm flange data;
it is a model-gap diagnostic, never an exact-assembly validation pass. No model
parameters are fitted. Reference temperature is handled by a declared frequency
rescaling which preserves kL; it does not supply missing thermoviscous losses.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from horn_core.acoustics import C0, RHO0
from horn_core.webster import compute_horn_transfer_tmm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    source = Path(args.reference_dir)
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((source/'manifest.json').read_text())
    constants = manifest['reference']['conditions']
    ref_c = constants['c_m_s']
    length, radius = constants['length_m'], constants['cylinder_inner_diameter_m']/2
    curves = [c for c in manifest['curves'] if c['reference_kind'] == 'measurement' and c['configuration'] == 'Brass_O' and not c.get('duplicate_of')]
    if not curves:
        raise ValueError('No measured brass open-cylinder curves')
    # Freeze comparison settings before reading response values.
    low, high, samples = 110., 3900., 20001
    frequency = np.geomspace(low, high, samples)
    model = compute_horn_transfer_tmm(frequency*C0/ref_c, lambda _: radius,
                                     length, radius, radius, radiation_model='unflanged')
    prediction = (model['z_real'] + 1j*model['z_imag'])/(RHO0*C0)
    pd.DataFrame({'frequency': frequency, 'z_real_normalized': prediction.real,
                  'z_imag_normalized': prediction.imag}).to_csv(out/'tmm.csv', index=False)

    import gmsh
    from horn_solver.solver import run_simulation_from_step
    step = out/'cylinder-air-volume.step'
    gmsh.initialize()
    try:
        gmsh.model.add('benchmark_cylinder')
        gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, length, radius)
        gmsh.model.occ.synchronize(); gmsh.write(str(step))
    finally:
        gmsh.finalize()
    fem_runs = []
    for mesh in (.004, .002):
        raw = out/f'fem-{mesh}.csv'
        run_simulation_from_step(str(step), (low*C0/ref_c, high*C0/ref_c), 241,
                                 {'length': length}, str(raw), high*C0/ref_c,
                                 mesh_size=mesh, radiation_model='unflanged_piston')
        frame = pd.read_csv(raw)
        freq = frame.frequency.to_numpy()*ref_c/C0
        z = (frame.z_real.to_numpy()+1j*frame.z_imag.to_numpy())/(RHO0*C0)
        analytic = compute_horn_transfer_tmm(frame.frequency.to_numpy(), lambda _: radius,
                                            length, radius, radius, radiation_model='unflanged')
        analytic_z = (analytic['z_real']+1j*analytic['z_imag'])/(RHO0*C0)
        error = np.abs(z-analytic_z)/(1+np.abs(analytic_z))
        fem_runs.append({'mesh_m': mesh, 'p95_scaled_complex_error_vs_tmm': float(np.percentile(error,95)),
                         'max_scaled_complex_error_vs_tmm': float(error.max())})
        pd.DataFrame({'frequency': freq, 'z_real_normalized': z.real,
                      'z_imag_normalized': z.imag}).to_csv(out/f'fem-{mesh}-reference-temperature.csv', index=False)

    fem_frequency, fem_impedance = freq.copy(), z.copy()
    rows, grouped = [], defaultdict(list)
    for curve in curves:
        path = source/curve['csv']
        if hashlib.sha256(path.read_bytes()).hexdigest() != curve['csv_sha256']:
            raise ValueError(f'Reference CSV checksum changed: {path}')
        frame = pd.read_csv(path)
        if frame.frequency.iloc[0] > low or frame.frequency.iloc[-1] < high:
            raise ValueError(f'Reference does not bracket fixed comparison interval: {path}')
        z = np.interp(frequency, frame.frequency, frame.z_real_normalized) + 1j*np.interp(frequency, frame.frequency, frame.z_imag_normalized)
        magnitude_error = np.abs(20*np.log10(np.maximum(np.abs(prediction),1e-15)/np.maximum(np.abs(z),1e-15)))
        complex_error = np.abs(prediction-z)/(1+np.abs(z))
        resonance = (frequency >= 300) & (frequency <= 600)
        first_peak = float(frequency[resonance][np.argmax(np.abs(z[resonance]))])
        row = {'source_member': curve['source_member'], 'operator': curve['operator'],
               'median_magnitude_error_db': float(np.median(magnitude_error)),
               'p95_magnitude_error_db': float(np.percentile(magnitude_error,95)),
               'max_magnitude_error_db': float(magnitude_error.max()),
               'p95_scaled_complex_error': float(np.percentile(complex_error,95)),
               'first_peak_frequency_hz': first_peak,
               'first_peak_height_error_db': float(20*np.log10(np.max(np.abs(prediction[resonance]))/np.max(np.abs(z[resonance])))),
               'magnitude_agreement_gate': bool(np.median(magnitude_error)<=2 and np.percentile(magnitude_error,95)<=4)}
        rows.append(row); grouped[curve['operator']].append(z)
    pd.DataFrame(rows).to_csv(out/'measurement-comparisons.csv', index=False)
    result = {'reference_url': manifest['reference']['source_url'],
              'reference_archive_sha256': manifest['reference']['sha256'],
              'reference_manifest_sha256': hashlib.sha256((source/'manifest.json').read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'status': 'model_gap_diagnostic', 'physical_validation_passed': False,
              'frequency_range_hz': [low, high], 'log_frequency_samples': samples,
              'measurement_curves': len(rows), 'operators': len(grouped),
              'frequency_mapping': 'solver_frequency = reference_frequency * 343 / 346.28592; normalized Z = specific Z / (rho_solver*c_solver)',
              'fit_parameters': [],
              'known_model_mismatches': ['No thermoviscous wall losses', 'Unflanged approximation versus measured 2 mm finite flange'],
              'magnitude_gate': {'median_db_max': 2, 'p95_db_max': 4, 'passing_curves': sum(r['magnitude_agreement_gate'] for r in rows)},
              'fem_vs_same_model_tmm': fem_runs,
              'comparisons': rows}
    (out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(12,4.5))
    for ax in axes:
        for operator, values in sorted(grouped.items()):
            median = np.median(np.abs(values),axis=0)
            ax.plot(frequency,20*np.log10(median),label=f'{operator}: median of {len(values)} measurements',lw=1.2)
        ax.plot(frequency,20*np.log10(np.abs(prediction)),'k--',label='Current lossless/unflanged TMM',lw=1.3)
        ax.plot(fem_frequency,20*np.log10(np.abs(fem_impedance)),'k.',ms=2,label='Production FEM, 2 mm mesh')
        ax.set_xlabel('Frequency at reference 25°C (Hz)'); ax.grid(True,alpha=.2)
        ax.set_ylabel('20 log10 |Z / (ρc/S)| (dB)')
    axes[0].set_xscale('log'); axes[0].set_xlim(low,high)
    axes[1].set_xlim(400,520); axes[1].set_ylim(5,65)
    axes[0].legend(fontsize=7)
    fig.suptitle('180 mm × 14 mm brass pipe: measured impedance versus current approximation')
    fig.text(.5,.01,'Ernoult et al., Zenodo 20024938 (CC BY 4.0) • No fitted parameters • Diagnostic, not a physical-validation pass',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.035,1,.94)); fig.savefig(out/'comparison.png',dpi=170); plt.close(fig)
    print(json.dumps({k:v for k,v in result.items() if k!='comparisons'},indent=2))


if __name__ == '__main__':
    main()
