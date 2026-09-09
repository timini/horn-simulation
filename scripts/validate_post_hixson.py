#!/usr/bin/env python3
"""Reproduce one measured exponential-horn load comparison (Post 1994).

The blue vector traces in published Figures 4.26/4.27 are digitized measurements,
not original analyzer samples. This does not validate electrical-drive output.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_URL = 'https://audioroundtable.com/misc/post_hixson_horns.pdf'
GEOMETRY = dict(profile='exponential', throat_radius=.0254, mouth_radius=.271, length=.559)
CASES = {'mesh24': (.024, 40), 'mesh18': (.018, 40), 'mesh14': (.014, 40)}
# Freeze before solving. Dimensionless ka avoids inventing the measured air temperature.
KA = np.geomspace(1., 5., 81)
LIMITS = dict(measurement_p95_magnitude_db=2., measurement_p95_phase_deg=10.,
              measurement_max_complex_normalized=.25,
              refinement_max_magnitude_db=.5, refinement_max_phase_deg=5.)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def source_hashes():
    files = list((ROOT/'packages').glob('*/src/**/*.py'))
    files += [Path(__file__), ROOT/'scripts/validate_candidate_resolution.py']
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(files)}


def coordinates(path):
    d = path.get('d')
    if re.search(r'[A-KN-Yac-z]', d):
        raise ValueError('Expected straight vector segments only')
    return np.array([float(v) for v in re.findall(r'[-+]?\d+(?:\.\d+)?', d)]).reshape(-1, 2)


def trace(svg):
    """Extract measured blue strokes; retain only the predeclared ka 1--5 band."""
    paths = ET.parse(svg).findall('.//{http://www.w3.org/2000/svg}path')
    frame = next(p for p in paths if p.get('stroke') == 'rgb(0%, 0%, 0%)')
    box = coordinates(frame)
    low, high = box.min(axis=0), box.max(axis=0)
    points = []
    for p in paths:
        if p.get('stroke') != 'rgb(0%, 0%, 100%)':
            continue
        if p.get('transform') != frame.get('transform'):
            raise ValueError('Unexpected plot transformation')
        xy = (coordinates(p)-low)/(high-low)*[10., 1.3]
        # The blue legend is at ka > 6, outside the declared comparison band.
        points.extend(xy[(xy[:, 0] >= .95) & (xy[:, 0] <= 5.05)])
    points = np.array(points)
    points = points[np.argsort(points[:, 0])]
    xs, inv = np.unique(points[:, 0], return_inverse=True)
    ys = np.array([np.mean(points[inv == i, 1]) for i in range(len(xs))])
    if len(xs) < 100 or xs[0] > KA[0] or xs[-1] < KA[-1] or np.max(np.diff(xs)) > .04:
        raise ValueError('Incomplete measured trace')
    return np.column_stack((xs, ys))


def prepare(pdf, out, image):
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True):
        raise ValueError('Commit the protocol before preparing evidence')
    if not re.fullmatch(r'sha256:[a-f0-9]{64}', image):
        raise ValueError('Use an immutable solver image ID')
    actual = subprocess.check_output(['docker', 'image', 'inspect', image, '--format', '{{.Id}}'], text=True).strip()
    if actual != image:
        raise ValueError('Image identity mismatch')
    out.mkdir(parents=True, exist_ok=False)
    for component, page in [('resistance', 120), ('reactance', 121)]:
        svg = out/f'{component}.svg'
        subprocess.run(['pdftocairo', '-f', str(page), '-l', str(page), '-svg', str(pdf), str(svg)], check=True)
        data = trace(svg)
        np.savetxt(out/f'measured_{component}.csv', data, delimiter=',', header='ka,normalized_value', comments='')
        # Source document/figure artwork is not republished; retain numeric observations.
        svg.unlink()
    protocol = dict(source_url=SOURCE_URL, source_pdf_sha256=sha(pdf),
                    source_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    source_hashes=source_hashes(), solver_image_id=image, geometry=GEOMETRY,
                    cases=CASES, ka=KA.tolist(), rho_kg_m3=1.225, c_m_s=343., limits=LIMITS,
                    measured={p.name: sha(p) for p in out.glob('measured_*.csv')},
                    figures={'resistance': '4.26, printed p.108 / PDF p.120',
                             'reactance': '4.27, printed p.109 / PDF p.121'},
                    measurement_scope='Digitized published complex normalized horn input impedance; no volts-to-SPL validation',
                    normalization='Z_acoustic * S_throat / (rho*c); plotted area correction 1.025 already applied, never applied twice',
                    digitization_resolution_normalized=.005,
                    geometry_limit='Nominal dimensions; actual throat area +2.5%, mouth tolerance 2 mm, remaining as-built contour unknown',
                    band_rationale='ka 1--5 is above sub-100 Hz calibration concerns and below documented ka > 6 wall resonances')
    write(out/'protocol.json', protocol)
    with (out/'source.tar').open('wb') as handle:
        subprocess.run(['git', 'archive', 'HEAD'], cwd=ROOT, stdout=handle, check=True)


def verify(out):
    p = json.loads((out/'protocol.json').read_text())
    if p['source_hashes'] != source_hashes() or p['limits'] != LIMITS or p['geometry'] != GEOMETRY:
        raise ValueError('Study source or protocol changed')
    if p['ka'] != KA.tolist() or p['cases'] != {k: list(v) for k, v in CASES.items()}:
        raise ValueError('Study cases changed')
    if any(sha(out/name) != value for name, value in p['measured'].items()):
        raise ValueError('Measurement observations changed')
    return p


def worker(out, case):
    sys.path[:0] = [str(p/'src') for p in (ROOT/'packages').glob('horn-*')]
    from horn_geometry.generator import create_horn
    from horn_solver.solver import create_mesh_from_step, run_simulation
    from validate_candidate_resolution import runtime_identity
    p = verify(out)
    h, sections = CASES[case]
    step = out/f'{case}.step'
    result = out/f'{case}.csv'
    if result.exists():
        raise FileExistsError(result)
    create_horn(**GEOMETRY, output_file=step, num_sections=sections)
    domain, tags = create_mesh_from_step(str(step), h, GEOMETRY['length'])
    frequencies = KA*p['c_m_s']/(2*np.pi*GEOMETRY['mouth_radius'])
    run_simulation(domain, tags, (frequencies[0], frequencies[-1]), len(KA),
                   {'length': GEOMETRY['length']}, str(result),
                   bc_mode='velocity', radiation_model='modal_baffled')
    write(out/f'{case}_runtime.json', dict(runtime=runtime_identity(),
          mesh_cells=domain.topology.index_map(domain.topology.dim).size_global))


def solve(out):
    p = verify(out)
    with (out/'execution.json').open('x') as f:
        json.dump(dict(protocol_sha256=sha(out/'protocol.json'), image=p['solver_image_id']), f)
    for case in CASES:
        command = ['docker', 'run', '--rm',
                   '-e', 'OPENBLAS_NUM_THREADS=1', '-e', 'OMP_NUM_THREADS=1',
                   '-e', 'PYTHONPATH=/usr/local/lib', '-v', f'{ROOT}:/workspace:ro',
                   '-v', f'{out}:/study', '-w', '/workspace', p['solver_image_id'],
                   'python3', 'scripts/validate_post_hixson.py', 'worker', '--out', '/study', '--case', case]
        with (out/f'{case}.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        print(f'Completed {case}', flush=True)
    verify(out)
    write(out/'outputs.json', {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()})


def metrics(prediction, reference):
    if not np.isfinite(prediction).all() or not np.isfinite(reference).all() or np.any(abs(reference) == 0):
        raise ValueError('Invalid impedance')
    level = np.abs(20*np.log10(abs(prediction)/abs(reference)))
    phase = np.abs(np.angle(prediction*np.conj(reference), deg=True))
    return dict(p95_magnitude_db=float(np.percentile(level, 95)),
                p95_phase_deg=float(np.percentile(phase, 95)),
                max_magnitude_db=float(max(level)), max_phase_deg=float(max(phase)),
                max_complex_normalized=float(max(abs(prediction-reference))))


def compare(out):
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p = json.loads((out/'protocol.json').read_text())
    seals = json.loads((out/'outputs.json').read_text())
    if any(sha(out/name) != value for name, value in seals.items()):
        raise ValueError('Executed inputs or outputs changed')
    measured = []
    for component in ('resistance', 'reactance'):
        curve = np.loadtxt(out/f'measured_{component}.csv', delimiter=',', skiprows=1)
        measured.append(np.interp(KA, curve[:, 0], curve[:, 1]))
    reference = measured[0]+1j*measured[1]
    calculated, health = {}, {}
    for case in CASES:
        df = pd.read_csv(out/f'{case}.csv')
        if len(df) != len(KA) or not np.allclose(df.frequency*2*np.pi*.271/343, KA, rtol=1e-10):
            raise ValueError('Wrong solver frequency grid')
        calculated[case] = (df.z_real.to_numpy()+1j*df.z_imag.to_numpy())/(p['rho_kg_m3']*p['c_m_s'])
        if np.any(calculated[case].real <= 0):
            raise ValueError('Nonpassive impedance')
        balance = np.abs(df.input_acoustic_power_w-df.mouth_acoustic_power_w)/np.maximum(np.abs(df.input_acoustic_power_w), 1e-30)
        health[case] = dict(relative_residual=float(df.relative_residual.abs().max()),
                            power_balance_relative_error=float(balance.max()))
        if health[case]['relative_residual'] > 1e-8 or health[case]['power_balance_relative_error'] > 1e-7:
            raise ValueError('Numerical health failed')
    comparisons = {case: metrics(value, reference) for case, value in calculated.items()}
    refinements = {a+'_'+b: metrics(calculated[b], calculated[a]) for a, b in [('mesh24','mesh18'), ('mesh18','mesh14')]}
    final = comparisons['mesh14']
    passed = (final['p95_magnitude_db'] <= LIMITS['measurement_p95_magnitude_db']
              and final['p95_phase_deg'] <= LIMITS['measurement_p95_phase_deg']
              and final['max_complex_normalized'] <= LIMITS['measurement_max_complex_normalized'])
    converged = all(v['max_magnitude_db'] <= LIMITS['refinement_max_magnitude_db']
                    and v['max_phase_deg'] <= LIMITS['refinement_max_phase_deg'] for v in refinements.values())
    summary = dict(scope=p['measurement_scope'], limits=LIMITS, measurement_comparison_pass=bool(passed),
                   mesh_refinement_pass=bool(converged), comparisons=comparisons, refinements=refinements,
                   health=health, frequency_hz=[float(KA[0]*343/(2*np.pi*.271)),float(KA[-1]*343/(2*np.pi*.271))],
                   physical_assembly_qualified=False)
    write(out/'comparison.json', summary)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    for ax, part, label in zip(axes, ('real', 'imag'), ('Resistance', 'Reactance')):
        ax.plot(KA, getattr(reference, part), color='black', label='Published measurement (digitized)')
        for case, curve in calculated.items():
            ax.plot(KA, getattr(curve, part), label=f'Production FEM {CASES[case][0]*1000:g} mm', alpha=.8)
        ax.set_ylabel(label+' / (rho c)'); ax.grid(alpha=.2)
    axes[0].legend(fontsize=8); axes[1].set_xlabel('ka (mouth radius 271 mm)')
    fig.suptitle('Post exponential horn: measured acoustic loading\nNominal geometry; infinite baffle; no fitted correction')
    fig.savefig(out/'comparison.png', dpi=160)
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare','solve','worker','compare'])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--pdf', type=Path)
    parser.add_argument('--image')
    parser.add_argument('--case', choices=CASES)
    args = parser.parse_args()
    out = args.out.resolve()
    if args.stage == 'prepare': prepare(args.pdf, out, args.image)
    elif args.stage == 'solve': solve(out)
    elif args.stage == 'worker': worker(out, args.case)
    else: compare(out)


if __name__ == '__main__':
    main()
