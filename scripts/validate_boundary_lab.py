#!/usr/bin/env python3
"""Reproduce an independent interior-FEM and ideal motor comparison.

Stages: prepare (native gmsh), horn (solver container), reference (pinned native
Boundary Lab), compare (numpy/pandas/scipy/meshio + horn-analysis).
No upstream source or fixture is copied into this repository.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

import numpy as np
from boundary_lab_fixture import case_definitions, frequencies, prepare_cases, write_json

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / 'data/validation/boundary_lab_protocol.json'
INPUT_NAMES = ('horn.step', 'horn.msh', 'project.blab.json', 'request.json', 'definition.json')


def bind_checkout_imports():
    """Use only the implementation covered by the source evidence."""
    packages = ('horn_core', 'horn_drivers', 'horn_analysis', 'horn_solver')
    roots = {name: ROOT/'packages'/name.replace('_', '-')/'src' for name in packages}
    # Fail if a caller has already imported another checkout or installed wheel.
    # Silently replacing loaded modules can retain stale classes/functions.
    for name, module in tuple(sys.modules.items()):
        package = name.split('.')[0]
        if package in roots:
            path = getattr(module, '__file__', None)
            if path is None or not Path(path).resolve().is_relative_to(roots[package].resolve()):
                raise ValueError(f'{name} was imported outside the hashed checkout')
    paths = [str(path) for path in roots.values()]
    sys.path[:] = paths + [path for path in sys.path if path not in paths]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_identity():
    paths = [p for directory in ('packages/horn-core/src', 'packages/horn-solver/src',
                                'packages/horn-analysis/src', 'packages/horn-drivers/src')
             for p in (ROOT/directory).rglob('*.py')]
    paths += [Path(__file__), ROOT/'scripts/boundary_lab_fixture.py', PROTOCOL]
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}


def protocol_inputs(root, protocol):
    return {f'{name}/{file}': sha(root/name/file)
            for name, _, _ in case_definitions(protocol) for file in INPUT_NAMES}


def verify_inputs(root):
    protocol = json.loads(PROTOCOL.read_text())
    frozen = json.loads((root/'protocol.json').read_text())
    if frozen['definition'] != protocol or frozen['inputs'] != protocol_inputs(root, protocol):
        raise ValueError('Inputs or protocol differ from the frozen reference definition')
    return protocol, frozen


def prepare(root):
    root.mkdir(parents=True, exist_ok=False)
    protocol = json.loads(PROTOCOL.read_text())
    prepare_cases(root, protocol)
    write_json(root/'protocol.json', dict(definition=protocol,
        prepared_at_utc=datetime.now(timezone.utc).isoformat(), inputs=protocol_inputs(root, protocol),
        source=source_identity()))


def verify_source(frozen):
    if frozen['source'] != source_identity():
        raise ValueError('Source changed since preparation; prepare a fresh run')


def seal_stage(root, stage, files, extra=None):
    write_json(root/f'{stage}-evidence.json', dict(
        completed_at_utc=datetime.now(timezone.utc).isoformat(),
        files={str(p.relative_to(root)): sha(p) for p in sorted(files)}, **(extra or {})))


def verify_stage(root, stage):
    evidence = json.loads((root/f'{stage}-evidence.json').read_text())
    if not evidence['files'] or any(sha(root/p) != digest for p, digest in evidence['files'].items()):
        raise ValueError(f'{stage} output evidence changed')
    return evidence


def horn(root):
    bind_checkout_imports()
    from dolfinx.io import gmshio
    from mpi4py import MPI
    from horn_core.duct import AirProperties
    from horn_solver.solver import run_simulation
    protocol, frozen = verify_inputs(root)
    verify_source(frozen)
    files = []
    for name, _, _ in case_definitions(protocol):
        out = root/name/'horn-fem.csv'
        if out.exists():
            raise FileExistsError(out)
        domain, _, facets = gmshio.read_from_msh(str(root/name/'horn.msh'), MPI.COMM_WORLD, 0, gdim=3)
        run_simulation(domain, facets,
            (protocol['frequency_min_hz'], protocol['frequency_max_hz']), protocol['frequency_count'],
            {}, str(out), bc_mode='velocity', inlet_velocity_rms=1.,
            radiation_model=protocol['radiation_model'], element_degree=protocol['element_degree'],
            air=AirProperties(rho=protocol['air_density_kg_m3'], c=protocol['sound_speed_m_s']))
        files.append(out)
    verify_inputs(root)
    verify_source(frozen)
    seal_stage(root, 'horn', files, dict(source=source_identity()))


def run_logged(command, log, timeout=600):
    """Bound the external solver and clean up its own process group on failure."""
    if os.name != 'posix':
        raise RuntimeError('External reference runner currently supports POSIX hosts only')
    with log.open('w') as stream:
        process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
            if code:
                raise subprocess.CalledProcessError(code, command)
        finally:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass
            finally:
                # The CLI can exit before its Julia worker. Waiting for the
                # parent alone does not reap or stop the remaining group.
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()


def reference(root, checkout, python, julia):
    protocol, frozen = verify_inputs(root)
    verify_source(frozen)
    checkout, python, julia = checkout.resolve(), python.absolute(), julia.resolve()
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    dirty = subprocess.check_output(['git', '-C', str(checkout), 'status', '--porcelain', '--untracked-files=no'], text=True)
    if revision != protocol['upstream_revision'] or dirty:
        raise ValueError('Reference checkout must be clean and pinned')
    probe = json.loads(subprocess.check_output([str(python), '-c',
        'import blab,json,sys; print(json.dumps({"version":list(sys.version_info[:2]),"module":blab.__file__}))'], text=True))
    if probe['version'] != protocol['python_minor'] or not Path(probe['module']).resolve().is_relative_to(checkout):
        raise ValueError('Python must import the pinned checkout with the declared version')
    version = subprocess.check_output([str(julia), '--version'], text=True).strip()
    if version != 'julia version '+protocol['julia_version']:
        raise ValueError('Unexpected Julia version')
    installed = subprocess.check_output([str(python), '-m', 'pip', 'freeze'], text=True)
    files = []
    for name, _, _ in case_definitions(protocol):
        directory = root/name
        output = directory/'boundary-lab'
        if output.exists():
            raise FileExistsError(output)
        run_logged([str(python), '-c', 'import sys; from blab.cli import main; sys.exit(main())', 'project', 'solve',
                    str(directory/'project.blab.json'), '--request', str(directory/'request.json'),
                    '--backend', 'beat_cpu', '--output', str(output), '--events', 'ndjson',
                    '--julia-executable', str(julia), '--julia-threads', '1'], directory/'boundary-lab.log')
        files.extend(p for p in output.rglob('*') if p.is_file())
        print('Completed independent reference:', name, flush=True)
    verify_inputs(root)
    verify_source(frozen)
    if subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip() != revision or subprocess.check_output(['git', '-C', str(checkout), 'diff', 'HEAD'], text=True):
        raise ValueError('Upstream source changed during the solve')
    seal_stage(root, 'reference', files, dict(runtime=dict(
        revision=revision, python=probe['version'], julia=version, packages=installed)))


def require_grid(actual, expected):
    actual = np.asarray(actual)
    if actual.shape != expected.shape or not np.isfinite(actual).all() or not np.allclose(actual, expected, rtol=1e-13, atol=0):
        raise ValueError('Incomplete or mismatched frequency grid')


def relative_error(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != expected.shape or actual.size == 0 or not np.isfinite(actual).all() or not np.isfinite(expected).all():
        raise ValueError('Comparison needs matching finite nonempty quantities')
    return float(np.max(np.abs(actual-expected)/np.maximum(np.abs(expected), 1e-30)))


def read_reference(directory, protocol):
    import meshio
    from scipy.spatial import cKDTree
    result = directory/'boundary-lab'
    manifest = json.loads((result/'manifest.json').read_text())
    expected = frequencies(protocol)
    require_grid(manifest['frequencies_hz'], expected)
    if (manifest['status'] != 'complete' or manifest['completion_mask'] != [True]*len(expected)
            or manifest['phasor_convention'] != 'exp(-i omega t)' or manifest['solve_kind'] != 'interior_fem'
            or manifest['excitation_port_ids'] != ['voltage']):
        raise ValueError('Independent solve is incomplete or uses different conventions')
    if len(manifest['meshes']) != 1 or manifest['meshes'][0]['sha256'] != sha(directory/'horn.msh'):
        raise ValueError('Independent solve used a different mesh')
    mesh = meshio.read(directory/'horn.msh')
    domains = json.loads((result/'domains.json').read_text())['domains']
    with np.load(result/'domains.npz', allow_pickle=False) as data:
        acoustic = next(d for d in domains if d['id'] == 'domain:fem-volume')
        motor = next(d for d in domains if d['id'] == 'components:electrodynamic-transducers')
        coords = data[acoustic['coordinates']['points_m']]
        area = float(data[motor['coordinates']['effective_area_m2']][0])
    if not np.isfinite(area) or area <= 0 or not np.isfinite(coords).all():
        raise ValueError('Invalid independent geometry')
    distance, mapping = cKDTree(coords).query(mesh.points)
    if max(distance) > 1e-10 or len(np.unique(mapping)) != len(mapping):
        raise ValueError('Independent nodes do not match the shared mesh')
    faces = mesh.cells_dict['triangle']
    tags = mesh.cell_data_dict['gmsh:physical']['triangle']
    def average(pressure, tag):
        triangles = faces[tags == tag]
        xyz = mesh.points[triangles]
        areas = np.linalg.norm(np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), axis=1)/2
        if not len(areas) or min(areas) <= 0:
            raise ValueError('Invalid boundary triangles')
        if tag == 2 and not np.isclose(sum(areas), area, rtol=1e-12, atol=0):
            raise ValueError('Piston area differs from its mesh boundary')
        return np.sum(areas*np.mean(pressure[mapping[triangles]], axis=1))/sum(areas)
    fields = {key: [] for key in ('z', 'mouth_per_velocity', 'velocity', 'current', 'inlet_pressure', 'mouth_pressure')}
    residuals = []
    for index, frequency in enumerate(expected):
        meta = json.loads((result/f'frequencies/{index:06d}.json').read_text())
        require_grid([meta['freq_hz']], np.array([frequency]))
        if meta['diagnostics']['transducer_reference_voltage_v'] != protocol['voltage_rms']:
            raise ValueError('Unexpected reference voltage')
        residual = meta['diagnostics']['relative_residual']
        if not np.isfinite(residual) or not 0 <= residual <= protocol['residual_limit']:
            raise ValueError('Independent linear solve failed its residual gate')
        residuals.append(residual)
        # Known fixed artifact paths, never a path supplied by result metadata.
        with np.load(result/f'frequencies/{index:06d}.npz', allow_pickle=False) as arrays:
            quantities = {q['quantity']: q for q in meta['quantities']}
            values = {}
            for quantity, unit, shape in [('fem_nodal_pressure', 'Pa', (1, len(coords))),
                                         ('diaphragm_velocity', 'm/s', (1, 1)),
                                         ('voice_coil_current', 'A', (1, 1))]:
                q = quantities[quantity]
                value = arrays[q['key']]
                if q['unit'] != unit or value.shape != shape or not np.isfinite(value).all():
                    raise ValueError('Independent quantity has invalid units, shape or values')
                values[quantity] = np.conj(value[0])  # exp(-i wt) -> exp(+i wt)
        velocity, current = values['diaphragm_velocity'][0], values['voice_coil_current'][0]
        if abs(velocity) <= 1e-30:
            raise ValueError('Cannot infer impedance from zero piston velocity')
        inlet = average(values['fem_nodal_pressure'], 2)
        mouth = average(values['fem_nodal_pressure'], 3)
        for key, value in [('z', inlet/velocity), ('mouth_per_velocity', mouth/velocity),
                           ('velocity', velocity), ('current', current),
                           ('inlet_pressure', inlet), ('mouth_pressure', mouth)]:
            fields[key].append(value)
    return area, {key: np.array(value) for key, value in fields.items()}, max(residuals)


def compare(root):
    bind_checkout_imports()
    import pandas as pd
    from horn_core.parameters import DriverParameters
    from horn_analysis.transfer_function import compute_driver_operating_point
    protocol, frozen = verify_inputs(root)
    verify_source(frozen)
    horn_evidence, reference_evidence = verify_stage(root, 'horn'), verify_stage(root, 'reference')
    rows = []
    for name, _, _ in case_definitions(protocol):
        directory = root/name
        frame = pd.read_csv(directory/'horn-fem.csv')
        frequency = frequencies(protocol)
        require_grid(frame.frequency, frequency)
        if (not np.isfinite(frame.select_dtypes(include='number')).all().all()
                or (frame.relative_residual > protocol['residual_limit']).any()
                or (frame.relative_residual < 0).any() or (frame.converged_reason <= 0).any()
                or set(frame.bc_mode) != {'velocity'} or set(frame.radiation_model) != {'plane_wave'}
                or set(frame.loss_model) != {'lossless'}
                or set(frame.air_c_m_s) != {protocol['sound_speed_m_s']}
                or set(frame.air_rho_kg_m3) != {protocol['air_density_kg_m3']}
                or set(frame.element_degree) != {protocol['element_degree']}):
            raise ValueError('Production solve failed health or model checks')
        area, fields, residual = read_reference(directory, protocol)
        np.testing.assert_allclose(frame.inlet_area_m2, area, rtol=1e-12, atol=0)
        np.testing.assert_allclose(frame.inlet_u_real, area, rtol=1e-12, atol=0)
        np.testing.assert_allclose(frame.inlet_u_imag, 0, atol=1e-15)
        d = protocol['driver']
        driver = DriverParameters(driver_id='synthetic-reference', manufacturer='Synthetic',
            model_name='Characterized piston', fs_hz=1/(2*np.pi*np.sqrt(d['mmd_kg']*d['cms_m_per_n'])),
            re_ohm=d['re_ohm'], le_h=d['le_h'], bl_tm=d['bl_n_per_a'], sd_m2=area,
            mms_kg=d['mmd_kg'], cms_m_per_n=d['cms_m_per_n'], rms_kg_per_s=d['rms_n_s_per_m'],
            mmd_kg=d['mmd_kg'], rear_load_mass_kg=0.)
        op = compute_driver_operating_point(driver, frequency, frame.z_real.to_numpy(),
            frame.z_imag.to_numpy(), area, protocol['voltage_rms'])
        mouth = frame.mouth_p_real.to_numpy()+1j*frame.mouth_p_imag.to_numpy()
        pairs = dict(acoustic_impedance=(frame.z_real.to_numpy()+1j*frame.z_imag.to_numpy(), fields['z']),
            mouth_pressure_per_velocity=(mouth, fields['mouth_per_velocity']),
            driver_current=(op['current_rms'], fields['current']),
            driver_velocity=(op['velocity_rms'], fields['velocity']),
            driver_throat_pressure=(op['throat_pressure'], fields['inlet_pressure']),
            coupled_mouth_pressure=(mouth*op['velocity_rms'], fields['mouth_pressure']))
        errors = {key: relative_error(a, b) for key, (a, b) in pairs.items()}
        rows.append(dict(case=name, frequencies=len(frequency), mesh_area_m2=area,
            max_complex_relative_errors=errors, max_upstream_residual=residual,
            passed=max(errors.values()) <= protocol['relative_complex_error_limit'],
            reference_values={key: [[float(x.real), float(x.imag)] for x in value]
                              for key, value in fields.items()}))
    # Recheck after reading arrays to avoid publishing mixed or modified evidence.
    verify_inputs(root)
    verify_stage(root, 'horn')
    verify_stage(root, 'reference')
    verify_source(frozen)
    output = dict(schema_version=1, scope=protocol['scope'], definition=protocol,
        protocol_sha256=sha(root/'protocol.json'), source=frozen['source'],
        horn_evidence_sha256=sha(root/'horn-evidence.json'),
        reference_evidence_sha256=sha(root/'reference-evidence.json'),
        runtime=reference_evidence['runtime'], cases=rows,
        passed=all(row['passed'] for row in rows))
    destination = root/'comparison.json'
    if destination.exists():
        raise FileExistsError(destination)
    write_json(destination, output)
    print(json.dumps({k: v for k, v in output.items() if k in ('scope', 'passed')}, indent=2))
    if not output['passed']:
        raise RuntimeError('Independent comparison failed its frozen limits')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'horn', 'reference', 'compare'])
    parser.add_argument('output', type=Path)
    parser.add_argument('--checkout', type=Path)
    parser.add_argument('--python', type=Path)
    parser.add_argument('--julia', type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    if args.stage == 'reference':
        if not all((args.checkout, args.python, args.julia)):
            parser.error('Reference stage requires --checkout, --python and --julia')
        reference(root, args.checkout, args.python, args.julia)
    else:
        {'prepare': prepare, 'horn': horn, 'compare': compare}[args.stage](root)


if __name__ == '__main__':
    main()
