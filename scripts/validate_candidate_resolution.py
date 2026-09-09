#!/usr/bin/env python3
"""Freeze, solve and check a generated candidate's numerical resolution.

This study does not certify a physical assembly or improve driver evidence.
Run prepare on the host, then solve in the source-mounted solver container.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    'mesh_10': dict(mesh_size=.010, sections=20, points=101),
    'mesh_6': dict(mesh_size=.006, sections=20, points=101),
    'mesh_4': dict(mesh_size=.004, sections=20, points=101),
    'loft_40': dict(mesh_size=.004, sections=40, points=101),
    'loft_80': dict(mesh_size=.004, sections=80, points=101),
    'frequency_201': dict(mesh_size=.004, sections=80, points=201),
}
LIMITS = dict(mesh_spl_db=.5, mesh_impedance_db=.5, mesh_phase_deg=5.,
              loft_spl_db=.5, loft_impedance_db=.5, loft_phase_deg=5.,
              frequency_spl_db=.2, frequency_impedance_db=.2,
              frequency_phase_deg=2., frequency_ripple_db=.2)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def source_identity():
    paths = list((ROOT/'packages').glob('*/src/**/*.py')) + [Path(__file__)]
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}


def bind_source():
    roots = {p.name.replace('-', '_'): p/'src' for p in (ROOT/'packages').glob('horn-*')}
    for name, module in tuple(sys.modules.items()):
        package = name.split('.')[0]
        if package in roots:
            path = getattr(module, '__file__', None)
            if path is None or not Path(path).resolve().is_relative_to(roots[package].resolve()):
                raise ValueError(f'{name} was imported outside the frozen checkout')
    paths = [str(p) for p in roots.values()]
    sys.path[:] = paths + [p for p in sys.path if p not in paths]


def prepare(out, ranking, driver, low, high, index):
    if not np.isfinite([low, high]).all() or not 0 < low < high:
        raise ValueError('A positive ordered target band is required')
    candidates = json.loads(ranking.read_text())
    if isinstance(candidates, dict):
        candidates = candidates['results']
    if index < 0:
        raise ValueError('Candidate index must be nonnegative')
    candidate = candidates[index]
    record = json.loads(driver.read_text())
    if record['driver_id'] != candidate['driver_id']:
        raise ValueError('Driver record does not match the selected candidate')
    if (candidate.get('loss_model') != 'lossless'
            or candidate.get('radiation_model') != 'flanged_piston'
            or candidate.get('element_degree') != 1):
        raise ValueError('This protocol supports the lossless flanged-piston P1 model only')
    for key in ('throat_radius', 'mouth_radius', 'length', 'drive_voltage_rms', 'observation_distance_m'):
        if not np.isfinite(candidate[key]) or candidate[key] <= 0:
            raise ValueError(f'Invalid candidate {key}')
    out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(ranking, out/'ranking.json')
    shutil.copyfile(driver, out/'driver.json')
    write_json(out/'protocol.json', dict(schema_version=1,
        prepared_at=datetime.now(timezone.utc).isoformat(),
        source_revision=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        source=source_identity(), inputs={p:sha(out/p) for p in ('ranking.json','driver.json')},
        candidate_index=index, candidate=candidate, target_band_hz=[low,high],
        simulation_band_hz=[low/np.sqrt(2),high*np.sqrt(2)], cases=CASES, limits=LIMITS,
        scope='Numerical resolution of this fixed ideal candidate; no physical qualification'))


def verify(out):
    p = json.loads((out/'protocol.json').read_text())
    if p['source'] != source_identity() or p['cases'] != CASES or p['limits'] != LIMITS:
        raise ValueError('Source or fixed protocol changed; prepare a fresh study')
    if any(sha(out/name) != digest for name,digest in p['inputs'].items()):
        raise ValueError('Frozen input changed')
    ranking=json.loads((out/'ranking.json').read_text())
    if isinstance(ranking,dict):
        ranking=ranking['results']
    if ranking[p['candidate_index']] != p['candidate']:
        raise ValueError('Frozen candidate changed')
    return p


def solve(out):
    bind_source()
    from horn_geometry.generator import create_horn
    from horn_solver.solver import run_simulation_from_step
    p=verify(out)
    c=p['candidate']
    for name,case in CASES.items():
        directory=out/name
        directory.mkdir(exist_ok=False)
        step=directory/'horn.step'
        create_horn(c['profile'],c['throat_radius'],c['mouth_radius'],c['length'],step,
                    num_sections=case['sections'])
        run_simulation_from_step(str(step),tuple(p['simulation_band_hz']),case['points'],
            {'length':c['length']},str(directory/'response.csv'),max(p['simulation_band_hz']),
            mesh_size=case['mesh_size'],element_degree=1,radiation_model='flanged_piston',loss_model='lossless')
        verify(out)
        print('Completed resolution case:',name,flush=True)
    files=list(out.glob('*/horn.step'))+list(out.glob('*/response.csv'))
    write_json(out/'solve-evidence.json',dict(protocol_sha256=sha(out/'protocol.json'),
        files={str(f.relative_to(out)):sha(f) for f in sorted(files)}))


def check_frame(frame, band, count):
    expected=np.geomspace(*band,count)
    if (len(frame)!=count or not np.isfinite(frame.select_dtypes(include='number')).all().all()
            or not np.allclose(frame.frequency,expected,rtol=1e-12,atol=0)):
        raise ValueError('Incomplete, nonfinite or incorrect frequency grid')
    for key,value in dict(schema_version=2,bc_mode='dirichlet',phasor_convention='exp(+iwt)_rms',
                          radiation_model='flanged_piston',loss_model='lossless',element_degree=1).items():
        if not (frame[key]==value).all():
            raise ValueError(f'Unexpected model contract: {key}')
    incoming=frame.input_acoustic_power_w.to_numpy()
    outgoing=(frame.mouth_acoustic_power_w+frame.viscous_wall_power_w+frame.thermal_wall_power_w).to_numpy()
    balance=np.abs(incoming-outgoing)/np.maximum(np.abs(incoming),1e-30)
    if ((frame.relative_residual<0).any() or (frame.relative_residual>1e-8).any()
            or (frame.converged_reason<=0).any() or (frame.z_real<0).any()
            or (incoming<0).any() or (outgoing<0).any() or np.max(balance)>1e-7):
        raise ValueError('Failed numerical health or power balance')
    return dict(mesh_cells=int(frame.mesh_cells.iloc[0]),max_residual=float(frame.relative_residual.max()),
                max_relative_power_balance_error=float(np.max(balance)))


def curve_change(coarse, fine):
    """Compare on the entire finer grid; interpolation cannot drop extra peaks."""
    f=fine['frequency']
    if coarse['frequency'][0]>f[0] or coarse['frequency'][-1]<f[-1]:
        raise ValueError('Comparison grids do not cover the same band')
    level=np.interp(f,coarse['frequency'],coarse['spl'])
    z=(np.interp(f,coarse['frequency'],coarse['z'].real)
       +1j*np.interp(f,coarse['frequency'],coarse['z'].imag))
    if np.any(np.abs(z)<=0) or np.any(np.abs(fine['z'])<=0):
        raise ValueError('Zero impedance cannot support a magnitude/phase comparison')
    return dict(spl_db=float(np.max(np.abs(level-fine['spl']))),
        impedance_db=float(np.max(np.abs(20*np.log10(np.abs(z)/np.abs(fine['z']))))),
        phase_deg=float(np.max(np.abs(np.angle(z*np.conj(fine['z']),deg=True)))))


def ripple(curve, band):
    f=curve['frequency']
    if f[0]>band[0] or f[-1]<band[1]:
        raise ValueError('Target band not covered')
    points=np.r_[band[0],f[(f>band[0]) & (f<band[1])],band[1]]
    return float(np.ptp(np.interp(points,f,curve['spl'])))


def compare(out):
    bind_source()
    import pandas as pd
    from horn_drivers.loader import _driver_from_dict
    from horn_analysis.evaluation import coupled_output
    p=verify(out)
    evidence=json.loads((out/'solve-evidence.json').read_text())
    expected={f'{name}/{file}' for name in CASES for file in ('horn.step','response.csv')}
    if (set(evidence['files'])!=expected or evidence['protocol_sha256']!=sha(out/'protocol.json')
            or any(sha(out/name)!=digest for name,digest in evidence['files'].items())):
        raise ValueError('Solve evidence changed or is incomplete')
    candidate=p['candidate']
    driver=_driver_from_dict(json.loads((out/'driver.json').read_text()))
    target=SimpleNamespace(voltage_rms=candidate['drive_voltage_rms'],
                           observation_distance_m=candidate['observation_distance_m'])
    curves,health={},{}
    for name,case in CASES.items():
        frame=pd.read_csv(out/name/'response.csv')
        health[name]=check_frame(frame,p['simulation_band_hz'],case['points'])
        levels,_,_=coupled_output(frame,driver,target)
        curves[name]=dict(frequency=frame.frequency.to_numpy(),spl=levels,
                         z=frame.z_real.to_numpy()+1j*frame.z_imag.to_numpy())
    results=[]
    for kind,a,b in [('mesh','mesh_10','mesh_6'),('mesh','mesh_6','mesh_4'),
                      ('loft','mesh_4','loft_40'),('loft','loft_40','loft_80'),
                      ('frequency','loft_80','frequency_201')]:
        change=curve_change(curves[a],curves[b])
        if kind=='frequency':
            change['ripple_db']=abs(ripple(curves[a],p['target_band_hz'])-ripple(curves[b],p['target_band_hz']))
        passed=all(np.isfinite(value) and value<=LIMITS[f'{kind}_{key}'] for key,value in change.items())
        results.append(dict(kind=kind,coarse=a,fine=b,changes=change,passed=bool(passed)))
    verify(out)
    if any(sha(out/name)!=digest for name,digest in evidence['files'].items()):
        raise ValueError('Solve evidence changed during comparison')
    destination=out/'comparison.json'
    if destination.exists():
        raise FileExistsError(destination)
    result=dict(scope=p['scope'],protocol_sha256=sha(out/'protocol.json'),
        solve_evidence_sha256=sha(out/'solve-evidence.json'),candidate=candidate,
        limits=LIMITS,health=health,comparisons=results,passed=all(r['passed'] for r in results),
        physical_validation_status='experimental_prediction')
    write_json(destination,result)
    print(json.dumps(result,indent=2))
    if not result['passed']:
        raise RuntimeError('Candidate failed fixed resolution gates; preserve results and investigate')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','solve','compare'])
    parser.add_argument('output',type=Path)
    parser.add_argument('--ranking',type=Path)
    parser.add_argument('--driver',type=Path)
    parser.add_argument('--f-low',type=float)
    parser.add_argument('--f-high',type=float)
    parser.add_argument('--candidate-index',type=int,default=0)
    args=parser.parse_args()
    out=args.output.resolve()
    if args.stage=='prepare':
        if any(v is None for v in (args.ranking,args.driver,args.f_low,args.f_high)):
            parser.error('prepare requires --ranking, --driver, --f-low and --f-high')
        prepare(out,args.ranking,args.driver,args.f_low,args.f_high,args.candidate_index)
    else:
        {'solve':solve,'compare':compare}[args.stage](out)


if __name__=='__main__':
    main()
