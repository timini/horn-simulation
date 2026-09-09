#!/usr/bin/env python3
"""Compare one ideal horn's infinite-baffle radiation with pinned modal reference.

This is a numerical comparison, never physical driver/assembly qualification.
The separately installed GPL toolbox is mounted read-only and not redistributed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REVISION = '24d533b191f299aafe1c6aae28cba294abc8f8e2'
GEOMETRY = dict(throat_radius=.032498,mouth_radius=.08870898640584517,length=.11790625,profile='hyperbolic')
FREQUENCIES = np.geomspace(800.,1600.,33).tolist()
CASES = [dict(id='m16_s250', modes=16, sections=250),
         dict(id='m32_s250', modes=32, sections=250),
         dict(id='m64_s250', modes=64, sections=250),
         dict(id='m64_s500', modes=64, sections=500),
         dict(id='m64_s1000', modes=64, sections=1000)]
# Freeze strict engineering comparison targets; failures remain failures.
LIMITS = dict(reference_level_db=.05, reference_phase_deg=1.,
              model_level_db=.5, model_impedance_db=.5, model_phase_deg=5.)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def source():
    paths = list((ROOT/'packages').glob('*/src/**/*.py'))
    paths += [Path(__file__), ROOT/'scripts/mmm_radiation_reference.m', ROOT/'scripts/mmm_reference.Dockerfile',
              ROOT/'scripts/validate_boundary_lab.py', ROOT/'scripts/boundary_lab_fixture.py']
    return {str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}


def clean_revision():
    def git(*args):
        return subprocess.check_output(['git','-C',str(ROOT),*args],text=True).strip()
    if git('status','--porcelain','--untracked-files=all') or set(source())-set(git('ls-files').splitlines()):
        raise ValueError('A clean tracked checkout is required')
    return git('rev-parse','HEAD')


def prepare(out):
    from scipy.io import savemat
    from scipy.special import jn_zeros
    revision=clean_revision()
    out.mkdir(parents=True,exist_ok=False)
    # Exact mathematical profile of the archived worked candidate, not a fit.
    geometry=GEOMETRY
    for count in (250,500,1000):
        z=np.linspace(0,geometry['length'],count+1)
        r=geometry['throat_radius']*np.cosh(np.arccosh(geometry['mouth_radius']/geometry['throat_radius'])*z/geometry['length'])
        np.savetxt(out/f'coordinates_{count}.csv',np.c_[z,r],delimiter=',')
    savemat(out/'MMM_besselzeros.mat',{'bz':np.r_[0.,jn_zeros(1,63)].reshape(-1,1)})
    # Octave 7 lacks this scalar string predicate used by the upstream API.
    (out/'contains.m').write_text("function value = contains(text, pattern)\nassert(ischar(text) && ischar(pattern));\nvalue = ~isempty(strfind(text, pattern));\nend\n")
    shutil.copyfile(ROOT/'scripts/mmm_radiation_reference.m',out/'reference.m')
    inputs={p.name:sha(p) for p in out.iterdir()}
    write(out/'protocol.json',dict(schema_version=1,source_revision=revision,source=source(),
        upstream_revision=REVISION,geometry=geometry,air=dict(rho=1.225,c=343.),
        frequencies_hz=FREQUENCIES,observer_distance_m=1.,
        excitation='uniform throat velocity 1 m/s RMS',phasor_convention='exp(+iwt)_rms',
        cases=CASES,limits=LIMITS,inputs=inputs,
        scope='Ideal axisymmetric horn in infinite baffle; no physical assembly qualification'))


def verify(out):
    p=json.loads((out/'protocol.json').read_text())
    if {f.name for f in out.glob('*.m')}!={'reference.m','contains.m'}:
        raise ValueError('Unexpected executable input in reference directory')
    if p['geometry']!=GEOMETRY or p['frequencies_hz']!=FREQUENCIES or p['air']!=dict(rho=1.225,c=343.) or p['observer_distance_m']!=1.:
        raise ValueError('Fixed physical specification changed')
    if p['source']!=source() or p['cases']!=CASES or p['limits']!=LIMITS or p['upstream_revision']!=REVISION:
        raise ValueError('Source or protocol changed; prepare a new study')
    if any(sha(out/name)!=digest for name,digest in p['inputs'].items()):
        raise ValueError('Frozen input changed')
    return p


def reference(out, checkout, image):
    from validate_boundary_lab import reference_checkout_clean, run_logged
    clean_revision();p=verify(out)
    revision=subprocess.check_output(['git','-C',str(checkout),'rev-parse','HEAD'],text=True).strip()
    if revision!=REVISION or not reference_checkout_clean(checkout):
        raise ValueError('Reference must be the pristine pinned checkout')
    if any((out/f"{c['id']}.csv").exists() for c in CASES):
        raise FileExistsError('Reference outputs already exist')
    image_id=subprocess.check_output(['docker','image','inspect',image,'--format','{{.Id}}'],text=True).strip()
    # Pinned image identity prevents a tag changing during this execution.
    command=['docker','run','--rm','--name',f'horn-mmm-{sha(out/"protocol.json")[:12]}',
             '-v',f'{checkout}:/reference:ro','-v',f'{out}:/study','-w','/study',image_id,
             'octave','--no-gui','--quiet','reference.m']
    try:
        run_logged(command,out/'reference.log',timeout=3600)
    finally:
        # Killing a Docker client does not reliably stop its container.
        subprocess.run(['docker','rm','-f',command[4]],capture_output=True)
    if not reference_checkout_clean(checkout):
        raise ValueError('Reference source changed during execution')
    clean_revision();verify(out)
    files={f"{c['id']}.csv":sha(out/f"{c['id']}.csv") for c in CASES}
    write(out/'reference-evidence.json',dict(protocol_sha256=sha(out/'protocol.json'),
        upstream_revision=revision,image_id=image_id,files=files))


def production(out):
    # Source-mounted solver container. The host verifies Git identity separately.
    roots=[str(p/'src') for p in (ROOT/'packages').glob('horn-*')]
    for name,module in tuple(sys.modules.items()):
        if name.startswith(('horn_core','horn_geometry','horn_solver')):
            path=getattr(module,'__file__',None)
            if path is None or not Path(path).resolve().is_relative_to(ROOT.resolve()):
                raise ValueError('Production module was loaded outside this checkout')
    sys.path[:0]=roots
    from horn_geometry.generator import create_horn
    from horn_solver.solver import run_simulation_from_step
    p=verify(out);g=p['geometry']
    if (out/'production.csv').exists():
        raise FileExistsError('Production output already exists')
    create_horn(g['profile'],g['throat_radius'],g['mouth_radius'],g['length'],out/'horn.step',num_sections=80)
    run_simulation_from_step(str(out/'horn.step'),(800.,1600.),33,{'length':g['length']},
        str(out/'production.csv'),1600.,mesh_size=.004,element_degree=1,
        radiation_model='flanged_piston',loss_model='lossless',bc_mode='velocity',inlet_velocity_rms=1.)
    verify(out)
    write(out/'production-evidence.json',dict(protocol_sha256=sha(out/'protocol.json'),
        files={name:sha(out/name) for name in ('horn.step','production.csv')}))


def change(a,b):
    if not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(np.abs(a)==0) or np.any(np.abs(b)==0):
        raise ValueError('Invalid complex comparison')
    return dict(level_db=float(np.max(np.abs(20*np.log10(np.abs(a)/np.abs(b))))),
                phase_deg=float(np.max(np.abs(np.angle(a*np.conj(b),deg=True)))))


def compare(out):
    import pandas as pd
    sys.path[:0]=[str(p/'src') for p in (ROOT/'packages').glob('horn-*')]
    from horn_core.acoustics import baffled_piston_on_axis
    clean_revision();p=verify(out)
    for name,expected in [('reference',{f"{c['id']}.csv" for c in CASES}),
                          ('production',{'horn.step','production.csv'})]:
        evidence=json.loads((out/f'{name}-evidence.json').read_text())
        if set(evidence['files'])!=expected or evidence['protocol_sha256']!=sha(out/'protocol.json') or any(sha(out/f)!=h for f,h in evidence['files'].items()):
            raise ValueError('Missing or modified solve evidence')
    curves={}
    for c in CASES:
        raw=np.loadtxt(out/f"{c['id']}.csv",delimiter=',')
        if raw.shape!=(33,7) or not np.isfinite(raw).all() or not np.allclose(raw[:,0],p['frequencies_hz'],rtol=1e-12):
            raise ValueError('Invalid reference frequency grid')
        curves[c['id']]=dict(z=raw[:,1]+1j*raw[:,2],pressure=raw[:,5]+1j*raw[:,6])
    refinement=[]
    for a,b in zip(CASES,CASES[1:]):
        row=dict(coarse=a['id'],fine=b['id'])
        row.update({key:change(curves[a['id']][key],curves[b['id']][key]) for key in ('z','pressure')})
        refinement.append(row)
    radial=change(raw[:,3]+1j*raw[:,4],raw[:,5]+1j*raw[:,6])
    frame=pd.read_csv(out/'production.csv')
    if len(frame)!=33 or not np.isfinite(frame.select_dtypes(include='number')).all().all() or not np.allclose(frame.frequency,p['frequencies_hz'],rtol=1e-12):
        raise ValueError('Invalid production grid')
    for key,value in dict(bc_mode='velocity',phasor_convention='exp(+iwt)_rms',radiation_model='flanged_piston',loss_model='lossless',element_degree=1).items():
        if not (frame[key]==value).all():raise ValueError(f'Unexpected {key}')
    if (frame.relative_residual<0).any() or (frame.relative_residual>1e-8).any() or (frame.converged_reason<=0).any():
        raise ValueError('Unhealthy production solve')
    area=np.pi*p['geometry']['throat_radius']**2
    # Equal nominal volume velocity corrects polygonal inlet-area error.
    flow=(frame.mouth_u_real+1j*frame.mouth_u_imag).to_numpy()*area/frame.mesh_inlet_area_m2.to_numpy()
    pressure=baffled_piston_on_axis(frame.frequency.to_numpy(),flow,np.pi*p['geometry']['mouth_radius']**2,p['observer_distance_m'])
    reference=curves[CASES[-1]['id']]
    model=dict(pressure=change(pressure,reference['pressure']),z=change((frame.z_real+1j*frame.z_imag).to_numpy(),reference['z']))
    # Mode 32->64 and section 500->1000 are the final independent refinements.
    resolution_rows=[refinement[1],refinement[3]]
    converged=all(row[k]['level_db']<=LIMITS['reference_level_db'] and row[k]['phase_deg']<=LIMITS['reference_phase_deg'] for row in resolution_rows for k in ('z','pressure')) and radial['level_db']<=LIMITS['reference_level_db'] and radial['phase_deg']<=LIMITS['reference_phase_deg']
    passed=converged and model['pressure']['level_db']<=LIMITS['model_level_db'] and model['z']['level_db']<=LIMITS['model_impedance_db'] and all(v['phase_deg']<=LIMITS['model_phase_deg'] for v in model.values())
    clean_revision();verify(out)
    if (out/'comparison.json').exists():raise FileExistsError('Comparison already exists')
    result=dict(scope=p['scope'],limits=LIMITS,reference_refinement=refinement,radial_refinement=radial,
        model_changes=model,reference_converged=bool(converged),passed=bool(passed),physical_validation_status='experimental_prediction')
    write(out/'comparison.json',result);print(json.dumps(result,indent=2))
    if not passed:raise RuntimeError('Frozen comparison gates failed; preserve and investigate')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','reference','production','compare'])
    parser.add_argument('output',type=Path)
    parser.add_argument('--checkout',type=Path)
    parser.add_argument('--image',default='horn-octave-reference:local')
    args=parser.parse_args();out=args.output.resolve()
    if args.stage=='reference':
        if args.checkout is None:parser.error('reference requires --checkout')
        reference(out,args.checkout.resolve(),args.image)
    else:{'prepare':prepare,'production':production,'compare':compare}[args.stage](out)


if __name__=='__main__':main()
