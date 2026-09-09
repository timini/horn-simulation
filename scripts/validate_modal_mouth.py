#!/usr/bin/env python3
"""Qualify the experimental nonlocal aperture against the archived MMM reference.

The input reference is immutable, independently converged, and uses the same
ideal axisymmetric infinite-baffle problem. No physical qualification is implied.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CASES={'m8_h4':(.004,8),'m16_h4':(.004,16),'m16_h3':(.003,16)}
LIMITS={'reference_level_db':.5,'reference_phase_deg':5.,'refinement_level_db':.05,'refinement_phase_deg':1.}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def source():
    paths=list((ROOT/'packages').glob('*/src/**/*.py'))+[Path(__file__)]
    return {str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}
def clean():
    if subprocess.check_output(['git','status','--porcelain','--untracked-files=all'],cwd=ROOT,text=True):
        raise ValueError('Clean checkout required')
    tracked=set(subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines())
    if set(source())-tracked:raise ValueError('Untracked executable source')
    return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
def verify(out):
    p=json.loads((out/'protocol.json').read_text())
    if p['source']!=source() or p['cases']!={k:list(v) for k,v in CASES.items()} or p['limits']!=LIMITS:
        raise ValueError('Source/protocol changed')
    if any(sha(out/k)!=v for k,v in p['inputs'].items()):raise ValueError('Input changed')
    return p

def prepare(out):
    revision=clean();directory=ROOT/'data/validation'
    manifest=json.loads((directory/'mmm_radiation_manifest.json').read_text())
    for name,digest in manifest['files'].items():
        if sha(directory/name)!=digest:raise ValueError('Reference archive identity changed')
    result=json.loads((directory/'mmm_radiation_reference.json').read_text())
    if not result['reference_converged']:raise ValueError('Reference has not converged')
    out.mkdir(parents=True,exist_ok=False)
    with tarfile.open(directory/'mmm_radiation_artifacts.tar.gz') as archive:
        for name in ('horn.step','m64_s1000.csv','protocol.json'):
            matches=[m for m in archive.getmembers() if m.isfile() and Path(m.name).name==name]
            if len(matches)!=1:raise ValueError('Ambiguous reference input')
            (out/('reference-'+name)).write_bytes(archive.extractfile(matches[0]).read())
    shutil.copyfile(directory/'mmm_radiation_manifest.json',out/'reference-manifest.json')
    shutil.copyfile(directory/'mmm_radiation_reference.json',out/'reference-comparison.json')
    p=json.loads((out/'reference-protocol.json').read_text())
    if p['geometry']!=dict(throat_radius=.032498,mouth_radius=.08870898640584517,length=.11790625,profile='hyperbolic'):
        raise ValueError('Unexpected reference geometry')
    write(out/'protocol.json',dict(source_revision=revision,source=source(),cases=CASES,limits=LIMITS,
        geometry=p['geometry'],frequencies_hz=p['frequencies_hz'],air=p['air'],observer_distance_m=p['observer_distance_m'],
        inputs={f.name:sha(f) for f in out.iterdir()},scope='Numerical aperture qualification only; infinite baffle, axisymmetric ideal horn'))

def solve(out):
    sys.path[:0]=[str(p/'src') for p in (ROOT/'packages').glob('horn-*')]
    from dolfinx import fem
    from petsc4py import PETSc
    import ufl
    from horn_solver.solver import create_mesh_from_step
    from horn_solver.modal_boundary import aperture_projection,solve_modal_aperture
    from horn_core.modal_radiation import modal_baffled_on_axis
    p=verify(out);g=p['geometry'];radius=g['mouth_radius'];nominal=np.pi*g['throat_radius']**2
    for name,(size,modes) in CASES.items():
        destination=out/f'{name}.json'
        if destination.exists():raise FileExistsError(destination)
        mesh,tags=create_mesh_from_step(str(out/'reference-horn.step'),size,g['length'])
        V=fem.functionspace(mesh,('Lagrange',1));ds=ufl.Measure('ds',domain=mesh,subdomain_data=tags)
        trial=ufl.TrialFunction(V);test=ufl.TestFunction(V);one=fem.Constant(mesh,PETSc.ScalarType(1))
        inlet=fem.assemble_scalar(fem.form(one*ds(2))).real
        B=aperture_projection(V,ds,3,radius,modes);rows=[]
        rho=p['air']['rho'];c=p['air']['c']
        for f in p['frequencies_hz']:
            k=2*np.pi*f/c
            a=(ufl.inner(ufl.grad(trial),ufl.grad(test))-k*k*ufl.inner(trial,test))*ufl.dx
            L=ufl.inner(fem.Constant(mesh,PETSc.ScalarType(1j*k*c*rho)),test)*ds(2)
            pressure,velocity,health=solve_modal_aperture(V,a,L,[],B,radius,k,rho=rho,c=c)
            pin=fem.assemble_scalar(fem.form(pressure*ds(2)))/inlet
            impedance=pin*nominal/inlet
            field=modal_baffled_on_axis(np.array([f]),velocity[None,:]*nominal/inlet,radius,p['observer_distance_m'],rho=rho,c=c)[0]
            power=float(np.real(pin*inlet))
            balance=abs(power-health['radiated_power_w'])/max(abs(power),1e-30)
            if power<=0 or health['radiated_power_w']<=0 or balance>1e-8:raise ValueError('Power balance failed')
            rows.append(dict(frequency=f,z_real=impedance.real,z_imag=impedance.imag,p_real=field.real,p_imag=field.imag,
                mesh_cells=mesh.topology.index_map(mesh.topology.dim).size_global,power_balance=balance,**health))
        write(destination,rows);print('Completed modal case',name,flush=True)
    verify(out)
    write(out/'solve-evidence.json',dict(protocol_sha256=sha(out/'protocol.json'),files={f'{name}.json':sha(out/f'{name}.json') for name in CASES}))

def change(a,b):
    if not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(abs(a)==0) or np.any(abs(b)==0):raise ValueError('Invalid comparison')
    return dict(level_db=float(np.max(abs(20*np.log10(abs(a)/abs(b))))),phase_deg=float(np.max(abs(np.angle(a*np.conj(b),deg=True)))))

def compare(out):
    clean();p=verify(out);e=json.loads((out/'solve-evidence.json').read_text())
    if e['protocol_sha256']!=sha(out/'protocol.json') or set(e['files'])!={f'{name}.json' for name in CASES} or any(sha(out/k)!=v for k,v in e['files'].items()):raise ValueError('Solve evidence changed')
    reference=np.loadtxt(out/'reference-m64_s1000.csv',delimiter=',')
    if reference.shape!=(33,7) or not np.allclose(reference[:,0],p['frequencies_hz'],rtol=1e-12):raise ValueError('Reference grid changed')
    curves={'reference':dict(z=reference[:,1]+1j*reference[:,2],pressure=reference[:,5]+1j*reference[:,6])};health={}
    for name in CASES:
        rows=json.loads((out/f'{name}.json').read_text())
        if len(rows)!=33 or not np.allclose([r['frequency'] for r in rows],p['frequencies_hz'],rtol=1e-12):raise ValueError('Incomplete solve grid')
        if not all(np.isfinite(list(r.values())).all() for r in rows):raise ValueError('Nonfinite solve')
        for key in ('power_balance','relative_residual','interface_relative_error'):
            if any(not 0<=r[key]<=1e-8 for r in rows):raise ValueError('Failed solve health')
        if any(r['radiated_power_w']<=0 for r in rows):raise ValueError('Nonpassive solve')
        curves[name]=dict(z=np.array([r['z_real']+1j*r['z_imag'] for r in rows]),pressure=np.array([r['p_real']+1j*r['p_imag'] for r in rows]))
        health[name]={key:max(r[key] for r in rows) for key in ('power_balance','relative_residual','interface_relative_error','mesh_cells')}
    if health['m16_h3']['mesh_cells']<=health['m16_h4']['mesh_cells']:raise ValueError('Mesh did not refine')
    comparisons=[]
    for kind,a,b in [('refinement','m8_h4','m16_h4'),('refinement','m16_h4','m16_h3'),('reference','m16_h3','reference')]:
        metrics={key:change(curves[a][key],curves[b][key]) for key in ('z','pressure')}
        passed=all(value<=LIMITS[f'{kind}_{metric}'] for item in metrics.values() for metric,value in item.items())
        comparisons.append(dict(kind=kind,coarse=a,fine=b,metrics=metrics,passed=bool(passed)))
    result=dict(scope=p['scope'],protocol_sha256=sha(out/'protocol.json'),solve_evidence_sha256=sha(out/'solve-evidence.json'),health=health,comparisons=comparisons,passed=all(r['passed'] for r in comparisons),physical_validation_status='experimental_prediction')
    if (out/'comparison.json').exists():raise FileExistsError('Comparison already exists')
    write(out/'comparison.json',result);print(json.dumps(result,indent=2))
    if not result['passed']:raise RuntimeError('Frozen qualification gates failed')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['prepare','solve','compare']);parser.add_argument('output',type=Path)
    args=parser.parse_args();globals()[args.stage](args.output.resolve())
