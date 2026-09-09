#!/usr/bin/env python3
"""Audit retained finite-baffle pilot outputs without declaring an accuracy pass."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'packages/horn-core/src'))
from horn_core.acoustics import baffled_piston_on_axis


def audit(root):
    nominal_area=np.pi*.032498**2
    mouth_area=np.pi*.08870898640584517**2
    configurations=[('exterior-pilot-v3','boundary-lab','float32',2),
                    ('exterior-pilot-v3','fp64-q4','float64',4),
                    ('exterior-pilot-v3','fp64-q6','float64',6),
                    ('exterior-pilot-v4','fp64-q4','float64',4),
                    ('exterior-pilot-v5','fp64-q4','float64',4)]
    rows=[]
    for case,variant,precision,quadrature in configurations:
        directory=root/case/variant
        domains=json.loads((directory/'domains.json').read_text())['domains']
        boundary=next(d for d in domains if d['id']=='domain:bem-boundary')
        with np.load(directory/'domains.npz') as arrays:
            xyz=arrays[boundary['coordinates']['points_m']]
            triangles=arrays[boundary['topology']['triangles']]
            tags=arrays[boundary['topology']['source_physical_tag']]
        points=xyz[triangles]
        cross=np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0])
        areas=np.linalg.norm(cross,axis=1)/2
        driven=tags==2;source_area=areas[driven].sum()
        if source_area<=0 or np.any(areas<=0) or not np.all(cross[driven,2]>0):
            raise ValueError('Invalid driven surface geometry')
        edges=np.sort(np.vstack([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]]),axis=1)
        if not np.all(np.unique(edges,axis=0,return_counts=True)[1]==2) or not np.allclose(cross.sum(axis=0),0,atol=1e-10):
            raise ValueError('Surface is not a closed oriented solid')
        definition=json.loads((root/case/'pilot-definition.json').read_text())
        for index,frequency in enumerate((800.,1200.,1600.)):
            meta=json.loads((directory/f'frequencies/{index:06}.json').read_text())
            diagnostics=meta['diagnostics']
            if meta['freq_hz']!=frequency or diagnostics['precision']!=precision or diagnostics['regular_quadrature_order']!=quadrature or diagnostics['formulation']!='exterior_burton_miller_neumann':
                raise ValueError('Unexpected reference solve contract')
            quantities={q['quantity']:q for q in meta['quantities']}
            if quantities['radiation_impedance']['unit']!='N*s/m' or quantities['exterior_pressure']['unit']!='Pa':
                raise ValueError('Reference units changed')
            with np.load(directory/f'frequencies/{index:06}.npz') as arrays:
                values={name:arrays[q['key']] for name,q in quantities.items()}
            if any(not np.isfinite(a).all() for a in values.values()):
                raise ValueError('Nonfinite reference output')
            force_z=values['radiation_impedance'].item()
            boundary_pressure=values['bem_boundary_pressure'][0]
            derivative=values['bem_boundary_neumann'][0]
            integrated_force=(boundary_pressure[triangles[driven]].mean(axis=1)*areas[driven]).sum()
            force_error=abs(integrated_force-force_z)/abs(force_z)
            velocity_error=np.max(abs(derivative[driven]/(1j*2*np.pi*frequency*1.225)-1))
            if force_error>1e-5 or velocity_error>1e-5 or np.any(derivative[~driven]!=0):
                raise ValueError('Source amplitude or force reduction does not match')
            # Convert exp(-iwt), mechanical load, and polygonal source area.
            z=force_z.conjugate()*nominal_area/source_area**2
            pressure=values['exterior_pressure'].item().conjugate()*nominal_area/source_area
            frame=pd.read_csv(root/f'exterior-pilot-v3/production-{frequency:g}.csv')
            row=frame.iloc[0]
            if row.frequency!=frequency or row.bc_mode!='velocity' or row.phasor_convention!='exp(+iwt)_rms' or row.relative_residual>1e-8 or row.converged_reason<=0:
                raise ValueError('Production pilot contract changed')
            flow=(row.mouth_u_real+1j*row.mouth_u_imag)*nominal_area/row.mesh_inlet_area_m2
            predicted=baffled_piston_on_axis(np.array([frequency]),np.array([flow]),mouth_area,1.)[0]
            predicted_z=row.z_real+1j*row.z_imag
            rows.append(dict(case=case,variant=variant,frequency_hz=frequency,
                baffle_radius_m=definition['baffle_radius'],faces=len(triangles),source_area_m2=float(source_area),
                source_velocity_relative_error=float(velocity_error),force_reduction_relative_error=float(force_error),
                reference_pressure_pa=[float(pressure.real),float(pressure.imag)],reference_impedance_pa_s_per_m=[float(z.real),float(z.imag)],
                output_difference_db=float(20*np.log10(abs(pressure/predicted))),
                impedance_difference_db=float(20*np.log10(abs(z/predicted_z))),
                output_phase_difference_deg=float(np.angle(pressure/predicted,deg=True))))
    return dict(scope='Method pilot at three frequencies; finite mounting differs from infinite-baffle approximation',
        accuracy_pass=False,physical_validation_status='experimental_prediction',
        notes=['No exterior linear-system residual is exported by this pinned backend.',
               'Excitation and integrated-force closures are consistency checks, not a residual or physical validation.',
               'Mesh, precision, quadrature and baffle-size changes are retained; no continuous-band bound is claimed.'],rows=rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    result=audit(args.root)
    result['raw_files_sha256']={str(p.relative_to(args.root)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(args.root.glob('exterior-pilot-v[345]/**/*')) if p.is_file()}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
