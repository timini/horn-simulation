"""Frozen, unfitted duct-loss validation against acquired independent references.

Run inside horn-solver with repository sources available. All model predictions
are saved before measurement values are read. Brass open cylinders are development
cases; other materials and cones are held out. This validates input impedance,
not a loudspeaker motor, exterior sound field, or complete horn recommendation.
"""
import argparse
from collections import defaultdict
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from horn_core.duct import AirProperties
from horn_core.webster import compute_horn_transfer_tmm

AIR=AirProperties(c=346.28592,rho=1.184490,gamma=1.402,viscosity=1.83183e-5,
                  conductivity=.02613336,heat_capacity=1004.16)
CASES={'Brass_O':(.007,.007,.002,'finite_flange'),
       'Wood_O':(.007,.007,.007,'finite_flange'),
       '3D_O':(.007,.007,.007,'finite_flange'),
       'Brass_C':(.007,.007,0.,'closed'),
       'Wood_C':(.007,.007,0.,'closed'),
       '3D_C':(.007,.007,0.,'closed'),
       'Cone_O':(.005,.0113,.0027,'finite_flange'),
       'Cone_C':(.005,.0113,0.,'closed')}


def predict(frequency, case, loss='boundary_layer', segments=400):
    throat,mouth,width,radiation=CASES[case]
    t=compute_horn_transfer_tmm(frequency,lambda z: throat+(mouth-throat)*z/.18,
                              .18,throat,mouth,n_segments=segments,air=AIR,
                              loss_model=loss,radiation_model=radiation,flange_width=width)
    return (t['z_real']+1j*t['z_imag'])/(AIR.rho*AIR.c)


def metrics(predicted, measured):
    error=np.abs(20*np.log10(np.maximum(np.abs(predicted),1e-15)/np.maximum(np.abs(measured),1e-15)))
    return {'median_db':float(np.median(error)), 'p95_db':float(np.percentile(error,95)),
            'max_db':float(error.max()),
            'p95_scaled_complex_error':float(np.percentile(np.abs(predicted-measured)/(1+np.abs(measured)),95)),
            'magnitude_agreement':bool(np.median(error)<=2 and np.percentile(error,95)<=4)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-dir',required=True);p.add_argument('--output-dir',required=True)
    p.add_argument('--with-fem',action='store_true')
    a=p.parse_args();root=Path(a.reference_dir);out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((root/'manifest.json').read_text())
    frequency=np.geomspace(110,3900,20001)
    predicted={}
    for case in CASES:
        z=predict(frequency,case);predicted[case]=z
        pd.DataFrame({'frequency':frequency,'z_real_normalized':z.real,'z_imag_normalized':z.imag}).to_csv(out/f'{case}-prediction.csv',index=False)
    protocol={'air':asdict(AIR),'cases':CASES,'frequency_range_hz':[110,3900],
              'loss_model':'boundary_layer','fit_parameters':[],
              'magnitude_limits_db':{'median':2,'p95':4},
              'prediction_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in out.glob('*-prediction.csv')},
              'source_archive_sha256':manifest['reference']['sha256'],
              'source_url':manifest['reference']['source_url'],
              'code_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),Path('packages/horn-core/src/horn_core/duct.py'),Path('packages/horn-core/src/horn_core/webster.py'),Path('packages/horn-solver/src/horn_solver/solver.py')]}}
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    fem={}
    if a.with_fem:
        import gmsh
        from horn_solver.solver import run_simulation_from_step
        for case in ['Brass_O','Wood_O','Cone_O']:
            throat,mouth,width,radiation=CASES[case];step=out/f'{case}.step'
            gmsh.initialize()
            try:
                gmsh.model.add(case)
                if throat==mouth:gmsh.model.occ.addCylinder(0,0,0,0,0,.18,throat)
                else:gmsh.model.occ.addCone(0,0,0,0,0,.18,throat,mouth)
                gmsh.model.occ.synchronize();gmsh.write(str(step))
            finally:gmsh.finalize()
            csv=out/f'{case}-fem.csv'
            run_simulation_from_step(str(step),(110,3900),121,{'length':.18},str(csv),3900,
                                     mesh_size=.004,element_degree=2,loss_model='boundary_layer',
                                     minimum_wall_scale=throat,radiation_model=radiation,
                                     flange_width=width,air=AIR)
            frame=pd.read_csv(csv);fem[case]=frame
    rows=[];measurements=defaultdict(list)
    for curve in manifest['curves']:
        if curve['reference_kind']!='measurement' or curve.get('duplicate_of'):continue
        case=curve['configuration'];path=root/curve['csv']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=curve['csv_sha256']:
            raise ValueError(f'Changed reference: {path}')
        frame=pd.read_csv(path)
        if frame.frequency.iloc[0]>110 or frame.frequency.iloc[-1]<3900:
            rows.append({'source_member':curve['source_member'],'case':case,'status':'outside_fixed_band'});continue
        measured=np.interp(frequency,frame.frequency,frame.z_real_normalized)+1j*np.interp(frequency,frame.frequency,frame.z_imag_normalized)
        result={'source_member':curve['source_member'],'case':case,'operator':curve['operator'],
                'role':'development' if case=='Brass_O' else 'held_out',**metrics(predicted[case],measured)}
        if case in fem:
            ff=fem[case].frequency.to_numpy();fz=(fem[case].z_real.to_numpy()+1j*fem[case].z_imag.to_numpy())/(AIR.rho*AIR.c)
            mz=np.interp(ff,frame.frequency,frame.z_real_normalized)+1j*np.interp(ff,frame.frequency,frame.z_imag_normalized)
            result['production_fem_at_121_frequencies']=metrics(fz,mz)
        rows.append(result);measurements[case].append(measured)
    summary={case:{'curves':len([r for r in rows if r['case']==case]),
                   'magnitude_passes':sum(r.get('magnitude_agreement',False) for r in rows if r['case']==case),
                   'median_p95_db':float(np.median([r['p95_db'] for r in rows if r['case']==case and 'p95_db' in r]))}
             for case in CASES}
    numerical={}
    # A separate operator's supplied numerical data; never call these measured.
    mapping={'Cylinder_closed':'Brass_C','Cylinder_finite_flanged_width2mm':'Brass_O','Cylinder_finite_flanged_width7mm':'Wood_O','Cone_closed':'Cone_C'}
    for c in manifest['curves']:
        if c['reference_kind']!='independent_simulation' or c['configuration'] not in mapping:continue
        if '_Operator-C_1DFEM_' not in c['source_member'] or ('Cone' in c['configuration'] and 'PlaneWave' not in c['source_member']):continue
        path=root/c['csv']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=c['csv_sha256']:raise ValueError('Changed numerical reference')
        d=pd.read_csv(path);case=mapping[c['configuration']];area=np.pi*CASES[case][0]**2
        z=(np.interp(frequency,d.frequency,d.z_real_pa_s_per_m3)+1j*np.interp(frequency,d.frequency,d.z_imag_pa_s_per_m3))*area/(AIR.rho*AIR.c)
        numerical[c['source_member']]=metrics(predicted[case],z)
    result={'protocol':protocol,'cases':summary,'measurements':rows,'independent_numerical_comparisons':numerical,
            'physical_horn_workflow_validated':False,
            'scope':'Nominal rigid pipe input impedance only; magnitude gates do not validate motor, output SPL or all specimen uncertainties.'}
    (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for ax,case in zip(axes,['Brass_O','Wood_O','Cone_O']):
        measured=np.median(np.abs(measurements[case]),axis=0)
        ax.semilogx(frequency,20*np.log10(measured),label='Median measurement')
        ax.semilogx(frequency,20*np.log10(np.abs(predicted[case])),'--',label='Unfitted loss/flange model')
        if case in fem:
            d=fem[case];z=(d.z_real+1j*d.z_imag)/(AIR.rho*AIR.c)
            ax.semilogx(d.frequency,20*np.log10(np.abs(z)),'.',ms=2,label='Production FEM')
        ax.set_title(case+' (development)' if case=='Brass_O' else case+' (held out)')
        ax.set_xlabel('Frequency (Hz)');ax.grid(alpha=.2);ax.set_xlim(110,3900)
    axes[0].set_ylabel('20 log10 |Z/(ρc/S)| (dB)');axes[0].legend(fontsize=7)
    fig.suptitle('Published pipe measurements versus predictions with losses and finite flanges')
    fig.tight_layout();fig.savefig(out/'validation.png',dpi=160);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
