#!/usr/bin/env python3
"""Exhaustive numerical search audit using 20 small horns and three manufacturer-sourced motors.

This tests acceleration against FEM on a finite grid, not physical accuracy or
commercial driver suitability. No benchmark results are used to fit the models.
Run inside horn-solver with the repository mounted at /workspace.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import hashlib
from horn_core.candidates import CandidateGeometry
from horn_drivers.loader import load_drivers
from horn_geometry.generator import create_horn
from horn_solver.solver import run_simulation_from_step
from horn_analysis.lem_prescreen import lem_prescreen_candidates
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec
from horn_analysis.search import screening_audit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True)
    p.add_argument('--radiation-model',choices=['flanged_piston','modal_baffled'],default='flanged_piston')
    a=p.parse_args();out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    candidates=[CandidateGeometry(f'grid_{i:02d}',profile,.0225,mouth,length)
                for i,(profile,mouth,length) in enumerate((p,m,l) for p in ('conical','exponential') for m in (.045,.07) for l in (.04,.055,.07,.085,.10))]
    drivers=load_drivers(str(Path(__file__).resolve().parents[1]/'data/drivers-curated'))
    # Frozen before this benchmark: no tuning of scores or gates from FEM results.
    bands=[(800,1600,'development'),(1000,1200,'held_out'),(1200,2000,'held_out'),(600,1000,'held_out')]
    targets=[TargetSpec(lo,hi) for lo,hi,_ in bands]
    band=(600/2**.5,2000*2**.5)
    started=time.perf_counter()
    screens=[lem_prescreen_candidates(candidates,drivers,t.f_low_hz,t.f_high_hz,band,num_frequencies=100,top_n=10,target=t,radiation_model=a.radiation_model) for t in targets]
    screening_seconds=time.perf_counter()-started
    protocol={'benchmark':'finite_grid_manufacturer_motors_v2', 'bands':bands,
              'geometry_count':len(candidates),'driver_count':len(drivers),
              'sim_band':band,'frequencies':100,'mesh_size':.008,
              'shortlist_budget':10,'score_regret_limit':.02,'top_ten_recall_limit':.9,
              'minimum_feasible_pairs_per_case':10,
              'physical_validation_passed':False,
              'radiation_model':a.radiation_model,
              'driver_scope':'Manufacturer parameters; interfaces and free-air mass substitution remain experimental'}
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2))
    rows=[[] for _ in targets]
    for c in candidates:
        step=out/f'{c.candidate_id}.step';csv=out/f'{c.candidate_id}_results.csv'
        create_horn(c.profile,c.throat_radius,c.mouth_radius,c.length,step)
        run_simulation_from_step(str(step),band,100,{'length':c.length},str(csv),band[1],mesh_size=.008,radiation_model=a.radiation_model)
        # Require a complete healthy sweep before admitting any ranking rows.
        import pandas as pd
        import numpy as np
        frame=pd.read_csv(csv)
        assert len(frame)==100 and np.allclose(frame.frequency,np.geomspace(*band,100),rtol=1e-12,atol=0)
        assert np.isfinite(frame.select_dtypes(include='number')).all().all()
        assert (frame.radiation_model==a.radiation_model).all() and (frame.relative_residual.between(0,1e-8)).all()
        assert (frame.converged_reason>0).all() and (frame.z_real>=0).all()
        assert (frame.input_acoustic_power_w>=0).all() and (frame.mouth_acoustic_power_w>=0).all()
        balance=np.abs(frame.input_acoustic_power_w-frame.mouth_acoustic_power_w)/np.maximum(np.abs(frame.input_acoustic_power_w),1e-30)
        assert balance.max()<=1e-7
        for i,target in enumerate(targets):
            rows[i].extend(rank_horn_drivers(str(csv),c.candidate_id,c.throat_radius,drivers,target,top_n=len(drivers)))
    cases=[]
    for target,(_,_,role),screen,exhaustive in zip(targets,bands,screens,rows):
        audit=screening_audit(exhaustive,screen['filtered_candidate_ids'])
        audit['feasible_pair_count']=sum(r['model_feasible'] for r in exhaustive)
        audit['sufficient_feasible_pairs']=audit['feasible_pair_count']>=10
        audit['passed']=audit['passed'] and audit['sufficient_feasible_pairs']
        cases.append({'target':asdict(target),'role':role,'screening':screen,'exhaustive':exhaustive,'audit':audit})
    result={**protocol,'candidates':[asdict(c) for c in candidates], 'drivers':[asdict(d) for d in drivers],
            'screening_seconds':screening_seconds,'total_seconds':time.perf_counter()-started,
            'cases':cases,'passed':all(c['audit']['passed'] for c in cases)}
    result['raw_sha256']={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(out.glob('grid_*')) if path.is_file()}
    (out/'search_benchmark.json').write_text(json.dumps(result,indent=2))
    print(json.dumps([c['audit'] for c in cases],indent=2))
    return 0 if result['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
