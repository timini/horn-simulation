#!/usr/bin/env python3
"""Exhaustive numerical search audits: the analytical shortlist against FEM.

``finite_grid`` (the default, protocol v2) uses 20 small horns and the curated
manufacturer motors over four bands. Two pipeline cases run the analytical
stage exactly as auto mode does (driver prescreen, geometry grid, LEM
shortlist) and then solve every grid geometry with FEM:

* ``wide_320_5000``: the 6NMB420 example over 320-5000 Hz with no throat
  radius given, 100-300 mm long;
* ``size_capped_800_1600``: the curated motors over 800-1600 Hz with only
  ``max_length = 0.1`` given.

This tests acceleration against FEM on a finite grid, not physical accuracy or
commercial driver suitability. No benchmark results are used to fit the models.
Run inside horn-solver with the repository mounted at /workspace. The
analytical stage (``analytical_case``) needs no solver and is unit tested.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import hashlib
from horn_core.candidates import CandidateGeometry
from horn_drivers.loader import load_drivers
from horn_analysis.lem_prescreen import lem_prescreen_candidates
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec
from horn_analysis.search import screening_audit


ROOT = Path(__file__).resolve().parents[1]

# Pipeline cases: auto-mode parameters, frozen before any FEM result.
PIPELINE_CASES = {
    'wide_320_5000': dict(drivers='examples/6nmb420-320-5000/driver-database.json', f_low=320, f_high=5000,
                          min_length=.1, max_length=.3, num_mouth_radii=3, num_lengths=3, top_n=8,
                          radiation_model='modal_baffled', frequencies=161),
    'size_capped_800_1600': dict(drivers='data/drivers-curated', f_low=800, f_high=1600,
                                 max_length=.1, num_mouth_radii=3, num_lengths=3, top_n=10,
                                 radiation_model='modal_baffled', frequencies=100),
}


def analytical_case(name, num_mouth_radii=None, num_lengths=None, radiation_model=None, front_chamber_volume_m3=None):
    """Run auto mode's analytical stage for a pipeline case.

    ``front_chamber_volume_m3`` supplies a chamber volume for every driver in
    place of the default (cone estimate); 0 models a phase-plugged adapter.
    Returns a dict with the target, the prescreen result, the generated
    candidates, the design summary, the LEM screen and its shortlist rows.
    """
    from dataclasses import replace
    from horn_analysis.prescreen import PrescreenConfig, prescreen_drivers
    from horn_core.geometry_designer import generate_auto_candidates
    case = PIPELINE_CASES[name]
    drivers = load_drivers(str(ROOT/case['drivers']))
    if front_chamber_volume_m3 is not None:
        drivers = [replace(d, front_chamber_volume_m3=front_chamber_volume_m3, front_chamber_area_m2=None) for d in drivers]
    target = TargetSpec(case['f_low'], case['f_high'])
    screen_drivers = prescreen_drivers(drivers, PrescreenConfig(case['f_low'], case['f_high']))
    candidates, design = generate_auto_candidates(
        case['f_low'], case['f_high'], screen_drivers.throat_radii_m,
        num_mouth_radii=num_mouth_radii or case['num_mouth_radii'], num_lengths=num_lengths or case['num_lengths'],
        min_length=case.get('min_length'), max_length=case.get('max_length'))
    screen = lem_prescreen_candidates(candidates, screen_drivers.drivers, case['f_low'], case['f_high'],
                                      design.sim_freq_range, num_frequencies=case['frequencies'],
                                      top_n=case['top_n'], target=target,
                                      radiation_model=radiation_model or case['radiation_model'])
    best = {}
    for row in screen['rankings']:
        best.setdefault(row['candidate_id'], row)
    shortlist = [best[i] for i in screen['filtered_candidate_ids']]
    return {'name': name, 'target': target, 'prescreen': screen_drivers, 'candidates': candidates,
            'design': design, 'screen': screen, 'shortlist': shortlist, 'drivers': screen_drivers.drivers}


def _solve(c, out, band, frequencies, radiation_model, mesh_size=.008):
    from horn_geometry.generator import create_horn
    from horn_solver.solver import run_simulation_from_step
    import numpy as np
    import pandas as pd
    step=out/f'{c.candidate_id}.step';csv=out/f'{c.candidate_id}_results.csv'
    create_horn(c.profile,c.throat_radius,c.mouth_radius,c.length,step)
    run_simulation_from_step(str(step),band,frequencies,{'length':c.length},str(csv),band[1],mesh_size=mesh_size,radiation_model=radiation_model)
    # Require a complete healthy sweep before admitting any ranking rows.
    frame=pd.read_csv(csv)
    assert len(frame)==frequencies and np.allclose(frame.frequency,np.geomspace(*band,frequencies),rtol=1e-12,atol=0)
    assert np.isfinite(frame.select_dtypes(include='number')).all().all()
    assert (frame.radiation_model==radiation_model).all() and (frame.relative_residual.between(0,1e-8)).all()
    assert (frame.converged_reason>0).all() and (frame.z_real>=0).all()
    assert (frame.input_acoustic_power_w>=0).all() and (frame.mouth_acoustic_power_w>=0).all()
    balance=np.abs(frame.input_acoustic_power_w-frame.mouth_acoustic_power_w)/np.maximum(np.abs(frame.input_acoustic_power_w),1e-30)
    assert balance.max()<=1e-7
    return csv


def run_pipeline_case(name, out):
    """Analytical shortlist for one pipeline case, audited against FEM on every grid geometry."""
    case = analytical_case(name)
    out = out/name; out.mkdir(parents=True, exist_ok=True)
    band = tuple(case['design'].sim_freq_range)
    spec = PIPELINE_CASES[name]
    exhaustive = []
    for c in case['candidates']:
        csv = _solve(c, out, band, spec['frequencies'], spec['radiation_model'])
        for row in rank_horn_drivers(str(csv), c.candidate_id, c.throat_radius, case['drivers'], case['target'],
                                     top_n=len(case['drivers']), length=c.length):
            exhaustive.append({**row, 'profile': c.profile})
    audit = screening_audit(exhaustive, case['screen']['filtered_candidate_ids'])
    audit['feasible_pair_count'] = sum(r['model_feasible'] for r in exhaustive)
    return {'case': name, 'parameters': spec, 'target': asdict(case['target']),
            'size_compromise': case['design'].size_compromise,
            'shortlist': case['screen']['filtered_candidate_ids'], 'audit': audit,
            'exhaustive': exhaustive, 'passed': audit['passed']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True)
    p.add_argument('--radiation-model',choices=['flanged_piston','modal_baffled'],default='flanged_piston')
    p.add_argument('--cases', default='finite_grid',
                   help="Comma-separated: finite_grid, " + ", ".join(PIPELINE_CASES))
    a=p.parse_args();out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    names=[n.strip() for n in a.cases.split(',') if n.strip()]
    unknown=set(names)-{'finite_grid',*PIPELINE_CASES}
    if unknown: p.error(f"Unknown cases: {sorted(unknown)}")
    passed=True
    if 'finite_grid' in names:
        passed=finite_grid(a, out)==0 and passed
    for name in names:
        if name in PIPELINE_CASES:
            result=run_pipeline_case(name, out)
            (out/f'{name}_benchmark.json').write_text(json.dumps(result,indent=2,default=str))
            print(json.dumps(result['audit'],indent=2))
            passed=passed and result['passed']
    return 0 if passed else 1


def finite_grid(a, out):
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
        csv=_solve(c,out,band,100,a.radiation_model)
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
