"""Refine a FEM shortlist inside the solver container with a fixed budget."""
import argparse
import json
from dataclasses import asdict
from pathlib import Path
import shutil
import numpy as np
from horn_core.candidates import CandidateGeometry
from horn_analysis.search import refine, geometry_key, search_objective
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec
from horn_core.frequency_grid import DEFAULT_MIN_POINTS_PER_OCTAVE, band_grids
from horn_drivers.loader import load_drivers

# Objective for a geometry that produced no metrics (CAD, profile or radiation
# domain rejection): below any measured infeasibility, but finite.
REJECTED_GEOMETRY_OBJECTIVE = -1e6
WIDENABLE_SOURCES = {"heuristic", "driver_prescreen", "size_cap_fraction"}


def select_seed(data, max_ripple_db=6.0):
    """Refinement seed: the best feasible FEM row, else the least-infeasible one.

    Returns (row, source) with source "best_feasible" or "least_infeasible",
    or (None, None) when no FEM row carries metrics (every geometry was
    rejected before solving).
    """
    if data.get('results'):
        return data['results'][0], 'best_feasible'
    measured = [row for row in data.get('rejected', [])
                if row.get('passband_ripple_db') is not None and row.get('simulation_eligible', True)
                and all(row.get(k) is not None for k in ('horn_label', 'profile', 'throat_radius', 'mouth_radius', 'length'))]
    if not measured:
        return None, None
    best = max(measured, key=lambda row: (search_objective(row, max_ripple_db), row['horn_label']))
    return best, 'least_infeasible'


def widening_policy(design, screen):
    """Bound sides refinement may widen, and the hard limits it never crosses.

    Only heuristic sides (and prescreen throat sampling) widen; user-set sides
    stay fixed. The throat never exceeds the prescreen ka cap at f_high.
    """
    sources = design.get('bound_sources') or {}
    widenable = {field: tuple(side for side in ('low', 'high') if sides.get(side) in WIDENABLE_SOURCES)
                 for field, sides in sources.items()}
    ka_max = design.get('throat_ka_max', screen.get('ka_max'))
    limits = {}
    if ka_max is not None:
        limits['throat_radius'] = (0.0, 343. * ka_max / (2 * np.pi * design['target_f_high']))
    return widenable, limits


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ranked-json',required=True);p.add_argument('--design-json',required=True)
    p.add_argument('--prescreen-json',required=True);p.add_argument('--drivers-db',required=True)
    p.add_argument('--solver-csvs',nargs='+',required=True);p.add_argument('--output-dir',required=True)
    p.add_argument('--top-n',type=int,default=10)
    p.add_argument('--budget',type=int,default=6);p.add_argument('--num-frequencies',type=int,default=100)
    p.add_argument('--min-points-per-octave',type=float,default=DEFAULT_MIN_POINTS_PER_OCTAVE)
    p.add_argument('--min-output',type=float,default=None);p.add_argument('--score-weights',default=None)
    p.add_argument('--num-bands',type=int,default=8);p.add_argument('--num-sections',type=int,default=20)
    p.add_argument('--mesh-size',type=float,default=.01);p.add_argument('--radiation-model',default='flanged_piston')
    p.add_argument('--voltage',type=float,default=2.83);p.add_argument('--distance',type=float,default=1.)
    p.add_argument('--max-ripple',type=float,default=6.);p.add_argument('--max-compression',type=float,default=10.)
    p.add_argument('--test-power',type=float,default=None)
    p.add_argument("--loss-model",choices=["lossless", "boundary_layer"],default="lossless")
    p.add_argument("--flange-width",type=float,default=0.)
    p.add_argument("--minimum-wall-scale",type=float)
    p.add_argument("--element-degree",type=int,default=1)
    p.add_argument("--max-mesh-nodes",type=int,default=None)
    a=p.parse_args()
    if a.num_bands < 1 or a.num_frequencies < 2 or a.budget < 0 or a.top_n < 1: p.error("Invalid simulation budget/grid")
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    data=json.loads(Path(a.ranked_json).read_text());design=json.loads(Path(a.design_json).read_text())
    for path in a.solver_csvs: shutil.copyfile(path,out/Path(path).name)
    seed_row, seed_source = select_seed(data, a.max_ripple)
    if a.budget==0 or seed_row is None:
        data['refinement']={'new_evaluations':0,'budget':a.budget,'reason':'disabled' if a.budget==0 else 'no_seed_with_metrics'}
    else:
        from horn_geometry.generator import create_horn, LoftAccuracyError
        from horn_core.parameters import HornParameters, FlareProfile
        from horn_solver.solver import run_simulation_from_step
        from horn_solver.meshing import rejection_reason
        screen=json.loads(Path(a.prescreen_json).read_text())
        drivers=[d for d in load_drivers(a.drivers_db) if d.driver_id in screen['drivers']]
        target=TargetSpec(design['target_f_low'],design['target_f_high'],voltage_rms=a.voltage,observation_distance_m=a.distance,max_ripple_db=a.max_ripple,max_compression_ratio=a.max_compression,min_output_db=a.min_output,score_weights=a.score_weights,test_power_w=a.test_power)
        r=seed_row
        seed=CandidateGeometry(candidate_id=r['horn_label'],profile=r['profile'],throat_radius=r['throat_radius'],mouth_radius=r['mouth_radius'],length=r['length'])
        # The search objective is the composite score for feasible designs and
        # minus the infeasibility for rejected ones, so a search that starts
        # from the least-infeasible FEM candidate still moves towards feasible.
        objective=lambda row: search_objective(row, a.max_ripple)
        cache={geometry_key(seed):max((objective(row) for row in data['results']+data.get('rejected',[])
                                       if row.get('horn_label')==r['horn_label'] and row.get('passband_ripple_db') is not None),
                                      default=objective(r))}
        all_rows=list(data['results']); rejected=list(data.get('rejected',[]))
        from horn_analysis.evaluation import radiation_domain_rejection
        low, high = design['sim_freq_range']
        grids = band_grids(low, high, a.num_bands, a.num_frequencies, a.min_points_per_octave)
        points = len(grids[0][2])
        simulation_frequencies = np.concatenate([f for _, _, f in grids])
        domain_rejected = []
        geometry_rejected = []
        def evaluate(c):
            key=geometry_key(c)
            if key in cache:return cache[key]
            rejection = radiation_domain_rejection(simulation_frequencies, c.mouth_radius, a.radiation_model, a.flange_width)
            if rejection:
                domain_rejected.append(c.candidate_id)
                rejected.extend({"horn_label": c.candidate_id, "driver_id": drv.driver_id,
                                 "profile": c.profile, "throat_radius": c.throat_radius,
                                 "mouth_radius": c.mouth_radius, "length": c.length,
                                 **rejection} for drv in drivers)
                cache[key] = REJECTED_GEOMETRY_OBJECTIVE
                return REJECTED_GEOMETRY_OBJECTIVE
            def reject(reason, detail):
                # Drop this geometry for every driver with one stated reason.
                geometry_rejected.append(c.candidate_id)
                rejected.extend({"horn_label": c.candidate_id, "driver_id": drv.driver_id,
                    "profile": c.profile, "throat_radius": c.throat_radius,
                    "mouth_radius": c.mouth_radius, "length": c.length,
                    "model_feasible": False, "simulation_eligible": False,
                    "composite_score": 0., "extension_octaves": None,
                    "passband_ripple_db": None, "avg_sensitivity_db": None,
                    "rejection_reasons": [reason],
                    "rejection_detail": detail, "evidence_gaps": []} for drv in drivers)
                cache[key] = REJECTED_GEOMETRY_OBJECTIVE
                return REJECTED_GEOMETRY_OBJECTIVE
            from horn_core.profiles import profile_geometry_error
            profile_error = profile_geometry_error(c.profile, c.throat_radius, c.mouth_radius, c.length)
            if profile_error:
                return reject("profile_geometry_infeasible", profile_error)
            step=out/f'{c.candidate_id}.step'
            try:
                create_horn(c.profile,c.throat_radius,c.mouth_radius,c.length,step,num_sections=a.num_sections)
            except LoftAccuracyError as error:
                return reject("loft_accuracy_unmet", str(error))
            csv=out/f'{c.candidate_id}_results.csv'
            from horn_analysis.merge import BandBoundaryMismatch, merge_bands
            band_dir=out/f'{c.candidate_id}_bands';band_dir.mkdir()
            paths=[]
            try:
                for index,(band_low,band_high,_) in enumerate(grids):
                    path=band_dir/f'results_{index}.csv'
                    run_simulation_from_step(str(step),(band_low,band_high),points,{'length':c.length},str(path),band_high,mesh_size=a.mesh_size,radiation_model=a.radiation_model,loss_model=a.loss_model,flange_width=a.flange_width,minimum_wall_scale=a.minimum_wall_scale,element_degree=a.element_degree,max_mesh_nodes=a.max_mesh_nodes)
                    paths.append(path)
            except Exception as error:
                # Mesh-size, throat-area and CAD rejections drop this candidate;
                # anything else is a real fault and stops the run.
                reason = rejection_reason(error)
                if reason is None:
                    raise
                return reject(reason, str(error))
            try:
                merge_bands(paths,num_bands=a.num_bands,min_freq=low,max_freq=high,points_per_band=points,output=csv)
            except BandBoundaryMismatch as error:
                # Bands that do not join mean an unresolved mesh: reject this candidate.
                return reject(error.reason, str(error))
            rows=rank_horn_drivers(str(csv),c.candidate_id,c.throat_radius,drivers,target,top_n=max(1,len(drivers)),length=c.length)
            for row in rows:
                row.update({'profile':c.profile,'throat_radius':c.throat_radius,'mouth_radius':c.mouth_radius,'length':c.length})
                (all_rows if row['model_feasible'] else rejected).append(row)
            score=max((objective(row) for row in rows),default=REJECTED_GEOMETRY_OBJECTIVE)
            score=max(score, REJECTED_GEOMETRY_OBJECTIVE)
            cache[key]=score
            return score
        bounds={'mouth_radius':tuple(design['mouth_radius_range']),'length':tuple(design['length_range']),
                'throat_radius':tuple(design.get('throat_radius_range') or (min(screen['throat_radii_m']),max(screen['throat_radii_m'])))}
        # The seed can lie outside the searched box (a least-infeasible seed from
        # a different tranche); include it so the box contains the start point.
        bounds={k:(min(lo,getattr(seed,k)),max(hi,getattr(seed,k))) for k,(lo,hi) in bounds.items()}
        widenable, hard_limits = widening_policy(design, screen)
        best,audit=refine(seed,evaluate,bounds,budget=a.budget,widenable=widenable,hard_limits=hard_limits)
        audit.update(domain_rejected_candidates=domain_rejected, geometry_rejected_candidates=geometry_rejected,
                     fem_evaluations=audit['new_evaluations']-len(domain_rejected)-len(geometry_rejected),
                     seed_candidate=r['horn_label'], seed_source=seed_source,
                     objective='composite score when feasible, else minus infeasibility (search.search_objective)')
        data.update(results=sorted(all_rows,key=lambda r:r['composite_score'],reverse=True),rejected=rejected,refinement=audit)
        data['status']='experimental_candidates' if all_rows else 'no_feasible_design'
        if isinstance(data.get('scoring'), dict):
            data['scoring']['refinement'] = {k: audit[k] for k in ('seed_candidate', 'seed_source', 'new_evaluations',
                                                                  'budget', 'widened_bounds', 'bound_sides_reached')}
        data['total_scored']+=audit['fem_evaluations']*len(drivers)
        data['total_candidates']+=audit['fem_evaluations']
    from horn_analysis.search import annotate_comparable_candidates
    data['results'] = annotate_comparable_candidates(sorted(data['results'], key=lambda row: row['composite_score'], reverse=True)[:a.top_n])
    (out/'search_audit.json').write_text(json.dumps(data['refinement'],indent=2))
    (out/'ranked_results.json').write_text(json.dumps(data,indent=2))


if __name__=='__main__':main()
