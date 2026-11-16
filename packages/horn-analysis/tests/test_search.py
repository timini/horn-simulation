from horn_core.candidates import CandidateGeometry
from horn_analysis.search import refine, screening_audit, geometry_key
import pytest


def test_refinement_respects_budget_bounds_and_preserves_best():
    seed = CandidateGeometry('seed','conical',.02,.10,.20)
    bounds = {'throat_radius':(.02,.02),'mouth_radius':(.08,.15),'length':(.1,.3)}
    calls = []
    def evaluate(c):
        calls.append(geometry_key(c))
        assert all(lo <= getattr(c,field) <= hi for field,(lo,hi) in bounds.items())
        return 1 - (c.mouth_radius-.14)**2 - (c.length-.24)**2
    best, audit = refine(seed, evaluate, bounds, budget=7)
    assert audit['new_evaluations'] <= 7
    assert len(calls) == len(set(calls)) == 1+audit['new_evaluations']
    assert audit['best_score'] >= evaluate(seed)
    again, repeated = refine(seed, lambda c: 1-(c.mouth_radius-.14)**2-(c.length-.24)**2,bounds,budget=7)
    assert best == again and audit == repeated


def test_screening_gate_counts_feasible_pairs_and_winner():
    rows = [{'horn_label':str(i),'model_feasible':True,'composite_score':1-i/100} for i in range(12)]
    rows.append({'horn_label':'bad','model_feasible':False,'composite_score':100})
    assert screening_audit(rows,[str(i) for i in range(10)])['passed']
    missed=screening_audit(rows,[str(i) for i in range(3,12)])
    assert not missed['passed'] and not missed['winner_retained']
    assert missed['score_regret'] == pytest.approx(.03)


def test_shortlist_preserves_winner_and_explores_near_ties_with_fixed_budget():
    from horn_analysis.search import shortlist_geometries
    rows = [dict(candidate_id=str(i), driver_id='motor', profile='conical',
                 throat_radius=.02, mouth_radius=.10+i*.00001, length=.1,
                 model_feasible=True, composite_score=.9-i*.001) for i in range(12)]
    rows[-1].update(profile='exponential', mouth_radius=.2, length=.2)
    selected = shortlist_geometries(rows,5)
    assert len(selected) == len(set(selected)) == 5
    # Best scores first; no profile takes more than half the budget (3 of 5)
    # while another profile still has candidates.
    assert selected[:3] == ['0','1','2']
    assert '11' in selected
    assert sum(not s.isdigit() or s != '11' for s in selected) == 4
    assert shortlist_geometries(list(reversed(rows)),5) == selected
    rows[-1]['model_feasible'] = False
    assert '11' not in shortlist_geometries(rows,5)
def test_refined_curve_ids_preserve_profile_grouping():
    from horn_analysis.html_report import _pick_representative_per_profile
    seed = CandidateGeometry('auto_exp_0001', 'exponential', .02, .10, .20)
    _, audit = refine(seed, lambda c: 1., {'throat_radius': (.02, .02),
                       'mouth_radius': (.08, .15), 'length': (.1, .3)}, budget=2)
    curves = {'auto_exp_0001': 'original.csv',
              **{row['candidate_id']: row['candidate_id'] + '.csv' for row in audit['history']}}
    assert list(_pick_representative_per_profile(curves)) == ['exponential']


def _row(cid, profile, score, feasible=True, ripple=3., reasons=None, **geometry):
    return dict(candidate_id=cid, driver_id='motor', profile=profile,
                throat_radius=geometry.get('throat_radius', .02), mouth_radius=geometry.get('mouth_radius', .1),
                length=geometry.get('length', .2), model_feasible=feasible,
                composite_score=score if feasible else 0., passband_ripple_db=ripple,
                rejection_reasons=[] if feasible else (reasons or ['target_band_ripple_exceeded']))


def test_search_objective_puts_feasible_first_then_least_infeasible():
    from horn_analysis.search import search_objective
    feasible = _row('a', 'os', .1)
    near = _row('b', 'os', 0, feasible=False, ripple=6.5)
    far = _row('c', 'os', 0, feasible=False, ripple=12.)
    two_reasons = _row('d', 'os', 0, feasible=False, ripple=6.5, reasons=['target_band_ripple_exceeded', 'compression_ratio_exceeded'])
    no_metrics = {**_row('e', 'os', 0, feasible=False), 'passband_ripple_db': None}
    values = [search_objective(r) for r in (feasible, near, far, two_reasons, no_metrics)]
    assert values[0] == .1
    assert values[0] > values[1] > values[2] and values[1] > values[3] > values[4]
    assert values[4] == float('-inf')


def test_shortlist_gives_every_profile_a_slot_and_caps_a_dominant_profile():
    """Live run B sent 8 of 10 FEM geometries to one profile."""
    from collections import Counter
    from horn_analysis.search import shortlist_geometries
    rows = [_row(f'hyp{i:02d}', 'hyperbolic', .9-i*.001, mouth_radius=.1+i*.001) for i in range(12)]
    rows += [_row(f'{p}{i}', p, .8-i*.001, mouth_radius=.1+i*.003) for p in ('conical', 'os', 'exponential') for i in range(3)]
    selected = shortlist_geometries(rows, 10)
    counts = Counter(r['profile'] for r in rows if r['candidate_id'] in selected)
    assert len(selected) == 10 and selected[0] == 'hyp00'
    assert set(counts) == {'hyperbolic', 'conical', 'os', 'exponential'}
    assert counts['hyperbolic'] == 5


def test_shortlist_profile_leaders_must_match_the_leading_feasibility_class():
    from horn_analysis.search import shortlist_geometries
    rows = [_row(f'hyp{i:02d}', 'hyperbolic', .9-i*.001, mouth_radius=.1+i*.001) for i in range(12)]
    rows.append(_row('con_bad', 'conical', 0, feasible=False, ripple=6.1))
    assert 'con_bad' not in shortlist_geometries(rows, 10)


def test_all_infeasible_shortlist_sends_the_least_infeasible_geometries():
    from horn_analysis.search import shortlist_geometries
    # Every score is 0, so the old order fell back to candidate ID.
    rows = [_row(f'c{i:02d}', 'conical', 0, feasible=False, ripple=20.-i, mouth_radius=.1+i*.001) for i in range(12)]
    selected = shortlist_geometries(rows, 3)
    assert selected[0] == 'c11'
    assert {'c11', 'c10'} <= set(selected)


def test_refinement_widens_a_heuristic_bound_it_ends_on_but_never_a_user_bound():
    seed = CandidateGeometry('seed', 'conical', .02, .10, .20)
    bounds = {'throat_radius': (.02, .02), 'mouth_radius': (.10, .10), 'length': (.1, .2)}
    longer_is_better = lambda c: c.length
    best, audit = refine(seed, longer_is_better, bounds, budget=6)
    assert best.length == pytest.approx(.2) and audit['widened_bounds'] == []
    assert 'length_high' in audit['bound_sides_reached']
    best, audit = refine(seed, longer_is_better, bounds, budget=6, widenable={'length': ('high',)},
                         hard_limits={'length': (0., .26)})
    assert .2 < best.length <= .26
    assert audit['widened_bounds'][0]['field'] == 'length' and audit['widened_bounds'][0]['side'] == 'high'
    assert audit['final_bounds']['length'][1] <= .26
    # A user-set side (not listed as widenable) stays put.
    best, audit = refine(seed, longer_is_better, bounds, budget=6, widenable={'length': ('low',)})
    assert best.length == pytest.approx(.2) and audit['widened_bounds'] == []


def test_refinement_seed_falls_back_to_least_infeasible_fem_row():
    from horn_analysis.refine import select_seed
    base = dict(profile='os', throat_radius=.03, mouth_radius=.15, length=.3, simulation_eligible=True)
    rejected = [
        {**base, 'horn_label': 'cad', 'passband_ripple_db': None, 'model_feasible': False,
         'rejection_reasons': ['cad_geometry_rejected'], 'simulation_eligible': False},
        {**base, 'horn_label': 'far', 'passband_ripple_db': 13., 'model_feasible': False,
         'rejection_reasons': ['target_band_ripple_exceeded']},
        {**base, 'horn_label': 'near', 'passband_ripple_db': 7., 'model_feasible': False,
         'rejection_reasons': ['target_band_ripple_exceeded']},
    ]
    row, source = select_seed({'results': [], 'rejected': rejected})
    assert row['horn_label'] == 'near' and source == 'least_infeasible'
    feasible = {**base, 'horn_label': 'ok', 'model_feasible': True, 'composite_score': .5, 'passband_ripple_db': 4.}
    assert select_seed({'results': [feasible], 'rejected': rejected}) == (feasible, 'best_feasible')
    assert select_seed({'results': [], 'rejected': rejected[:1]}) == (None, None)


def test_refinement_widens_only_heuristic_sides_and_keeps_throat_under_ka_cap():
    import math
    from horn_analysis.refine import widening_policy
    design = {'target_f_high': 5000, 'throat_ka_max': 3.8317,
              'bound_sources': {'mouth_radius': {'low': 'heuristic', 'high': 'heuristic'},
                                'length': {'low': 'user', 'high': 'user'},
                                'throat_radius': {'low': 'driver_prescreen', 'high': 'driver_prescreen'}}}
    widenable, limits = widening_policy(design, {'ka_max': 3.8317})
    assert widenable == {'mouth_radius': ('low', 'high'), 'length': (), 'throat_radius': ('low', 'high')}
    assert limits['throat_radius'][1] == pytest.approx(343*3.8317/(2*math.pi*5000))
