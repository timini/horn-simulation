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
    assert selected[:4] == ['0','1','2','3']
    assert '11' in selected
    assert shortlist_geometries(list(reversed(rows)),5) == selected
    rows[-1]['model_feasible'] = False
    assert '11' not in shortlist_geometries(rows,5)
