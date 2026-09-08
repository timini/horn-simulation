"""Deterministic bounded refinement and exhaustive-screening audit metrics."""
from dataclasses import replace
import math


def geometry_key(c):
    return (c.profile, round(c.throat_radius,8), round(c.mouth_radius,8), round(c.length,8))


def refine(seed, evaluate, bounds, budget=6, relative_step=.1, min_improvement=.001):
    """Coordinate search with a hard evaluation budget; callback returns a score.

    The baseline is evaluated once, outside the new-geometry budget. Callers
    should cache it. No randomness; shrinking steps stop at 1% of seed values.
    """
    if budget < 0 or relative_step <= 0 or min_improvement < 0:
        raise ValueError("Invalid refinement controls")
    best=seed;best_score=float(evaluate(seed));seen={geometry_key(seed)};history=[]
    if not math.isfinite(best_score): raise ValueError("Non-finite baseline score")
    step=relative_step
    while len(history)<budget and step>=.01:
        improved=False
        for field in ('mouth_radius','length','throat_radius'):
            for sign in (-1,1):
                value=getattr(best,field)*(1+sign*step)
                lo,hi=bounds[field]
                if not lo<=value<=hi: continue
                trial=replace(best,**{field:round(value,8)},candidate_id=f'refine_{best.profile}_{len(history):04d}')
                if trial.throat_radius>=trial.mouth_radius or geometry_key(trial) in seen: continue
                seen.add(geometry_key(trial));score=float(evaluate(trial))
                if not math.isfinite(score): raise ValueError("Non-finite refinement score")
                history.append({'candidate_id':trial.candidate_id,'score':score,'geometry':geometry_key(trial)})
                if score>best_score+min_improvement:
                    best,best_score,improved=trial,score,True
                if len(history)>=budget: break
            if len(history)>=budget: break
        if not improved: step/=2
    at_bounds=[key for key,(lo,hi) in bounds.items() if math.isclose(getattr(best,key),lo,rel_tol=1e-6) or math.isclose(getattr(best,key),hi,rel_tol=1e-6)]
    return best, {'best_score':best_score,'new_evaluations':len(history),'budget':budget,
                  'budget_exhausted':len(history)>=budget,'at_bounds':at_bounds,'history':history}


def screening_audit(exhaustive_rows, shortlisted_ids, top_k=10):
    feasible=sorted((r for r in exhaustive_rows if r['model_feasible']),key=lambda r:r['composite_score'],reverse=True)
    if not feasible:
        return {'status':'no_feasible_reference','passed':False}
    selected=[r for r in feasible if r['horn_label'] in set(shortlisted_ids)]
    top=feasible[:top_k]
    recall=sum(r['horn_label'] in set(shortlisted_ids) for r in top)/len(top)
    regret=feasible[0]['composite_score']-selected[0]['composite_score'] if selected else None
    return {'status':'evaluated','feasible_top_k_recall':recall,'score_regret':regret,
            'winner_retained':feasible[0]['horn_label'] in set(shortlisted_ids),
            'passed':regret is not None and regret<=.02 and recall>=.9}


def shortlist_geometries(rows, budget=10, score_margin=.02):
    """Keep the best scores plus two geometrically diverse near-tied candidates.

    The score margin is a conservative search policy, not a measured uncertainty
    interval. A fixed total budget is never exceeded. Hard-rejected geometries
    cannot displace feasible ones; all-rejected searches still diagnose a sample.
    """
    if budget < 1 or score_margin < 0:
        raise ValueError("Invalid shortlist controls")
    rows = [r for r in rows if r.get("simulation_eligible", True)]
    best = {}
    for row in sorted(rows, key=lambda r: (-r['composite_score'], r['candidate_id'], r['driver_id'])):
        best.setdefault(row['candidate_id'], row)
    ordered = sorted(best.values(), key=lambda r: (not r['model_feasible'], -r['composite_score'], r['candidate_id']))
    if len(ordered) <= budget:
        return [r['candidate_id'] for r in ordered]
    # Preserve at least 80% by score and the highest-scoring candidate always.
    elite_count = max(1, math.ceil(.8*budget))
    selected = ordered[:elite_count]
    remaining = ordered[elite_count:]
    fields = ('throat_radius', 'mouth_radius', 'length')
    spans = {k: max(math.log(r[k]) for r in ordered)-min(math.log(r[k]) for r in ordered) for k in fields}
    while len(selected) < budget and remaining:
        cutoff = remaining[0]['composite_score']-score_margin
        pool = [r for r in remaining if r['model_feasible']==remaining[0]['model_feasible'] and r['composite_score'] >= cutoff]
        def distance(r):
            return min((r['profile'] != s['profile']) + sum(((math.log(r[k]/s[k])/spans[k]) if spans[k] else 0)**2 for k in fields) for s in selected)
        chosen = max(pool, key=lambda r: (distance(r), r['composite_score']))
        selected.append(chosen)
        remaining.remove(chosen)
    return [r['candidate_id'] for r in selected]


def annotate_comparable_candidates(rows, score_margin=.02):
    """Mark near-ties without implying calibrated physical uncertainty."""
    best = max((r['composite_score'] for r in rows if r.get('model_feasible')), default=0.)
    return [{**r, 'score_gap_from_best':best-r['composite_score'],
             'comparison_status':('near_tie' if r.get('model_feasible') and best-r['composite_score']<=score_margin else 'lower_ranked'),
             'comparison_scope':'0.02 score comparison margin; physical uncertainty unquantified'} for r in rows]
