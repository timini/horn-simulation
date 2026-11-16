"""Deterministic bounded refinement and exhaustive-screening audit metrics."""
from dataclasses import replace
import math


# Reasons whose size can be measured from a ranking row: (field, limit source).
_MEASURED_REASONS = {"target_band_ripple_exceeded": "passband_ripple_db"}


def infeasibility(row, max_ripple_db=6.0):
    """How far a rejected row is from feasible; 0 for a feasible row.

    Each failed constraint counts 1, plus its measured excess where the row
    holds it (ripple: dB over the limit divided by the limit). Rows without
    metrics (geometry or domain rejections) are infinitely infeasible.
    """
    if row.get("model_feasible"):
        return 0.0
    reasons = row.get("rejection_reasons") or []
    if row.get("passband_ripple_db") is None or not reasons:
        return math.inf
    total = float(len(reasons))
    if "target_band_ripple_exceeded" in reasons:
        total += max(0.0, row["passband_ripple_db"] - max_ripple_db) / max_ripple_db
    return total


def search_objective(row, max_ripple_db=6.0):
    """Score for search decisions: the composite score when feasible, else -infeasibility.

    Every feasible design (score >= 0) beats every infeasible one, and among
    infeasible designs the one nearest to feasible wins, so a search with no
    feasible candidate still moves towards one.
    """
    if row.get("model_feasible"):
        return float(row["composite_score"])
    return -infeasibility(row, max_ripple_db)


def geometry_key(c):
    return (c.profile, round(c.throat_radius,8), round(c.mouth_radius,8), round(c.length,8))


def refine(seed, evaluate, bounds, budget=6, relative_step=.1, min_improvement=.001,
           widenable=None, hard_limits=None, widen_factor=1.25, max_widenings=2):
    """Coordinate search with a hard evaluation budget; callback returns a score.

    The baseline is evaluated once, outside the new-geometry budget. Callers
    should cache it. No randomness; shrinking steps stop at 1% of seed values.

    ``widenable`` maps a field to the sides of its bound that came from a
    heuristic, e.g. ``{'length': ('high',)}``. When the search stops with the
    best design on such a side and budget remains, that side moves outward by
    ``widen_factor`` (never past ``hard_limits``) and the search resumes, up to
    ``max_widenings`` times. User-set sides are never moved.
    """
    if budget < 0 or relative_step <= 0 or min_improvement < 0 or widen_factor <= 1 or max_widenings < 0:
        raise ValueError("Invalid refinement controls")
    bounds = {k: tuple(v) for k, v in bounds.items()}
    widenable = widenable or {}
    hard_limits = hard_limits or {}
    best=seed;best_score=float(evaluate(seed));seen={geometry_key(seed)};history=[]
    if not math.isfinite(best_score): raise ValueError("Non-finite baseline score")
    widenings = []
    def at_bounds(candidate):
        return [(key, side) for key, (lo, hi) in bounds.items()
                for side, edge in (("low", lo), ("high", hi)) if math.isclose(getattr(candidate, key), edge, rel_tol=1e-6)]
    while True:
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
        if len(history) >= budget or len(widenings) >= max_widenings:
            break
        moved = []
        for key, side in at_bounds(best):
            if side not in widenable.get(key, ()):
                continue
            lo, hi = bounds[key]
            limit_lo, limit_hi = hard_limits.get(key, (0.0, math.inf))
            if side == "high":
                new = (lo, min(hi * widen_factor, limit_hi))
            else:
                new = (max(lo / widen_factor, limit_lo), hi)
            if new != (lo, hi):
                bounds[key] = new
                moved.append({"field": key, "side": side, "from": [lo, hi], "to": list(new)})
        if not moved:
            break
        widenings.extend(moved)
    return best, {'best_score':best_score,'new_evaluations':len(history),'budget':budget,
                  'budget_exhausted':len(history)>=budget,
                  'at_bounds':sorted({key for key, _ in at_bounds(best)}),
                  'bound_sides_reached':[f"{key}_{side}" for key, side in at_bounds(best)],
                  'widened_bounds':widenings,'final_bounds':{k:list(v) for k,v in bounds.items()},
                  'history':history}


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


SHORTLIST_POLICY = ("best design of every profile in the leading feasibility class, then by search objective "
                    "with at most half the budget per profile while other profiles have candidates of the same feasibility, last 20% "
                    "geometrically diverse within 0.02 score; fixed total budget")


def shortlist_geometries(rows, budget=10, score_margin=.02, max_ripple_db=6.0, profile_share=.5):
    """Pick the geometries that go to FEM, with a fixed total budget.

    Order is by ``search_objective``: feasible designs by score, then
    infeasible ones nearest to feasible, so an all-infeasible screen still
    sends its most promising geometries. Profiles are explored fairly: the
    best design of each profile goes first (only profiles whose best design is
    in the same feasibility class as the overall best), and no profile takes
    more than ``profile_share`` of the budget while another profile still has
    candidates. The last 20% of slots prefer geometrically different
    near-ties. The score margin is a conservative search policy, not a
    measured uncertainty interval. Hard-rejected geometries without metrics
    (simulation_eligible False) are never sent.
    """
    if budget < 1 or score_margin < 0 or not 0 < profile_share <= 1:
        raise ValueError("Invalid shortlist controls")
    rows = [r for r in rows if r.get("simulation_eligible", True)]
    objective = lambda r: search_objective(r, max_ripple_db)
    best = {}
    for row in sorted(rows, key=lambda r: (-objective(r), r['candidate_id'], r['driver_id'])):
        best.setdefault(row['candidate_id'], row)
    ordered = sorted(best.values(), key=lambda r: (-objective(r), r['candidate_id']))
    if len(ordered) <= budget:
        return [r['candidate_id'] for r in ordered]
    lead_class = ordered[0]['model_feasible']
    selected, profiles = [], set()
    for row in ordered:
        if row['profile'] not in profiles and row['model_feasible'] == lead_class and len(selected) < budget:
            selected.append(row)
            profiles.add(row['profile'])
    remaining = [r for r in ordered if r not in selected]
    cap = max(1, math.ceil(profile_share * budget))

    def allowed(row):
        count = sum(s['profile'] == row['profile'] for s in selected)
        if count < cap:
            return True
        # Infeasible geometries never displace feasible ones for diversity.
        return not any(r['profile'] != row['profile'] and r['model_feasible'] == row['model_feasible']
                       and sum(s['profile'] == r['profile'] for s in selected) < cap for r in remaining)

    elite_count = max(len(selected), math.ceil(.8*budget))
    while len(selected) < elite_count and remaining:
        chosen = next((r for r in remaining if allowed(r)), remaining[0])
        selected.append(chosen)
        remaining.remove(chosen)
    fields = ('throat_radius', 'mouth_radius', 'length')
    spans = {k: max(math.log(r[k]) for r in ordered)-min(math.log(r[k]) for r in ordered) for k in fields}
    while len(selected) < budget and remaining:
        head = next((r for r in remaining if allowed(r)), remaining[0])
        cutoff = objective(head)-score_margin
        pool = [r for r in remaining if r['model_feasible'] == head['model_feasible'] and objective(r) >= cutoff and allowed(r)] or [head]
        def distance(r):
            return min((r['profile'] != s['profile']) + sum(((math.log(r[k]/s[k])/spans[k]) if spans[k] else 0)**2 for k in fields) for s in selected)
        chosen = max(pool, key=lambda r: (distance(r), objective(r)))
        selected.append(chosen)
        remaining.remove(chosen)
    return [r['candidate_id'] for r in sorted(selected, key=lambda r: (-objective(r), r['candidate_id']))]


def annotate_comparable_candidates(rows, score_margin=.02):
    """Mark near-ties without implying calibrated physical uncertainty."""
    best = max((r['composite_score'] for r in rows if r.get('model_feasible')), default=0.)
    return [{**r, 'score_gap_from_best':best-r['composite_score'],
             'comparison_status':('near_tie' if r.get('model_feasible') and best-r['composite_score']<=score_margin else 'lower_ranked'),
             'comparison_scope':'0.02 score comparison margin; physical uncertainty unquantified'} for r in rows]
