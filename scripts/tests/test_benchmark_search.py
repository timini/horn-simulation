"""Analytical stage of the search benchmark cases (no FEM)."""
from collections import Counter
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('benchmark_search', ROOT/'scripts/benchmark_search.py')
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


@pytest.fixture(scope='module')
def wide():
    # The wide-search example specification: 5 x 5 mouth/length grid, 8 FEM slots, no throat radius.
    parameters = json.loads((ROOT/'examples/6nmb420-320-5000/wide-search-specification.json').read_text())['parameters']
    assert parameters['throat_radius'] is None
    assert (parameters['num_mouth_radii'], parameters['num_lengths'], parameters['lem_top_n']) == (5, 5, 8)
    return benchmark.analytical_case('wide_320_5000', num_mouth_radii=5, num_lengths=5)


def test_wide_search_throats_stay_below_the_first_radial_mode(wide):
    cap = 343 * 3.8317 / (2 * math.pi * 5000)
    assert max(wide['prescreen'].throat_radii_m) == pytest.approx(cap, rel=1e-6)
    assert all(c.throat_radius <= cap * (1 + 1e-9) for c in wide['candidates'])


def test_wide_search_shortlist_holds_a_design_within_6_db(wide):
    """The broad search's 8 FEM finalists all failed the ripple limit (11-12 dB)."""
    shortlist = wide['shortlist']
    assert len(shortlist) == 8
    # Without adapter data the cone estimate is applied and flagged; the
    # designs meet 6 dB without it and are kept, not rejected.
    assert all(r['front_chamber_estimated'] for r in shortlist)
    assert all(r['model_feasible'] for r in shortlist)
    assert min(r['without_front_chamber']['passband_ripple_db'] for r in shortlist) <= 6.0
    plugged = benchmark.analytical_case('wide_320_5000', num_mouth_radii=5, num_lengths=5, front_chamber_volume_m3=0.)
    best = min(r['passband_ripple_db'] for r in plugged['shortlist'])
    assert best <= 6.0
    assert not any(r['front_chamber_estimated'] for r in plugged['shortlist'])


def test_wide_search_shortlist_explores_profiles_fairly(wide):
    counts = Counter(r['profile'] for r in wide['shortlist'])
    assert len(counts) >= 5
    assert max(counts.values()) <= 4


def test_size_capped_search_returns_designs_within_the_cap_with_a_stated_shortfall():
    case = benchmark.analytical_case('size_capped_800_1600')
    assert case['candidates'] and case['shortlist']
    assert all(c.length <= .1 + 1e-12 for c in case['candidates'])
    assert min(c.length for c in case['candidates']) == pytest.approx(.04)
    assert any(r['model_feasible'] for r in case['shortlist'])
    compromise = case['design'].size_compromise
    assert compromise['expected_f_low_hz'] == pytest.approx(343 / .4)
    assert compromise['shortfall_octaves'] == pytest.approx(math.log2(857.5 / 800))


def test_benchmark_lists_the_pipeline_cases():
    assert set(benchmark.PIPELINE_CASES) == {'wide_320_5000', 'size_capped_800_1600'}
    assert benchmark.PIPELINE_CASES['size_capped_800_1600']['max_length'] == .1
    assert 'min_length' not in benchmark.PIPELINE_CASES['size_capped_800_1600']
