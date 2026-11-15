"""Direct tests of the production ranking score in evaluate_response."""
import numpy as np
import pytest

from horn_analysis.evaluation import evaluate_response
from horn_analysis.kpi import format_kpi
from horn_analysis.scoring import TargetSpec, describe_scoring, exchange_rates, parse_score_weights

F = np.geomspace(500 / np.sqrt(2), 4000 * np.sqrt(2), 97)
TARGET = TargetSpec(500, 4000)


def response(level=100., slope_below_db_per_oct=12., slope_above_db_per_oct=12., ripple_db=0.):
    """Flat in-band level, linear roll-off in dB per octave outside the band.

    ``ripple_db`` adds a smooth in-band undulation of that peak-to-peak size.
    """
    octave = np.log2(F)
    y = np.full_like(F, level)
    below, above = F < 500, F > 4000
    y[below] -= slope_below_db_per_oct * (np.log2(500) - octave[below])
    y[above] -= slope_above_db_per_oct * (octave[above] - np.log2(4000))
    inside = ~below & ~above
    phase = (octave[inside] - np.log2(500)) / 3
    y[inside] += ripple_db / 2 * np.cos(2 * np.pi * phase)
    return y


def test_below_band_extension_changes_score_of_equal_in_band_designs():
    steep = evaluate_response(F, response(slope_below_db_per_oct=24), TARGET)
    gentle = evaluate_response(F, response(slope_below_db_per_oct=6), TARGET)
    assert steep["model_feasible"] and gentle["model_feasible"]
    assert steep["passband_ripple_db"] == pytest.approx(gentle["passband_ripple_db"])
    assert steep["avg_sensitivity_db"] == pytest.approx(gentle["avg_sensitivity_db"])
    # 3 dB below the band mean: 24 dB/oct reaches it 1/8 octave below 500 Hz,
    # 6 dB/oct half an octave below (the sweep edge).
    assert steep["f3_low_hz"] == pytest.approx(500 * 2**-0.125, rel=1e-3)
    assert gentle["f3_low_hz"] == pytest.approx(500 * 2**-0.5, rel=1e-3)
    assert gentle["extension_octaves"] > steep["extension_octaves"]
    assert gentle["composite_score"] > steep["composite_score"]


@pytest.mark.parametrize("flatter_by_db", [1., 3.])
def test_much_more_sensitive_design_is_not_beaten_by_slightly_flatter_one(flatter_by_db):
    flat = evaluate_response(F, response(level=90., ripple_db=1.), TARGET)
    loud = evaluate_response(F, response(level=105., ripple_db=1. + flatter_by_db), TARGET)
    assert flat["model_feasible"] and loud["model_feasible"]
    assert loud["passband_ripple_db"] - flat["passband_ripple_db"] == pytest.approx(flatter_by_db, abs=0.05)
    assert loud["composite_score"] > flat["composite_score"]


def test_default_exchange_rates_are_documented_values():
    rates = exchange_rates(TARGET.score_weights)
    assert rates["sensitivity_db_per_db_of_ripple"] == pytest.approx(2.)
    assert rates["sensitivity_db_per_octave_of_extension"] == pytest.approx(6.)
    assert rates["sensitivity_db_per_doubling_of_horn_volume"] == 0.
    # 1 dB of extra ripple costs exactly as much score as 2 dB of sensitivity.
    a = evaluate_response(F, response(level=100., ripple_db=2.), TARGET)
    b = evaluate_response(F, response(level=98., ripple_db=1.), TARGET)
    assert a["composite_score"] == pytest.approx(b["composite_score"], abs=1e-4)


def test_order_follows_the_weights():
    flat = response(level=95., ripple_db=0.5)
    loud = response(level=100., ripple_db=5.)
    for weights, winner in [("ripple=1,sensitivity=0.1", "flat"), ("ripple=0.1,sensitivity=1", "loud")]:
        target = TargetSpec(500, 4000, score_weights=weights)
        scores = {"flat": evaluate_response(F, flat, target)["composite_score"],
                  "loud": evaluate_response(F, loud, target)["composite_score"]}
        assert max(scores, key=scores.get) == winner


def test_size_term_ranks_smaller_horn_only_when_weighted():
    y = response()
    for weights, expect_smaller_wins in [(None, False), ("size=0.2", True)]:
        target = TargetSpec(500, 4000, score_weights=weights)
        short = evaluate_response(F, y, target, length=0.2, mouth_radius=0.1)
        long = evaluate_response(F, y, target, length=0.4, mouth_radius=0.1)
        if expect_smaller_wins:
            assert short["composite_score"] > long["composite_score"]
        else:
            assert short["composite_score"] == long["composite_score"]


def test_f3_uses_band_mean_not_a_narrow_peak():
    y = response(slope_below_db_per_oct=12, slope_above_db_per_oct=12)
    peak = np.argmin(abs(F - 1500))
    y[peak] += 4.
    result = evaluate_response(F, y, TargetSpec(500, 4000, max_ripple_db=6))
    # A peak-referenced f3 would collapse onto the 1.5 kHz lobe.
    assert result["f3_low_hz"] < 500 and result["f3_high_hz"] > 4000
    assert result["f3_reference"] == "target_band_mean_minus_3db"


def test_response_already_down_at_band_edge_gives_inband_f3_and_negative_extension():
    y = np.where(F < 1000, 100 - 12 * (np.log2(1000) - np.log2(F)), 100.)
    result = evaluate_response(F, y, TargetSpec(500, 4000, max_ripple_db=20))
    assert 500 < result["f3_low_hz"] < 1000
    assert not result["f3_low_is_bound"]
    assert result["extension_octaves"] < 0


def test_flat_sweep_reports_f3_as_sweep_bounds():
    result = evaluate_response(F, np.full_like(F, 100.), TARGET)
    assert result["f3_low_is_bound"] and result["f3_high_is_bound"]
    assert format_kpi(result, "f3_low_hz") == "≤354"
    assert result["extension_octaves"] == pytest.approx(1.)


def test_min_output_floor_is_a_hard_constraint():
    y = response(level=95., ripple_db=2.)
    assert evaluate_response(F, y, TargetSpec(500, 4000, min_output_db=93.5))["model_feasible"]
    rejected = evaluate_response(F, y, TargetSpec(500, 4000, min_output_db=94.5))
    assert rejected["rejection_reasons"] == ["minimum_output_not_met"]
    assert rejected["composite_score"] == 0


def test_score_weights_parse_and_normalise():
    assert parse_score_weights(None) == pytest.approx(TARGET.score_weights)
    weights = parse_score_weights("sensitivity=1, ripple=1, extension=0")
    assert weights == pytest.approx({"sensitivity": .5, "ripple": .5, "extension": 0., "size": 0.})
    assert sum(parse_score_weights({"size": 1.}).values()) == pytest.approx(1.)
    for bad in ["coverage=1", "ripple=-1", "ripple", "sensitivity=0,ripple=0,extension=0,size=0"]:
        with pytest.raises(ValueError):
            parse_score_weights(bad)


def test_scoring_description_separates_user_and_heuristic_bounds():
    target = TargetSpec(500, 4000, max_ripple_db=3., min_output_db=95., score_weights="size=0.1")
    design = {"mouth_radius_range": [.1, .2], "length_range": [.15, .3], "throat_radius_range": [.02, .03]}
    summary = describe_scoring(target, design, {"max_length": .3, "mouth_radius": None})
    assert summary["weights_source"] == "user"
    assert summary["constraints"]["max_ripple_db"] == {"value": 3., "source": "user"}
    assert summary["constraints"]["max_compression_ratio"]["source"] == "pipeline_default"
    assert summary["constraints"]["min_output_db"]["source"] == "user"
    assert summary["geometry_bounds"]["length_m"]["source"] == "heuristic_with_user_limit"
    assert summary["geometry_bounds"]["mouth_radius_m"]["source"] == "heuristic"
    assert summary["geometry_bounds"]["throat_radius_m"]["source"] == "driver_prescreen"
    assert "clip" in summary["formula"]
