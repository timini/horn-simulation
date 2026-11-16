"""An estimated front chamber shapes the score but cannot by itself reject a design."""
from dataclasses import replace
import json
from pathlib import Path

import pytest

from horn_core.candidates import CandidateGeometry
from horn_core.driver_interface import estimate_cone_front_chamber_m3
from horn_drivers.loader import load_drivers, load_driver
from horn_analysis.lem_prescreen import lem_prescreen_candidates
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec, weighted_score

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "examples/6nmb420-320-5000"
# A 320-5000 Hz CD horn at the ka = 3.83 throat cap. Analytically it has about
# 5.7 dB ripple without a front chamber and about 7.1 dB with the 139 cm3 cone
# estimate, so only the estimate breaks the 6 dB limit.
HORN = CandidateGeometry("cd_cap", "cd", 0.0418346, 0.1446, 0.3)


def _screen(driver):
    target = TargetSpec(320, 5000)
    result = lem_prescreen_candidates([HORN], [driver], 320, 5000, (320 / 2**.5, 5000 * 2**.5),
                                      num_frequencies=161, top_n=1, radiation_model="modal_baffled", target=target)
    return result["rankings"][0], target


def _driver():
    return load_drivers(str(EXAMPLE / "driver-database.json"))[0]


def test_estimated_chamber_failure_is_flagged_not_rejected():
    row, target = _screen(_driver())
    assert row["front_chamber_estimated"] is True
    assert row["passband_ripple_db"] > 6.0  # levels and score include the estimate
    assert row["without_front_chamber"]["passband_ripple_db"] < 6.0
    assert row["estimated_chamber_only_reasons"] == ["target_band_ripple_exceeded"]
    assert row["rejection_reasons"] == [] and row["model_feasible"] is True
    assert "response_limited_by_estimated_front_chamber" in row["evidence_gaps"]
    assert "front_chamber_volume_estimated" in row["evidence_gaps"]
    assert row["composite_score"] == pytest.approx(weighted_score(target, row["score_terms"]))
    assert row["composite_score"] > 0


def test_supplied_chamber_volume_overrides_the_estimate_and_can_reject():
    driver = _driver()
    supplied = replace(driver, front_chamber_volume_m3=estimate_cone_front_chamber_m3(driver.sd_m2),
                       front_chamber_area_m2=driver.sd_m2)
    row, _ = _screen(supplied)
    assert row["front_chamber_estimated"] is False
    assert row["model_feasible"] is False
    assert row["rejection_reasons"] == ["target_band_ripple_exceeded"]
    estimated, _ = _screen(driver)
    assert row["passband_ripple_db"] == pytest.approx(estimated["passband_ripple_db"], abs=1e-9)
    plugged, _ = _screen(replace(driver, front_chamber_volume_m3=0.))
    assert plugged["front_chamber_estimated"] is False and plugged["model_feasible"] is True
    assert plugged["passband_ripple_db"] == pytest.approx(estimated["without_front_chamber"]["passband_ripple_db"], abs=1e-9)


def test_failure_without_the_chamber_still_rejects():
    driver = _driver()
    row, _ = _screen(driver)
    strict = TargetSpec(320, 5000, max_ripple_db=3.0)
    result = lem_prescreen_candidates([HORN], [driver], 320, 5000, (320 / 2**.5, 5000 * 2**.5),
                                      num_frequencies=161, top_n=1, radiation_model="modal_baffled", target=strict)
    strict_row = result["rankings"][0]
    assert strict_row["model_feasible"] is False
    assert strict_row["rejection_reasons"] == ["target_band_ripple_exceeded"]
    assert strict_row["estimated_chamber_only_reasons"] == []


def test_fem_ranking_applies_the_same_policy_to_the_preserved_800_1600_example():
    """The integration branch reported this example flipping to rejected under the estimate."""
    example = ROOT / "examples/800-1600-modal"
    preserved = json.loads((example / "ranking.json").read_text())[0]
    driver = load_driver(str(ROOT / "data/drivers-curated"), preserved["driver_id"])
    row = rank_horn_drivers(str(example / "horn-impedance.csv"), preserved["horn_label"], preserved["throat_radius"],
                            [driver], TargetSpec(800, 1600), top_n=1, length=preserved["length"])[0]
    assert row["front_chamber_estimated"] is True
    assert row["passband_ripple_db"] > 6.0
    assert row["without_front_chamber"]["passband_ripple_db"] <= 6.0
    assert row["model_feasible"] is True
    assert row["estimated_chamber_only_reasons"] == ["target_band_ripple_exceeded"]


def test_wide_throat_is_flagged_above_first_asymmetric_mode():
    row, _ = _screen(replace(_driver(), front_chamber_volume_m3=0.))
    assert "throat_above_first_asymmetric_mode" in row["evidence_gaps"]
    small = CandidateGeometry("small", "conical", 0.01, 0.1446, 0.3)
    result = lem_prescreen_candidates([small], [replace(_driver(), front_chamber_volume_m3=0.)], 320, 5000,
                                      (320 / 2**.5, 5000 * 2**.5), num_frequencies=161, top_n=1,
                                      radiation_model="modal_baffled", target=TargetSpec(320, 5000))
    assert "throat_above_first_asymmetric_mode" not in result["rankings"][0]["evidence_gaps"]
