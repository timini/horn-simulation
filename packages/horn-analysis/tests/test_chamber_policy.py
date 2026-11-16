"""An estimated cone front chamber is modelled but does not by itself reject a design."""
import json
from dataclasses import replace
from pathlib import Path

import pytest

from horn_core.driver_interface import estimate_cone_front_chamber_m3
from horn_drivers.loader import load_driver
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "examples/800-1600-modal"


@pytest.fixture(scope="module")
def example():
    preserved = json.loads((EXAMPLE / "ranking.json").read_text())[0]
    driver = load_driver(str(ROOT / "data/drivers-curated"), preserved["driver_id"])
    assert driver.driver_type == "cone" and driver.front_chamber_volume_m3 is None
    return preserved, driver


def _rank(preserved, driver):
    return rank_horn_drivers(str(EXAMPLE / "horn-impedance.csv"), preserved["horn_label"], preserved["throat_radius"],
                             [driver], TargetSpec(800, 1600), top_n=1, length=preserved["length"])[0]


def test_estimated_chamber_flags_instead_of_rejecting(example):
    row = _rank(*example)
    # The estimate is applied: ripple is over the 6 dB limit with the chamber.
    assert row["passband_ripple_db"] > 6.0
    assert row["front_chamber_estimated"] is True
    assert row["driver_interface"]["front_chamber_source"] == "estimated"
    assert row["estimated_chamber_reasons"] == ["target_band_ripple_exceeded"]
    assert row["rejection_reasons"] == []
    assert row["model_feasible"] is True
    assert row["composite_score"] > 0
    assert "front_chamber_volume_estimated" in row["evidence_gaps"]


def test_estimated_chamber_still_lowers_the_score(example):
    preserved, driver = example
    estimated = _rank(preserved, driver)
    ideal = _rank(preserved, replace(driver, front_chamber_volume_m3=0.0))
    assert ideal["passband_ripple_db"] == pytest.approx(1.485, abs=1e-3)
    assert ideal["front_chamber_estimated"] is False
    assert estimated["composite_score"] < ideal["composite_score"]


def test_supplied_chamber_rejects_normally(example):
    preserved, driver = example
    supplied = replace(driver, front_chamber_volume_m3=estimate_cone_front_chamber_m3(driver.sd_m2),
                       front_chamber_area_m2=driver.sd_m2)
    row = _rank(preserved, supplied)
    assert row["front_chamber_estimated"] is False
    assert row["rejection_reasons"] == ["target_band_ripple_exceeded"]
    assert row["model_feasible"] is False and row["composite_score"] == 0.0
