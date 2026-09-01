"""Tests for LEM/Webster prescreening orchestrator."""

import json
from pathlib import Path

import numpy as np
import pytest

from horn_core.candidates import CandidateGeometry
from horn_core.parameters import DriverParameters
from horn_analysis.lem_prescreen import (
    lem_prescreen_candidates,
    _spl_from_pressure,
    _load_candidates_csv,
    _write_filtered_csv,
)


def _make_driver(
    driver_id="drv1",
    fs_hz=500.0,
    qes=0.4,
    qms=5.0,
    sd_m2=0.0008,
    re_ohm=6.0,
    bl_tm=8.0,
    mms_kg=0.003,
    le_h=0.0005,
):
    return DriverParameters(
        driver_id=driver_id,
        manufacturer="Test",
        model_name=driver_id,
        fs_hz=fs_hz,
        re_ohm=re_ohm,
        bl_tm=bl_tm,
        sd_m2=sd_m2,
        mms_kg=mms_kg,
        le_h=le_h,
        qms=qms,
        qes=qes,
    )


def _make_candidate(cid, profile="conical", throat=0.025, mouth=0.15, length=0.3):
    return CandidateGeometry(
        candidate_id=cid,
        profile=profile,
        throat_radius=throat,
        mouth_radius=mouth,
        length=length,
    )


class TestSplFromPressure:
    def test_reference_pressure(self):
        """20 µPa should give 0 dB SPL."""
        p = np.array([20e-6])
        spl = _spl_from_pressure(p)
        assert spl[0] == pytest.approx(0.0, abs=0.01)

    def test_1_pa(self):
        """1 Pa should give 94 dB SPL."""
        p = np.array([1.0])
        spl = _spl_from_pressure(p)
        assert spl[0] == pytest.approx(94.0, abs=0.1)


class TestLemPrescreenCandidates:
    @pytest.fixture
    def driver(self):
        return _make_driver()

    @pytest.fixture
    def candidates(self):
        return [
            _make_candidate("c1", "conical", 0.025, 0.15, 0.3),
            _make_candidate("c2", "exponential", 0.025, 0.15, 0.3),
            _make_candidate("c3", "hyperbolic", 0.025, 0.15, 0.3),
        ]

    def test_output_keys(self, candidates, driver):
        result = lem_prescreen_candidates(
            candidates=candidates,
            drivers=[driver],
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=30,
            top_n=2,
        )
        assert "total_evaluated" in result
        assert "total_pairs" in result
        assert "top_n" in result
        assert "filtered_candidate_ids" in result
        assert "rankings" in result

    def test_top_n_respected(self, candidates, driver):
        result = lem_prescreen_candidates(
            candidates=candidates,
            drivers=[driver],
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=30,
            top_n=2,
        )
        assert len(result["filtered_candidate_ids"]) <= 2

    def test_total_evaluated_matches_input(self, candidates, driver):
        result = lem_prescreen_candidates(
            candidates=candidates,
            drivers=[driver],
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=30,
            top_n=10,
        )
        assert result["total_evaluated"] == 3

    def test_total_pairs(self, candidates, driver):
        """3 candidates × 1 driver = 3 pairs."""
        result = lem_prescreen_candidates(
            candidates=candidates,
            drivers=[driver],
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=30,
            top_n=10,
        )
        assert result["total_pairs"] == 3

    def test_multiple_drivers(self, candidates):
        drivers = [_make_driver("d1"), _make_driver("d2")]
        result = lem_prescreen_candidates(
            candidates=candidates,
            drivers=drivers,
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=30,
            top_n=10,
        )
        assert result["total_pairs"] == 6  # 3 candidates × 2 drivers

    def test_filtered_ids_are_valid(self, candidates, driver):
        result = lem_prescreen_candidates(
            candidates=candidates,
            drivers=[driver],
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=30,
            top_n=10,
        )
        valid_ids = {c.candidate_id for c in candidates}
        for cid in result["filtered_candidate_ids"]:
            assert cid in valid_ids

    def test_rankings_sorted_by_score(self, candidates, driver):
        result = lem_prescreen_candidates(
            candidates=candidates,
            drivers=[driver],
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=30,
            top_n=10,
        )
        scores = [r["composite_score"] for r in result["rankings"]]
        assert scores == sorted(scores, reverse=True)

    def test_ranking_entries_have_expected_keys(self, candidates, driver):
        result = lem_prescreen_candidates(
            candidates=candidates,
            drivers=[driver],
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=30,
            top_n=10,
        )
        for entry in result["rankings"]:
            assert "candidate_id" in entry
            assert "driver_id" in entry
            assert "composite_score" in entry
            assert "profile" in entry


class TestKnownGoodGeometry:
    """A well-matched horn should rank higher than a poorly-matched one."""

    def test_good_beats_bad(self):
        driver = _make_driver(fs_hz=300, sd_m2=0.0008)
        good = _make_candidate("good", "exponential", 0.025, 0.15, 0.3)
        bad = _make_candidate("bad", "conical", 0.025, 0.03, 0.05)  # tiny horn

        result = lem_prescreen_candidates(
            candidates=[good, bad],
            drivers=[driver],
            target_f_low=500,
            target_f_high=4000,
            sim_freq_range=(354, 5657),
            num_frequencies=50,
            top_n=2,
        )
        # The well-matched horn should appear first in filtered IDs
        assert result["filtered_candidate_ids"][0] == "good"


class TestCsvIO:
    def test_roundtrip(self, tmp_path):
        candidates = [
            _make_candidate("c1", "conical"),
            _make_candidate("c2", "exponential"),
            _make_candidate("c3", "hyperbolic"),
        ]
        csv_path = str(tmp_path / "candidates.csv")
        # Write
        import csv
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["candidate_id", "profile", "throat_radius", "mouth_radius", "length"])
            for c in candidates:
                writer.writerow([c.candidate_id, c.profile, c.throat_radius, c.mouth_radius, c.length])

        loaded = _load_candidates_csv(csv_path)
        assert len(loaded) == 3
        assert loaded[0].candidate_id == "c1"

    def test_filtered_csv(self, tmp_path):
        candidates = [
            _make_candidate("c1", "conical"),
            _make_candidate("c2", "exponential"),
            _make_candidate("c3", "hyperbolic"),
        ]
        out_path = str(tmp_path / "filtered.csv")
        _write_filtered_csv(candidates, ["c1", "c3"], out_path)

        loaded = _load_candidates_csv(out_path)
        assert len(loaded) == 2
        ids = {c.candidate_id for c in loaded}
        assert ids == {"c1", "c3"}
