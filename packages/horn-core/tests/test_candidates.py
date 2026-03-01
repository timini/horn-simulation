"""Tests for geometry candidate generation."""

import csv

import pytest

from horn_core.candidates import (
    CandidateGeometry,
    generate_candidates,
    write_candidates_csv,
)


class TestGenerateCandidates:
    def test_generates_non_empty_list(self):
        candidates = generate_candidates(
            target_f_low=500, target_f_high=4000,
            max_length=0.5, max_mouth_radius=0.2,
        )
        assert len(candidates) > 0

    def test_all_candidates_valid_geometry(self):
        candidates = generate_candidates(
            target_f_low=500, target_f_high=4000,
            max_length=0.5, max_mouth_radius=0.2,
        )
        for c in candidates:
            assert c.mouth_radius > c.throat_radius, (
                f"{c.candidate_id}: mouth_radius must exceed throat_radius"
            )
            assert c.length > 0
            assert c.length <= 0.5
            assert c.mouth_radius <= 0.2

    def test_respects_max_constraints(self):
        candidates = generate_candidates(
            target_f_low=500, target_f_high=4000,
            max_length=0.3, max_mouth_radius=0.1,
        )
        for c in candidates:
            assert c.length <= 0.3 + 1e-9
            assert c.mouth_radius <= 0.1 + 1e-9

    def test_all_profiles_represented(self):
        candidates = generate_candidates(
            target_f_low=500, target_f_high=4000,
            max_length=0.5, max_mouth_radius=0.2,
        )
        profiles = {c.profile for c in candidates}
        assert profiles == {"conical", "exponential", "hyperbolic", "tractrix", "os", "lecleach", "cd"}

    def test_custom_profiles(self):
        candidates = generate_candidates(
            target_f_low=500, target_f_high=4000,
            max_length=0.5, max_mouth_radius=0.2,
            profiles=["conical", "exponential"],
        )
        profiles = {c.profile for c in candidates}
        assert "hyperbolic" not in profiles

    def test_unique_candidate_ids(self):
        candidates = generate_candidates(
            target_f_low=500, target_f_high=4000,
            max_length=0.5, max_mouth_radius=0.2,
        )
        ids = [c.candidate_id for c in candidates]
        assert len(ids) == len(set(ids))

    def test_expected_grid_size(self):
        """With default throat radii (3) and 4x4 grid across 7 profiles,
        maximum would be 7*3*4*4 = 336 but invalid combos are skipped."""
        candidates = generate_candidates(
            target_f_low=500, target_f_high=4000,
            max_length=0.5, max_mouth_radius=0.2,
        )
        # Some combos filtered (mouth <= throat), so count should be < 336
        assert len(candidates) <= 336
        assert len(candidates) > 100  # but we should still have plenty


class TestWriteCandidatesCsv:
    def test_writes_valid_csv(self, tmp_path):
        candidates = generate_candidates(
            target_f_low=500, target_f_high=4000,
            max_length=0.5, max_mouth_radius=0.2,
        )
        out = tmp_path / "candidates.csv"
        write_candidates_csv(candidates, str(out))
        assert out.exists()

        with open(out) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == len(candidates)
        assert set(rows[0].keys()) == {
            "candidate_id", "profile", "throat_radius", "mouth_radius", "length",
        }
