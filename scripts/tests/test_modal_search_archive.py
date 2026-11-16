"""Reassess the modal shortlist against every sealed FEM curve in the grid."""
import hashlib
import json
from pathlib import Path
import tarfile

from horn_core.candidates import CandidateGeometry
from horn_core.parameters import DriverParameters
from horn_analysis.lem_prescreen import lem_prescreen_candidates
from horn_analysis.rank_pipeline import rank_horn_drivers
from horn_analysis.scoring import TargetSpec
from horn_analysis.search import screening_audit

DATA=Path(__file__).resolve().parents[2]/'data/validation'


def test_modal_exhaustive_search_replays_from_frozen_responses(tmp_path):
    identity=json.loads((DATA/'modal_search_manifest.json').read_text())
    for name,digest in identity['files'].items():
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest()==digest
    with tarfile.open(DATA/'modal_search_artifacts.tar.gz') as archive:
        def read(name):return archive.extractfile('modal-search/'+name).read()
        host=json.loads(read('host-execution.json'))
        inputs=json.loads(read('execution-input.json'))
        result=json.loads(read('search_benchmark.json'))
        runtime=json.loads(read('runtime.json'))
        assert host['exit_code']==0
        assert host['image_id']==inputs['image_id']==identity['solver_image_id']
        assert inputs['source_revision']==identity['reproduction_source_commit']
        assert host['image_id'] in host['command']
        assert runtime['petsc_scalar']=='complex128'
        for name,digest in host['files'].items():
            assert hashlib.sha256(read(name)).hexdigest()==digest
        assert result['passed'] and result['physical_validation_passed'] is False
        assert result['radiation_model']=='modal_baffled'
        assert result['geometry_count']==20 and result['driver_count']==3
        assert len(result['cases'])==4
        assert len(result['raw_sha256'])==40
        # The frozen shortlist and audit figures were recorded with the ideal
        # area transformer (no front chamber); replay under that same model.
        drivers=[DriverParameters(**{**record,'front_chamber_volume_m3':0.}) for record in result['drivers']]
        for name,digest in result['raw_sha256'].items():
            raw=read(name);assert hashlib.sha256(raw).hexdigest()==digest
            if name.endswith('.csv'):(tmp_path/name).write_bytes(raw)
        candidates=[CandidateGeometry(**record) for record in result['candidates']]
        for case in result['cases']:
            target=TargetSpec(**case['target'])
            # The frozen shortlist was chosen with the score of its day. Re-screen
            # with the current score (analytical, no solves) so the audit compares
            # one scoring rule against the frozen FEM responses.
            screen=lem_prescreen_candidates(candidates,drivers,target.f_low_hz,target.f_high_hz,
                tuple(result['sim_band']),num_frequencies=result['frequencies'],
                top_n=result['shortlist_budget'],target=target,radiation_model=result['radiation_model'])
            rows=[]
            for candidate in result['candidates']:
                rows.extend(rank_horn_drivers(str(tmp_path/(candidate['candidate_id']+'_results.csv')),
                    candidate['candidate_id'],candidate['throat_radius'],drivers,
                    target,top_n=3))
            assert len(rows)==60
            assert all(row['output_metric']=='modal_baffled_on_axis' for row in rows)
            audit=screening_audit(rows,screen['filtered_candidate_ids'])
            assert audit['passed'] and audit['winner_retained']
            assert audit['feasible_top_k_recall']==case['audit']['feasible_top_k_recall']==1.
            assert audit['score_regret']==case['audit']['score_regret']==0.
            assert sum(row['model_feasible'] for row in rows)==case['audit']['feasible_pair_count']>=10
