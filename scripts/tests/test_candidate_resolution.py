"""Regression gates for candidate resolution, including unsampled peaks."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

PATH=Path(__file__).resolve().parents[1]/'validate_candidate_resolution.py'
spec=importlib.util.spec_from_file_location('candidate_resolution',PATH)
v=importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def test_extra_fine_grid_peak_cannot_disappear_through_downsampling():
    coarse=dict(frequency=np.array([1.,3.]),spl=np.zeros(2),z=np.ones(2,dtype=complex))
    fine=dict(frequency=np.array([1.,2.,3.]),spl=np.array([0.,4.,0.]),z=np.array([1.,2j,1.]))
    change=v.curve_change(coarse,fine)
    assert change['spl_db']==4.
    assert change['phase_deg']==90.
    assert change['impedance_db']==pytest.approx(20*np.log10(2))
    assert v.ripple(fine,[1.5,2.5])==pytest.approx(4*np.log(2.5/2)/np.log(3/2))


def test_partial_band_and_zero_impedance_cannot_pass():
    full=dict(frequency=np.array([1.,3.]),spl=np.zeros(2),z=np.ones(2,dtype=complex))
    partial=dict(full,frequency=np.array([2.,3.]))
    with pytest.raises(ValueError,match='same band'):
        v.curve_change(partial,full)
    with pytest.raises(ValueError,match='not covered'):
        v.ripple(partial,[1.,3.])
    with pytest.raises(ValueError,match='Zero impedance'):
        v.curve_change(full,dict(full,z=np.zeros(2)))


def frame():
    return pd.DataFrame(dict(frequency=[1.,2.],spl=[0.,0.],z_real=[1.,1.],
        schema_version=2,bc_mode='dirichlet',phasor_convention='exp(+iwt)_rms',
        radiation_model='flanged_piston',loss_model='lossless',element_degree=1,
        input_acoustic_power_w=1.,mouth_acoustic_power_w=1.,viscous_wall_power_w=0.,
        thermal_wall_power_w=0.,relative_residual=1e-12,converged_reason=4,mesh_cells=100))


@pytest.mark.parametrize('column,value', [('frequency',1.),('spl',np.nan),
    ('relative_residual',-1.),('relative_residual',1e-3),('converged_reason',0),
    ('mouth_acoustic_power_w',.8),('bc_mode','velocity'),('element_degree',2)])
def test_failed_health_or_changed_contract_cannot_pass(column,value):
    good=frame()
    assert v.check_frame(good,[1.,2.],2)['mesh_cells']==100
    good[column]=value
    with pytest.raises(ValueError):
        v.check_frame(good,[1.,2.],2)


def test_preparation_freezes_candidate_driver_source_and_limits(tmp_path, monkeypatch):
    monkeypatch.setattr(v, 'clean_source_revision', lambda: 'test-revision')
    monkeypatch.setattr(v, 'origin_files', lambda *args: {})
    candidate=dict(driver_id='test',loss_model='lossless',radiation_model='flanged_piston',
        element_degree=1,throat_radius=.01,mouth_radius=.03,length=.1,
        drive_voltage_rms=2.83,observation_distance_m=1.,profile='conical')
    ranking=tmp_path/'input-ranking.json'
    driver=tmp_path/'input-driver.json'
    ranking.write_text(json.dumps([candidate]))
    driver.write_text(json.dumps(dict(driver_id='wrong')))
    out=tmp_path/'study'
    with pytest.raises(ValueError,match='does not match'):
        v.prepare(out,ranking,driver,800.,1600.,0)
    assert not out.exists()
    driver.write_text(json.dumps(dict(driver_id='test')))
    v.prepare(out,ranking,driver,800.,1600.,0)
    assert v.verify(out)['candidate']==candidate
    driver_copy=out/'driver.json'
    driver_copy.write_text('{}')
    with pytest.raises(ValueError,match='Frozen input'):
        v.verify(out)


def test_archived_workflow_and_resolution_evidence_have_recorded_identity():
    import hashlib
    directory=v.ROOT/'data/validation'
    for name in ('candidate_resolution_manifest.json','worked_example_800_1600_manifest.json'):
        manifest=json.loads((directory/name).read_text())
        for filename,digest in manifest['files'].items():
            assert hashlib.sha256((directory/filename).read_bytes()).hexdigest()==digest
    result=json.loads((directory/'candidate_resolution_reference.json').read_text())
    assert result['passed'] and len(result['health']) in (6,7)
    assert len(result['comparisons']) in (5,6) and all(row['passed'] for row in result['comparisons'])
    assert result['physical_validation_status']=='experimental_prediction'


def test_wavelength_cap_cannot_collapse_refinement():
    assert v.check_mesh_schedule(1600.) == [.01, .006, .004]
    for high in (7000., 12000.):
        with pytest.raises(ValueError, match='collapses'):
            v.check_mesh_schedule(high)


def test_ripple_uses_production_logarithmic_band_edges():
    curve = dict(frequency=np.array([100., 1000., 10000.]), spl=np.array([0., 10., 0.]))
    assert v.ripple(curve, [np.sqrt(100*1000), np.sqrt(1000*10000)]) == pytest.approx(5.)


@pytest.mark.parametrize('change', ['modified', 'staged', 'untracked', 'ignored'])
def test_dirty_source_cannot_be_advertised_as_a_reproducible_commit(tmp_path, monkeypatch, change):
    import subprocess
    def git(*args):
        subprocess.run(['git', '-C', str(tmp_path), *args], check=True, capture_output=True)
    git('init')
    git('config', 'user.name', 'Test')
    git('config', 'user.email', 'test@example.invalid')
    source = tmp_path/'packages/horn-core/src/horn_core'
    source.mkdir(parents=True)
    module = source/'example.py'
    module.write_text('original = True')
    (tmp_path/'.gitignore').write_text('ignored.py\n')
    git('add', '.')
    git('commit', '-m', 'Initial')
    monkeypatch.setattr(v, 'ROOT', tmp_path)
    monkeypatch.setattr(v, 'source_identity', lambda: {str(p.relative_to(tmp_path)): v.sha(p) for p in source.glob('*.py')})
    assert len(v.clean_source_revision()) == 40
    if change in ('modified', 'staged'):
        module.write_text('original = False')
        if change == 'staged':
            git('add', '.')
    else:
        (source/f'{change}.py').write_text('extra = True')
    with pytest.raises(ValueError, match='clean|Ignored'):
        v.clean_source_revision()


def origin_fixture(tmp_path):
    import tarfile,io,hashlib
    run=tmp_path/'run';report=run/'outputs/auto/report';refine=run/'outputs/auto/refinement'
    report.mkdir(parents=True);refine.mkdir()
    candidate=dict(driver_id='test',horn_label='example',drive_voltage_rms=2.83,observation_distance_m=1.)
    ranking=report/'auto_ranking.json';ranking.write_text(json.dumps([candidate]))
    driver=tmp_path/'driver.json';driver.write_text(json.dumps(dict(driver_id='test',parameters={'re_ohm':6.})))
    raw=driver.read_bytes();digest=hashlib.sha256(raw).hexdigest();name='data/drivers/test.json'
    manifest=dict(status='completed',exit_code=0,input_sha256={'--drivers_db':{name:digest}},source_sha256={name:digest,**{k:h for k,h in v.source_identity().items() if k.startswith("packages/")}})
    (run/'manifest.json').write_text(json.dumps(manifest))
    parameters=dict(mesh_size=.01,num_sections=20,num_intervals=101,target_f_low=800.,target_f_high=1600.,element_degree=1,radiation_model='flanged_piston',loss_model='lossless',voltage_rms=2.83,observation_distance=1.)
    (run/'outputs/resolved_specification.json').write_text(json.dumps({'parameters':parameters}))
    with tarfile.open(run/'source.tar.gz','w:gz') as archive:
        member=tarfile.TarInfo(name);member.size=len(raw);archive.addfile(member,io.BytesIO(raw))
    (refine/'example.step').write_text('original geometry')
    (refine/'example_results.csv').write_text('original response')
    manifest['output_sha256']={str(p.relative_to(run)):v.sha(p) for p in (run/'outputs').rglob('*') if p.is_file()}
    (run/'manifest.json').write_text(json.dumps(manifest))
    return run,ranking,driver,candidate


@pytest.mark.parametrize('field,value',[('num_sections',2),('mesh_size',.02),('num_intervals',20)])
def test_originating_coarse_settings_cannot_be_silently_replaced(tmp_path,field,value):
    run,ranking,driver,candidate=origin_fixture(tmp_path)
    assert len(v.origin_files(run,ranking,driver,candidate,800.,1600.))==5
    path=run/'outputs/resolved_specification.json';p=json.loads(path.read_text())
    p['parameters'][field]=value;path.write_text(json.dumps(p))
    with pytest.raises(ValueError,match='resolution/specification'):
        v.origin_files(run,ranking,driver,candidate,800.,1600.)


def test_same_driver_id_with_refreshed_parameters_is_not_the_ranked_driver(tmp_path):
    run,ranking,driver,candidate=origin_fixture(tmp_path)
    driver.write_text(json.dumps(dict(driver_id='test',parameters={'re_ohm':8.})))
    with pytest.raises(ValueError,match='Driver bytes'):
        v.origin_files(run,ranking,driver,candidate,800.,1600.)


def test_spatial_cases_use_the_fine_grid_and_preserve_original_geometry():
    assert all(case['points']==201 for name,case in v.CASES.items() if name!='frequency_101')
    assert v.CASES['frequency_101']['sections']==v.CASES['loft_80']['sections']


@pytest.mark.parametrize('name', ['outputs/auto/refinement/example.step','outputs/auto/refinement/example_results.csv','outputs/auto/report/auto_ranking.json'])
def test_changed_completed_output_cannot_enter_study(tmp_path,name):
    run,ranking,driver,candidate=origin_fixture(tmp_path)
    path=run/name;path.write_text(path.read_text()+' ')
    with pytest.raises(ValueError,match='completion-time digest'):
        v.origin_files(run,ranking,driver,candidate,800.,1600.)


def test_changed_package_source_cannot_be_conflated_with_resolution(tmp_path):
    run,ranking,driver,candidate=origin_fixture(tmp_path)
    path=run/'manifest.json';manifest=json.loads(path.read_text())
    key=next(k for k in manifest['source_sha256'] if k.startswith('packages/'))
    manifest['source_sha256'][key]='different';path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='package source differs'):
        v.origin_files(run,ranking,driver,candidate,800.,1600.)


def test_originating_band_grid_can_differ_from_global_geometric_grid():
    grid=np.unique(np.r_[np.geomspace(1.,2.5,3),np.geomspace(2.5,4.,3)])
    good=pd.concat([frame().iloc[:1]]*len(grid),ignore_index=True)
    good['frequency']=grid
    assert v.check_frame(good,[1.,4.],len(grid),expected=grid)['mesh_cells']==100
    with pytest.raises(ValueError,match='frequency grid'):
        v.check_frame(good,[1.,4.],len(grid))
