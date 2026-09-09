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
    origin=tmp_path/'origin.json'
    origin.write_text(json.dumps({'containers':{'horn-solver':{'id':'sha256:'+'0'*64}}}))
    monkeypatch.setattr(v, 'origin_files', lambda *args: {'origin_manifest':origin})
    monkeypatch.setattr(v, 'inspect_image', lambda value: {'Id':value})
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


@pytest.mark.parametrize('modal',[False,True])
def test_archived_workflow_and_resolution_evidence_have_recorded_identity(modal):
    import hashlib
    directory=v.ROOT/'data/validation'
    prefix='modal_' if modal else ''
    workflow='worked_example_modal_800_1600_manifest.json' if modal else 'worked_example_800_1600_manifest.json'
    for name in (prefix+'candidate_resolution_manifest.json',workflow):
        manifest=json.loads((directory/name).read_text())
        for filename,digest in manifest['files'].items():
            assert hashlib.sha256((directory/filename).read_bytes()).hexdigest()==digest
    result=json.loads((directory/(prefix+'candidate_resolution_reference.json')).read_text())
    assert result['passed'] and len(result['health']) == 7
    assert len(result['comparisons']) == 6 and all(row['passed'] for row in result['comparisons'])
    assert result['physical_validation_status']=='experimental_prediction'
    assert result['comparisons'][0]['kind']=='frequency'
    assert result['comparisons'][0]['coarse']=='original_ranking'
    assert 'ripple_db' in result['comparisons'][0]['changes']
    import tarfile
    with tarfile.open(directory/(prefix+'candidate_resolution_artifacts.tar.gz')) as archive:
        root='modal-candidate-resolution/' if modal else 'candidate-resolution/'
        def read(name):return archive.extractfile(root+name).read()
        protocol=json.loads(read('protocol.json'))
        evidence=json.loads(read('solve-evidence.json'))
        execution=json.loads(read('host-execution.json'))
        identity=json.loads((directory/(prefix+'candidate_resolution_manifest.json')).read_text())
        assert protocol['source_revision']==identity['reproduction_source_commit']
        assert protocol['limits']==v.LIMITS and protocol['cases']==v.CASES
        assert result['protocol_sha256']==evidence['protocol_sha256']==execution['protocol_sha256']==hashlib.sha256(read('protocol.json')).hexdigest()
        assert result['comparison_grid']=='union_log_frequency' and result['historical_reanalysis']
        assert result['analysis_revision']==identity['analysis_revision']
        assert result['previous_comparison_sha256']==hashlib.sha256(read('comparison.json')).hexdigest()
        import io
        with tarfile.open(fileobj=io.BytesIO(read('analysis_source.tar.gz'))) as source:
            for name,digest in result['analysis_source'].items():
                assert hashlib.sha256(source.extractfile(name).read()).hexdigest()==digest
        assert result['solve_evidence_sha256']==execution['solve_evidence_sha256']==hashlib.sha256(read('solve-evidence.json')).hexdigest()
        assert result['host_execution_sha256']==hashlib.sha256(read('host-execution.json')).hexdigest()
        assert execution['exit_code']==0
        assert protocol['solver_image_id']==evidence['solver_image_id']==execution['image_id']==result['solver_image_id']
        assert execution['image_id'] in execution['command']
        assert evidence['runtime']['petsc_scalar']=='complex128'
        assert all(evidence['runtime'][key] for key in ('python','numpy','scipy','dolfinx','gmsh','mpi_library','petsc','numpy_configuration'))
        for name,digest in {**protocol['inputs'],**evidence['files']}.items():
            assert hashlib.sha256(read(name)).hexdigest()==digest


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
    candidate=dict(driver_id='test',horn_label='example',drive_voltage_rms=2.83,observation_distance_m=1.,radiation_model='flanged_piston')
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


def test_modal_health_contract_is_not_conflated_with_local_radiation():
    good=frame();good['radiation_model']='modal_baffled'
    assert v.check_frame(good,[1.,2.],2,radiation_model='modal_baffled')['mesh_cells']==100
    with pytest.raises(ValueError,match='radiation_model'):v.check_frame(good,[1.,2.],2)


def test_refinement_uses_logarithmic_frequency_interpolation():
    coarse=dict(frequency=np.array([1.,4.]),spl=np.array([0.,2.]),z=np.array([1.+0j,3.+0j]))
    fine=dict(frequency=np.array([1.,2.,4.]),spl=np.array([0.,1.,2.]),z=np.array([1.+0j,2.+0j,3.+0j]))
    assert all(value==pytest.approx(0.) for value in v.curve_change(coarse,fine).values())


def test_non_nested_original_grid_peak_survives_in_both_directions():
    original=dict(frequency=np.array([1.,1.5,3.]),spl=np.array([0.,5.,0.]),z=np.array([1.+0j,4j,1.+0j]))
    baseline=dict(frequency=np.array([1.,2.,3.]),spl=np.zeros(3),z=np.ones(3,dtype=complex))
    for a,b in ((original,baseline),(baseline,original)):
        change=v.curve_change(a,b)
        assert change['spl_db']==5.
        assert change['impedance_db']==pytest.approx(20*np.log10(4))
        assert change['phase_deg']==90.


def test_historical_reanalysis_allows_only_a_preserved_comparator_change(tmp_path,monkeypatch):
    import hashlib,tarfile,io
    harness='scripts/validate_candidate_resolution.py';package='packages/horn-core/src/example.py'
    old=b'original comparator';old_digest=hashlib.sha256(old).hexdigest()
    monkeypatch.setattr(v,'source_identity',lambda:{harness:'new comparator',package:'same model'})
    protocol=dict(source={harness:old_digest,package:'same model'},cases=v.CASES,limits=v.LIMITS,
                  inputs={},solver_image_id='image',candidate_index=0,candidate={})
    (tmp_path/'protocol.json').write_text(json.dumps(protocol))
    (tmp_path/'origin_manifest').write_text(json.dumps({'containers':{'horn-solver':{'id':'image'}}}))
    (tmp_path/'ranking.json').write_text('[{}]')
    with tarfile.open(tmp_path/'study_source.tar.gz','w:gz') as archive:
        member=tarfile.TarInfo(harness);member.size=len(old);archive.addfile(member,io.BytesIO(old))
    with pytest.raises(ValueError,match='Source or fixed protocol'):v.verify(tmp_path)
    assert v.verify(tmp_path,analysis_only=True)==protocol
    monkeypatch.setattr(v,'source_identity',lambda:{harness:'new comparator',package:'changed model'})
    with pytest.raises(ValueError,match='new solve is required'):v.verify(tmp_path,analysis_only=True)
    monkeypatch.setattr(v,'source_identity',lambda:{harness:'new comparator',package:'same model'})
    protocol['source'][harness]='unpreserved old comparator'
    (tmp_path/'protocol.json').write_text(json.dumps(protocol))
    with pytest.raises(ValueError,match='not preserved'):v.verify(tmp_path,analysis_only=True)


@pytest.mark.parametrize('image',['horn-solver:latest','--format','sha256:bad'])
def test_mutable_or_invalid_image_identifiers_are_rejected_before_docker(image):
    with pytest.raises(ValueError,match='immutable'):
        v.inspect_image(image)


def test_host_executes_frozen_image_and_preserves_failure_without_a_success_seal(tmp_path,monkeypatch):
    from types import SimpleNamespace
    image_id='sha256:'+'1'*64
    (tmp_path/'protocol.json').write_text('{}')
    monkeypatch.setattr(v,'clean_source_revision',lambda:'source')
    monkeypatch.setattr(v,'verify',lambda out:{'solver_image_id':image_id})
    monkeypatch.setattr(v,'inspect_image',lambda image:{'Id':image,'Architecture':'amd64','Os':'linux'})
    commands=[]
    def execute(command,**kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=7 if command[:2]==['docker','run'] else 0)
    monkeypatch.setattr(v.subprocess,'run',execute)
    with pytest.raises(RuntimeError,match='exited 7'):
        v.solve(tmp_path,jobs=3)
    assert image_id in commands[0] and not any('latest' in part for part in commands[0])
    assert commands[-1][:3]==['docker','rm','-f']
    evidence=json.loads((tmp_path/'host-execution.json').read_text())
    assert evidence['image_id']==image_id and evidence['exit_code']==7
    assert evidence['solve_evidence_sha256'] is None
    claim=json.loads((tmp_path/'execution-claim.json').read_text())
    assert commands[0][commands[0].index('--name')+1]=='horn-resolution-'+claim['invocation']


def test_overlapping_execution_cannot_launch_or_remove_another_container(tmp_path,monkeypatch):
    image_id='sha256:'+'1'*64
    (tmp_path/'protocol.json').write_text('{}')
    (tmp_path/'execution-claim.json').write_text('{"invocation":"already-running"}')
    monkeypatch.setattr(v,'clean_source_revision',lambda:'source')
    monkeypatch.setattr(v,'verify',lambda out:{'solver_image_id':image_id})
    monkeypatch.setattr(v,'inspect_image',lambda image:{'Id':image,'Architecture':'amd64','Os':'linux'})
    monkeypatch.setattr(v.subprocess,'run',lambda *a,**kw:pytest.fail('Overlapping invocation must not touch Docker'))
    with pytest.raises(FileExistsError):v.solve(tmp_path)
    assert not (tmp_path/'host-execution.json').exists()
    assert json.loads((tmp_path/'execution-claim.json').read_text())['invocation']=='already-running'
