#!/usr/bin/env python3
"""Freeze, solve and check a generated candidate's numerical resolution.

This study does not certify a physical assembly or improve driver evidence.
Run prepare and solve on the host; solve launches the frozen solver container.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import re
import os
import uuid
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    'mesh_10': dict(mesh_size=.010, sections=20, points=201),
    'mesh_6': dict(mesh_size=.006, sections=20, points=201),
    'mesh_4': dict(mesh_size=.004, sections=20, points=201),
    'loft_40': dict(mesh_size=.004, sections=40, points=201),
    'loft_80': dict(mesh_size=.004, sections=80, points=201),
    'frequency_101': dict(mesh_size=.004, sections=80, points=101),
}
LIMITS = dict(mesh_spl_db=.5, mesh_impedance_db=.5, mesh_phase_deg=5.,
              loft_spl_db=.5, loft_impedance_db=.5, loft_phase_deg=5.,
              frequency_spl_db=.2, frequency_impedance_db=.2,
              frequency_phase_deg=2., frequency_ripple_db=.2)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def inspect_image(image_id):
    if not re.fullmatch(r'sha256:[0-9a-f]{64}',image_id):
        raise ValueError('An immutable solver image ID is required')
    image=json.loads(subprocess.check_output(['docker','image','inspect',image_id],text=True))[0]
    if image['Id']!=image_id:raise ValueError('Solver image identity changed')
    return image


def source_identity():
    paths = list((ROOT/'packages').glob('*/src/**/*.py')) + [Path(__file__)]
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}


def clean_source_revision():
    """Host-side provenance check; container solves verify the frozen bytes."""
    def git(*args):
        return subprocess.check_output(['git', '-C', str(ROOT), *args], text=True)
    if git('status', '--porcelain', '--untracked-files=all'):
        raise ValueError('A clean tracked checkout with no untracked files is required')
    tracked = set(git('ls-files').splitlines())
    if set(source_identity()) - tracked:
        raise ValueError('Ignored or untracked Python source cannot identify a revision')
    return git('rev-parse', 'HEAD').strip()


def check_mesh_schedule(high):
    # The production solver uses c/(6*f_max); preserve all three refinements.
    cap = 343.0 / (6.0 * high * np.sqrt(2))
    effective = [min(CASES[name]['mesh_size'], cap) for name in ('mesh_10', 'mesh_6', 'mesh_4')]
    if not effective[0] > effective[1] > effective[2]:
        raise ValueError('Wavelength cap collapses the mesh refinement schedule')
    return effective


def bind_source():
    roots = {p.name.replace('-', '_'): p/'src' for p in (ROOT/'packages').glob('horn-*')}
    for name, module in tuple(sys.modules.items()):
        package = name.split('.')[0]
        if package in roots:
            path = getattr(module, '__file__', None)
            if path is None or not Path(path).resolve().is_relative_to(roots[package].resolve()):
                raise ValueError(f'{name} was imported outside the frozen checkout')
    paths = [str(p) for p in roots.values()]
    sys.path[:] = paths + [p for p in sys.path if p not in paths]


def origin_files(run_dir, ranking, driver, candidate, low, high):
    """Bind the experiment to the actual completed ranking and original inputs."""
    if run_dir is None:
        raise ValueError('The originating completed run directory is required')
    run_dir=Path(run_dir)
    manifest_path=run_dir/'manifest.json'
    resolved_path=run_dir/'outputs/resolved_specification.json'
    manifest=json.loads(manifest_path.read_text())
    parameters=json.loads(resolved_path.read_text())['parameters']
    if manifest.get('status')!='completed' or manifest.get('exit_code')!=0:
        raise ValueError('The originating run did not complete')
    # Documentation/harness revisions may differ; all executing package bytes
    # must be identical so this study varies resolution, not implementation.
    current={name:digest for name,digest in source_identity().items() if name.startswith('packages/')}
    original={name:digest for name,digest in manifest['source_sha256'].items()
              if name.startswith('packages/') and '/src/' in name and name.endswith('.py')}
    if original!=current:
        raise ValueError('Originating package source differs from the resolution study')
    expected=dict(mesh_size=.01,num_sections=20,num_intervals=101,
                  target_f_low=low,target_f_high=high,element_degree=1,
                  radiation_model=candidate.get('radiation_model','flanged_piston'),loss_model='lossless',
                  voltage_rms=candidate['drive_voltage_rms'],
                  observation_distance=candidate['observation_distance_m'])
    if any(parameters.get(key)!=value for key,value in expected.items()):
        raise ValueError('Originating resolution/specification differs from this fixed protocol')
    if sha(ranking)!=sha(run_dir/'outputs/auto/report/auto_ranking.json'):
        raise ValueError('Ranking is not from the originating run')
    label=candidate['horn_label']
    if not re.fullmatch(r'[A-Za-z0-9_-]+',label):
        raise ValueError('Invalid original horn label')
    refined=run_dir/'outputs/auto/refinement'
    if (refined/f'{label}.step').exists():
        step=refined/f'{label}.step';response=refined/f'{label}_results.csv'
    else:
        step=run_dir/f'outputs/auto/geometry/horn_{label}.step'
        response=run_dir/f'outputs/auto/{label}_results.csv'
    if not step.is_file() or not response.is_file():
        raise ValueError('Original ranked STEP/response is required')
    for path in (run_dir/'outputs/auto/report/auto_ranking.json',resolved_path,step,response):
        if manifest.get('output_sha256',{}).get(str(path.relative_to(run_dir)))!=sha(path):
            raise ValueError('Originating output has no matching completion-time digest')
    source=run_dir/'source.tar.gz'
    expected_hashes={digest for name,digest in manifest['input_sha256']['--drivers_db'].items() if name.endswith('.json')}
    records={};found=set()
    with tarfile.open(source,'r:gz') as archive:
        for member in archive.getmembers():
            digest=manifest['source_sha256'].get(member.name)
            if not member.isfile() or not member.name.endswith('.json') or digest not in expected_hashes:
                continue
            raw=archive.extractfile(member).read()
            if hashlib.sha256(raw).hexdigest()!=digest:
                raise ValueError('Originating driver snapshot changed')
            record=json.loads(raw)
            identifier=record.get('driver_id') if isinstance(record,dict) else None
            if not identifier:
                raise ValueError('This study requires individually archived driver records')
            records.setdefault(identifier,set()).add(digest);found.add(digest)
    if found!=expected_hashes or not expected_hashes:
        raise ValueError('Original driver database is not completely source-snapshotted')
    if any(len(digests)!=1 for digests in records.values()):
        raise ValueError('Ambiguous duplicate driver IDs in the originating database')
    if records.get(candidate['driver_id'])!={sha(driver)}:
        raise ValueError('Driver bytes differ from the originating ranking')
    return dict(origin_manifest=manifest_path,origin_specification=resolved_path,
                origin_source=source,original_step=step,original_response=response)


def prepare(out, ranking, driver, low, high, index, run_dir=None):
    if not np.isfinite([low, high]).all() or not 0 < low < high:
        raise ValueError('A positive ordered target band is required')
    candidates = json.loads(ranking.read_text())
    if isinstance(candidates, dict):
        candidates = candidates['results']
    if index < 0:
        raise ValueError('Candidate index must be nonnegative')
    candidate = candidates[index]
    record = json.loads(driver.read_text())
    if record['driver_id'] != candidate['driver_id']:
        raise ValueError('Driver record does not match the selected candidate')
    if (candidate.get('loss_model') != 'lossless'
            or candidate.get('radiation_model') not in {'flanged_piston','modal_baffled'}
            or candidate.get('element_degree') != 1):
        raise ValueError('This protocol supports lossless local-piston or modal-baffled P1 only')
    for key in ('throat_radius', 'mouth_radius', 'length', 'drive_voltage_rms', 'observation_distance_m'):
        if not np.isfinite(candidate[key]) or candidate[key] <= 0:
            raise ValueError(f'Invalid candidate {key}')
    origin=origin_files(run_dir,ranking,driver,candidate,low,high)
    check_mesh_schedule(high)
    revision = clean_source_revision()
    solver_image=json.loads(origin['origin_manifest'].read_text())['containers']['horn-solver']['id']
    inspected=inspect_image(solver_image)
    if inspected['Id']!=solver_image:
        raise ValueError('Originating solver image is unavailable')
    out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(ranking, out/'ranking.json')
    shutil.copyfile(driver, out/'driver.json')
    for name,path in origin.items():
        shutil.copyfile(path,out/name)
    write_json(out/'protocol.json', dict(schema_version=1,
        prepared_at=datetime.now(timezone.utc).isoformat(),
        source_revision=revision,
        solver_image_id=solver_image,
        source=source_identity(), inputs={p:sha(out/p) for p in ('ranking.json','driver.json',*origin)},
        candidate_index=index, candidate=candidate, target_band_hz=[low,high],
        simulation_band_hz=[low/np.sqrt(2),high*np.sqrt(2)], cases=CASES, limits=LIMITS,
        scope='Numerical resolution of this fixed ideal candidate; no physical qualification'))


def verify(out, *, analysis_only=False):
    p = json.loads((out/'protocol.json').read_text())
    current=source_identity()
    if analysis_only:
        # Reassess sealed solves with a new comparator, never a different model.
        harness='scripts/validate_candidate_resolution.py'
        if set(p['source'])!=set(current) or any(p['source'][k]!=current[k] for k in current if k!=harness):
            raise ValueError('Executed package source differs; a new solve is required')
        with tarfile.open(out/'study_source.tar.gz') as archive:
            if hashlib.sha256(archive.extractfile(harness).read()).hexdigest()!=p['source'][harness]:
                raise ValueError('Original study harness source is not preserved')
    elif p['source'] != current:
        raise ValueError('Source or fixed protocol changed; prepare a fresh study')
    if p['cases'] != CASES or p['limits'] != LIMITS:
        raise ValueError('Source or fixed protocol changed; prepare a fresh study')
    if any(sha(out/name) != digest for name,digest in p['inputs'].items()):
        raise ValueError('Frozen input changed')
    origin=json.loads((out/'origin_manifest').read_text())
    if p['solver_image_id']!=origin['containers']['horn-solver']['id']:
        raise ValueError('Solver image differs from the originating execution')
    ranking=json.loads((out/'ranking.json').read_text())
    if isinstance(ranking,dict):
        ranking=ranking['results']
    if ranking[p['candidate_index']] != p['candidate']:
        raise ValueError('Frozen candidate changed')
    return p


def solve_case(arguments):
    out,name=arguments
    bind_source()
    from horn_geometry.generator import create_horn
    from horn_solver.solver import run_simulation_from_step
    p=verify(out);c=p['candidate'];case=CASES[name]
    runtime=runtime_identity()
    if runtime['petsc_scalar']!='complex128':raise ValueError('Complex double PETSc is required')
    directory=out/name
    directory.mkdir(exist_ok=False)
    step=directory/'horn.step'
    if case['sections']==20:
        shutil.copyfile(out/'original_step',step)
    else:
        create_horn(c['profile'],c['throat_radius'],c['mouth_radius'],c['length'],step,
                    num_sections=case['sections'])
    run_simulation_from_step(str(step),tuple(p['simulation_band_hz']),case['points'],
        {'length':c['length']},str(directory/'response.csv'),max(p['simulation_band_hz']),
        mesh_size=case['mesh_size'],element_degree=1,radiation_model=c['radiation_model'],loss_model='lossless')
    verify(out)
    print('Completed resolution case:',name,flush=True)
    return runtime


def runtime_identity():
    import platform, io, contextlib
    import scipy, pandas, dolfinx, gmsh, mpi4py
    from mpi4py import MPI
    from petsc4py import PETSc
    stream=io.StringIO()
    with contextlib.redirect_stdout(stream):np.show_config()
    return dict(python=sys.version,architecture=platform.machine(),numpy=np.__version__,
                scipy=scipy.__version__,pandas=pandas.__version__,dolfinx=dolfinx.__version__,
                gmsh=gmsh.__version__,mpi4py=mpi4py.__version__,mpi_library=MPI.Get_library_version(),
                petsc=list(PETSc.Sys.getVersion()),petsc_scalar=str(np.dtype(PETSc.ScalarType)),
                numpy_configuration=stream.getvalue(),
                openblas_threads=os.environ.get('OPENBLAS_NUM_THREADS'),omp_threads=os.environ.get('OMP_NUM_THREADS'))


def solve(out,jobs=1):
    """Host-owned execution pins the actual image, not a caller-supplied label."""
    clean_source_revision();p=verify(out)
    image=inspect_image(p['solver_image_id'])
    if image['Id']!=p['solver_image_id']:raise ValueError('Solver image identity changed')
    if (out/'host-execution.json').exists() or (out/'solve-evidence.json').exists():
        raise FileExistsError('Study execution already exists')
    # Claim the output directory atomically before launching anything. Keep the
    # claim on failure: a failed study is evidence, not a directory to overwrite.
    token=uuid.uuid4().hex
    with (out/'execution-claim.json').open('x') as claim:
        json.dump(dict(invocation=token,protocol_sha256=sha(out/'protocol.json')),claim)
    name='horn-resolution-'+token
    command=['docker','run','--rm','--name',name,
             '-e','OPENBLAS_NUM_THREADS=1','-e','OMP_NUM_THREADS=1','-e','PYTHONPATH=/usr/local/lib',
             '-e','HORN_STUDY_IMAGE_ID='+image['Id'],
             '-v',f'{ROOT}:/workspace:ro','-v',f'{out}:/study','-w','/workspace',image['Id'],
             'python3','scripts/validate_candidate_resolution.py','_solve','/study','--jobs',str(jobs)]
    code=None
    try:
        code=subprocess.run(command,timeout=7200).returncode
        if code!=0:raise RuntimeError(f'Solver container exited {code}')
        verify(out);clean_source_revision()
        evidence=json.loads((out/'solve-evidence.json').read_text())
        if evidence['solver_image_id']!=image['Id']:raise ValueError('Solver evidence image differs')
    finally:
        subprocess.run(['docker','rm','-f',name],capture_output=True)
        write_json(out/'host-execution.json',dict(command=command,image_id=image['Id'],
            image_architecture=image['Architecture'],image_os=image['Os'],exit_code=code,
            protocol_sha256=sha(out/'protocol.json'),
            solve_evidence_sha256=sha(out/'solve-evidence.json') if code==0 and (out/'solve-evidence.json').exists() else None))


def solve_in_container(out, jobs=1):
    if not isinstance(jobs,int) or not 1<=jobs<=len(CASES):
        raise ValueError('Jobs must be between one and the number of cases')
    p=verify(out)
    if os.environ.get('HORN_STUDY_IMAGE_ID')!=p['solver_image_id']:
        raise ValueError('Use the host solve command to enforce the frozen image')
    arguments=[(out,name) for name in CASES]
    if jobs==1:
        runtimes=[solve_case(argument) for argument in arguments]
    else:
        # Fresh processes isolate Gmsh, MPI/PETSc and source imports per worker.
        from multiprocessing import get_context
        from concurrent.futures import ProcessPoolExecutor
        # Pool.__exit__ terminates even successful workers. PETSc intercepts
        # SIGTERM and can hang in MPI_Abort; orderly executor shutdown lets
        # each completed worker finalize MPI normally. Crashes break the pool.
        with ProcessPoolExecutor(max_workers=jobs,mp_context=get_context('spawn')) as pool:
            runtimes=list(pool.map(solve_case,arguments))
    runtime=runtimes[0]
    if any(item!=runtime for item in runtimes):raise ValueError('Runtime changed between solver cases')
    verify(out)
    files=list(out.glob('*/horn.step'))+list(out.glob('*/response.csv'))
    write_json(out/'solve-evidence.json',dict(protocol_sha256=sha(out/'protocol.json'),jobs=jobs,
        solver_image_id=p['solver_image_id'],runtime=runtime,
        files={str(f.relative_to(out)):sha(f) for f in sorted(files)}))


def check_frame(frame, band, count, *, expected=None, radiation_model='flanged_piston'):
    expected=np.geomspace(*band,count) if expected is None else np.asarray(expected)
    if (len(frame)!=count or not np.isfinite(frame.select_dtypes(include='number')).all().all()
            or not np.allclose(frame.frequency,expected,rtol=1e-12,atol=0)):
        raise ValueError('Incomplete, nonfinite or incorrect frequency grid')
    for key,value in dict(schema_version=2,bc_mode='dirichlet',phasor_convention='exp(+iwt)_rms',
                          radiation_model=radiation_model,loss_model='lossless',element_degree=1).items():
        if not (frame[key]==value).all():
            raise ValueError(f'Unexpected model contract: {key}')
    incoming=frame.input_acoustic_power_w.to_numpy()
    outgoing=(frame.mouth_acoustic_power_w+frame.viscous_wall_power_w+frame.thermal_wall_power_w).to_numpy()
    balance=np.abs(incoming-outgoing)/np.maximum(np.abs(incoming),1e-30)
    if ((frame.relative_residual<0).any() or (frame.relative_residual>1e-8).any()
            or (frame.converged_reason<=0).any() or (frame.z_real<0).any()
            or (incoming<0).any() or (outgoing<0).any() or np.max(balance)>1e-7):
        raise ValueError('Failed numerical health or power balance')
    return dict(mesh_cells=int(frame.mesh_cells.iloc[0]),max_residual=float(frame.relative_residual.max()),
                max_relative_power_balance_error=float(np.max(balance)))


def curve_change(coarse, fine):
    """Preserve every sample from both grids, including non-nested band grids."""
    if not np.allclose(coarse['frequency'][[0,-1]],fine['frequency'][[0,-1]],rtol=1e-12,atol=0):
        raise ValueError('Comparison grids do not cover the same band')
    f=np.unique(np.r_[coarse['frequency'],fine['frequency']])
    levels=[]; impedances=[]
    for curve in (coarse,fine):
        levels.append(np.interp(np.log(f),np.log(curve['frequency']),curve['spl']))
        impedances.append(np.interp(np.log(f),np.log(curve['frequency']),curve['z'].real)
                          +1j*np.interp(np.log(f),np.log(curve['frequency']),curve['z'].imag))
    z,other=impedances
    if np.any(np.abs(z)<=0) or np.any(np.abs(other)<=0):
        raise ValueError('Zero impedance cannot support a magnitude/phase comparison')
    return dict(spl_db=float(np.max(np.abs(levels[0]-levels[1]))),
        impedance_db=float(np.max(np.abs(20*np.log10(np.abs(z)/np.abs(other))))),
        phase_deg=float(np.max(np.abs(np.angle(z*np.conj(other),deg=True)))))


def ripple(curve, band):
    f=curve['frequency']
    if f[0]>band[0] or f[-1]<band[1]:
        raise ValueError('Target band not covered')
    points=np.r_[band[0],f[(f>band[0]) & (f<band[1])],band[1]]
    return float(np.ptp(np.interp(np.log(points),np.log(f),curve['spl'])))


def compare(out, *, reanalyze=False):
    analysis_revision=clean_source_revision()
    analysis_source=source_identity()
    bind_source()
    import pandas as pd
    from horn_drivers.loader import _driver_from_dict
    from horn_analysis.evaluation import coupled_output, evaluate_response
    from horn_analysis.scoring import TargetSpec
    p=verify(out,analysis_only=reanalyze)
    previous_comparison=sha(out/'comparison.json') if reanalyze else None
    evidence=json.loads((out/'solve-evidence.json').read_text())
    execution=json.loads((out/'host-execution.json').read_text())
    if (execution['exit_code']!=0 or execution['image_id']!=p['solver_image_id']
            or execution['protocol_sha256']!=sha(out/'protocol.json')
            or execution['solve_evidence_sha256']!=sha(out/'solve-evidence.json')
            or evidence.get('solver_image_id')!=p['solver_image_id']
            or evidence.get('runtime',{}).get('petsc_scalar')!='complex128'):
        raise ValueError('Missing or changed solver runtime/execution identity')
    expected={f'{name}/{file}' for name in CASES for file in ('horn.step','response.csv')}
    if (set(evidence['files'])!=expected or evidence['protocol_sha256']!=sha(out/'protocol.json')
            or any(sha(out/name)!=digest for name,digest in evidence['files'].items())):
        raise ValueError('Solve evidence changed or is incomplete')
    candidate=p['candidate']
    driver=_driver_from_dict(json.loads((out/'driver.json').read_text()))
    target=SimpleNamespace(voltage_rms=candidate['drive_voltage_rms'],
                           observation_distance_m=candidate['observation_distance_m'])
    curves,health={},{}
    original=pd.read_csv(out/'original_response')
    parameters=json.loads((out/'origin_specification').read_text())['parameters']
    bands=parameters['num_bands'];points=max(2,int(np.ceil(parameters['num_intervals']/bands)))
    low,high=p['simulation_band_hz'];width=(high-low)/bands
    original_grid=np.unique(np.concatenate([np.geomspace(low+i*width,low+(i+1)*width,points) for i in range(bands)]))
    health['original_ranking']=check_frame(original,p['simulation_band_hz'],len(original_grid),expected=original_grid,radiation_model=candidate['radiation_model'])
    original_levels,_,_=coupled_output(original,driver,target)
    metrics=evaluate_response(original.frequency,original_levels,TargetSpec(*p['target_band_hz']))
    for key in ('passband_ripple_db','avg_sensitivity_db'):
        if not np.isclose(metrics[key],candidate[key],atol=1e-8,rtol=0):
            raise ValueError('Original response/driver does not reproduce ranking metrics')
    curves['original_ranking']=dict(frequency=original.frequency.to_numpy(),spl=original_levels,z=(original.z_real+1j*original.z_imag).to_numpy())
    for name,case in CASES.items():
        frame=pd.read_csv(out/name/'response.csv')
        health[name]=check_frame(frame,p['simulation_band_hz'],case['points'],radiation_model=candidate['radiation_model'])
        levels,_,_=coupled_output(frame,driver,target)
        curves[name]=dict(frequency=frame.frequency.to_numpy(),spl=levels,
                         z=frame.z_real.to_numpy()+1j*frame.z_imag.to_numpy())
    check_mesh_schedule(p['target_band_hz'][1])
    counts = [health[name]['mesh_cells'] for name in ('mesh_10', 'mesh_6', 'mesh_4')]
    if not counts[0] < counts[1] < counts[2]:
        raise ValueError('Mesh cases did not produce strictly increasing cell counts')
    results=[]
    for kind,a,b in [('frequency','original_ranking','mesh_10'),('mesh','mesh_10','mesh_6'),('mesh','mesh_6','mesh_4'),
                      ('loft','mesh_4','loft_40'),('loft','loft_40','loft_80'),
                      ('frequency','frequency_101','loft_80')]:
        change=curve_change(curves[a],curves[b])
        if kind=='frequency':
            change['ripple_db']=abs(ripple(curves[a],p['target_band_hz'])-ripple(curves[b],p['target_band_hz']))
        passed=all(np.isfinite(value) and value<=LIMITS[f'{kind}_{key}'] for key,value in change.items())
        results.append(dict(kind=kind,coarse=a,fine=b,changes=change,passed=bool(passed)))
    verify(out,analysis_only=reanalyze)
    clean_source_revision()
    if source_identity()!=analysis_source:
        raise ValueError('Analysis source changed during comparison')
    if any(sha(out/name)!=digest for name,digest in evidence['files'].items()):
        raise ValueError('Solve evidence changed during comparison')
    destination=out/('reanalysis.json' if reanalyze else 'comparison.json')
    if destination.exists():
        raise FileExistsError(destination)
    result=dict(scope=p['scope'],protocol_sha256=sha(out/'protocol.json'),
        solve_evidence_sha256=sha(out/'solve-evidence.json'),candidate=candidate,
        host_execution_sha256=sha(out/'host-execution.json'),solver_image_id=p['solver_image_id'],
        analysis_revision=analysis_revision,analysis_source=analysis_source,
        comparison_grid='union_log_frequency',historical_reanalysis=reanalyze,
        previous_comparison_sha256=previous_comparison,
        limits=LIMITS,health=health,comparisons=results,passed=all(r['passed'] for r in results),
        physical_validation_status='experimental_prediction')
    write_json(destination,result)
    print(json.dumps(result,indent=2))
    if not result['passed']:
        raise RuntimeError('Candidate failed fixed resolution gates; preserve results and investigate')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','solve','_solve','compare','reanalyze'])
    parser.add_argument('output',type=Path)
    parser.add_argument('--ranking',type=Path)
    parser.add_argument('--run-dir',type=Path)
    parser.add_argument('--jobs',type=int,default=1)
    parser.add_argument('--driver',type=Path)
    parser.add_argument('--f-low',type=float)
    parser.add_argument('--f-high',type=float)
    parser.add_argument('--candidate-index',type=int,default=0)
    args=parser.parse_args()
    out=args.output.resolve()
    if args.stage=='prepare':
        if any(v is None for v in (args.ranking,args.driver,args.f_low,args.f_high,args.run_dir)):
            parser.error('prepare requires --run-dir, --ranking, --driver, --f-low and --f-high')
        prepare(out,args.ranking,args.driver,args.f_low,args.f_high,args.candidate_index,args.run_dir)
    elif args.stage=='solve':
        solve(out,args.jobs)
    elif args.stage=='_solve':
        solve_in_container(out,args.jobs)
    else:
        compare(out,reanalyze=args.stage=='reanalyze')


if __name__=='__main__':
    main()
