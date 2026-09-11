"""Automated 500–6500 Hz, verified 6.5-inch, annular-horn screening study.

A finite-grid plane-wave screen; NOT FEM, a cone/phase-plug model, or a
prediction of the summed three-way coaxial loudspeaker.
"""
from __future__ import annotations
import argparse
from dataclasses import replace
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import study_t90a_annular as core

ROOT = core.ROOT
INPUT = ROOT / 'data/drivers'
LOW, HIGH = 500., 6500.


def load_candidates(path=INPUT):
    from horn_drivers.loader import load_drivers_raw
    records = [r for r in load_drivers_raw(str(path))
               if r.get('catalogue_status') == 'manufacturer_verified'
               and r.get('nominal_diameter_verified') is True
               and r.get('nominal_diameter') == '6.5in' and r.get('driver_type') == 'cone']
    if not records or len({r['driver_id'] for r in records}) != len(records):
        raise ValueError('A nonempty set of uniquely identified verified 6.5-inch drivers is required')
    for r in records:
        if not r.get('parameter_source'):
            raise ValueError('Manufacturer source required')
    return records, {r['driver_id']: core._driver_from_dict(r) for r in records}


SEALED_IDS = set()


def response(g, f, *, driver, **kwargs):
    # Published resonance/compliance of factory-sealed drivers already includes
    # their rear loading. Adding the study's 2 L box would count it twice.
    if driver.driver_id in SEALED_IDS:
        kwargs['rear_litres'] = None
    return core.response(g, f, driver=driver, **kwargs)


def metrics(f, result):
    if f[0] > LOW or f[-1] < HIGH or not np.all(np.diff(f) > 0):
        raise ValueError('Increasing frequency grid must cover the entire target band')
    # Include exact endpoints even for externally supplied grids.
    fb = np.unique(np.r_[LOW, f[(f > LOW) & (f < HIGH)], HIGH])
    spl = np.interp(fb, f, result['spl'])
    midf = np.geomspace(800, 2000, 101)
    mid = float(np.mean(np.interp(midf, f, result['spl'])))
    return dict(ripple_db=float(np.ptp(spl)), mean_spl_db=float(np.mean(spl)),
                spl_500_db=float(spl[0]), spl_6500_db=float(spl[-1]),
                level_6500_relative_mid_db=float(spl[-1]-mid),
                max_excursion_mm=float(np.max(np.interp(fb, f, result['x_mm']))),
                max_throat_speed_m_s=float(np.max(np.interp(fb, f, result['throat_speed']))))


def sort_key(r):
    return (r['mean_spl_db'] < 100, r['ripple_db'], -r['mean_spl_db'],
            r['mouth_diameter_m']**2*r['length_m'])


def geometry(row):
    return core.Geometry(**{k: row[k] for k in core.Geometry.__dataclass_fields__})


def write_response(path, f, r):
    core.write_csv(path, [dict(frequency_hz=float(freq),spl_db=float(r['spl'][i]),
        pressure_real_pa=float(r['pressure'][i].real),pressure_imag_pa=float(r['pressure'][i].imag),
        impedance_real_ohm=float(r['z_e'][i].real),impedance_imag_ohm=float(r['z_e'][i].imag),
        excursion_peak_mm=float(r['x_mm'][i]),throat_velocity_rms_m_s=float(r['throat_speed'][i]))
        for i,freq in enumerate(f)])


def numerical_checks(g, f, driver):
    coarse = response(g,f,driver=driver,segments=300)
    fine = response(g,f,driver=driver,segments=600)
    densef = np.unique(np.r_[np.geomspace(300,10000,1201),LOW,HIGH])
    dense = response(g,densef,driver=driver,segments=600)
    errors = dict(section_doubling_max_db=float(np.max(abs(coarse['spl']-fine['spl']))),
        frequency_doubling_interpolation_max_db=float(np.max(abs(dense['spl']-np.interp(densef,f,fine['spl'])))),
        relative_real_power_error=float(np.max(abs(fine['radiation_power']-fine['horn_power'])/np.maximum(fine['horn_power'],1e-20))))
    assert errors['section_doubling_max_db'] < .25
    assert errors['frequency_doubling_interpolation_max_db'] < .25
    assert errors['relative_real_power_error'] < 1e-8
    return errors


def plot(out, f, best, drivers, selected, sensitivity):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig, axes = plt.subplots(2,1,figsize=(11,9),layout='constrained')
    for row in best:
        if not row['nominal_band_covers_target'] and '6p5-006' not in row['scenario_id']:
            continue
        r=response(geometry(row),f,driver=drivers[row['scenario_id']],segments=600)
        axes[0].semilogx(f,r['spl'],label=row['label'],ls='--' if row['sensitivity_only'] else '-')
    for ax in axes:
        ax.axvspan(LOW,HIGH,color='#e6edf3',alpha=.4)
        ax.axvline(LOW,color='grey',ls=':');ax.axvline(HIGH,color='grey',ls=':')
        ax.grid(alpha=.2);ax.set(xlim=(300,10000),ylabel='Predicted dB SPL at 2.83 V / 1 m')
    axes[0].set_title('Nominal band candidates + REDCATT comparison — reduced-order screen')
    axes[0].legend(fontsize=8,ncol=2)
    for cc in (10,20,40,60):
        r=response(geometry(selected),f,driver=drivers[selected['scenario_id']],front_cc=cc,segments=600)
        axes[1].semilogx(f,r['spl'],label=f'{cc} cm³ front chamber')
    axes[1].set(title=f"{selected['label']}: chamber-volume sensitivity",xlabel='Frequency (Hz)')
    axes[1].legend(fontsize=9)
    fig.savefig(out/'responses.png',dpi=160);plt.close(fig)
    g=geometry(selected); z=np.linspace(0,g.length_m,401)*1000
    radius=np.array([g.radius(v/1000) for v in z])*1000
    fig,ax=plt.subplots(figsize=(9,7))
    fig.subplots_adjust(left=.15,right=.75,bottom=.23,top=.87)
    ax.fill_between(z,radius,radius+5,color='#354a5f');ax.fill_between(z,-radius,-radius-5,color='#354a5f')
    ax.fill_between([0,100],-32,32,color='#b8a282',label='Reserved central HF housing, Ø64 × 100 mm')
    ax.plot([100,100],[-32,32],color='#8a4c1b',lw=3)
    ax.annotate('HF face position\n(no HF source simulated)',xy=(100,0),xytext=(112,40),arrowprops={'arrowstyle':'->'},fontsize=9)
    fig.text(.10,.06,'Inlet at z=0. Cone adapter, rear chamber and outer 15-inch horn are not drawn.\nWall thickness shown is illustrative; STEP contains the acoustic air volume.',fontsize=9)
    ax.set(aspect='equal',xlabel='Distance from annular inlet (mm)',ylabel='Radius (mm)',
           title=f"{selected['label']}: {g.profile}, mouth {g.mouth_diameter_m*1000:g} mm, length {g.length_m*1000:g} mm")
    ax.legend(loc='upper left',fontsize=8);ax.grid(alpha=.2)
    fig.savefig(out/'section.png',dpi=160);plt.close(fig)


def report(out, summary, records):
    s=summary; w=s['selected']; c=next(r for r in s['best_by_scenario'] if r['scenario_id']=='celestion-cf0617m')
    lines=['# 6.5-inch midrange search: 500 Hz–6.5 kHz','',
      '**Completed automated reduced-order screen. No physical crossover or complete coaxial system is validated.**','',
      f"Screened **{s['search']['geometry_count']} geometries × {s['search']['verified_driver_count']} verified 6.5-inch drivers**, plus one separately labelled Celestion inductance scenario: **{s['search']['scenario_count']} calculations**. The target remains exactly 500–6500 Hz.",'',
      '![Predicted responses](responses.png)','',
      '## Results','',
      'All figures are unfiltered, on-axis, 2.83 V RMS at 1 m. These are model rankings, not sound-quality rankings. Each row uses its own best horn.','',
      '| Driver/scenario | Profile | Mouth / length (mm) | Equivalent / outer throat (mm) | Variation (dB) | Mean SPL (dB) | SPL at 6.5 kHz (dB) |',
      '|---|---|---:|---:|---:|---:|---:|']
    for r in s['best_by_scenario']:
        lines.append(f"| {r['label']} | {r['profile']} | {r['mouth_diameter_m']*1000:g} / {r['length_m']*1000:g} | {r['equivalent_throat_diameter_m']*1000:g} / {r['outer_throat_diameter_mm']:.1f} | {r['ripple_db']:.2f} | {r['mean_spl_db']:.2f} | {r['spl_6500_db']:.2f} |")
    lines += ['',f"**Selected nominal model: {w['label']}**, {w['profile']}, {w['mouth_diameter_m']*1000:g} mm mouth and {w['length_m']*1000:g} mm acoustic length. Open annular throat {w['throat_area_cm2']:.2f} cm², compression {w['compression_ratio']:.2f}:1, radial gap {w['radial_gap_mm']:.2f} mm around the 64 mm housing.",'',
      f"The selected model {'meets' if s['passes_nominal_screen'] else 'does not meet'} the declared screen of ≤6 dB peak-to-peak variation and ≥100 dB mean output. Neither threshold establishes maximum output, distortion or physical suitability.",'',
      'Drivers with unknown published bands or upper limits below 6.5 kHz remain exploratory comparisons and cannot win the nominal shortlist. The Celestion 0.29 mH scenario is a sensitivity case, not a sixth driver and not eligible for nominal selection.','',
      '![Selected acoustic section](section.png)','',
      '[Selected acoustic-air STEP, mm](selected/acoustic-air-mm.step) · [Housing envelope STEP, mm](selected/housing-concept-mm.step) · [Celestion nominal geometry STEP, mm](celestion/acoustic-air-mm.step) · [Full search CSV](search.csv) · [Driver inputs](drivers.json) · [Machine-readable study](study.json)','',
      '## Search and selection','',
      'Four profiles: conical, exponential, hyperbolic and OS. Equivalent throat diameters: 42, 47, 52, 60, 68 and 76 mm. Mouth diameters: 190, 210, 245 and 280 mm. Acoustic lengths: 125, 150, 200, 250 and 300 mm. Every geometry retains a 64 mm diameter × 100 mm central body. Outer throat diameter is sqrt(equivalent² + 64²), never a circular hole smaller than the tweeter.','',
      'Nominal front volume is 20 cm³ and external sealed rear volume 2 litres for open-basket drivers. Factory-sealed Beyma 6MCF200Nd and PHL 1660NdM-SQ2 use their published combined driver/rear-load parameters with no additional rear chamber. Neither is a demonstrated cone adapter or packaged rear chamber. Selection first requires nominal manufacturer band coverage and excludes sensitivity scenarios; among those it prefers mean SPL ≥100 dB, then minimum unfiltered peak-to-peak variation. If none meet the mean level, it selects minimum variation and reports the miss. All 480 geometries are evaluated at 300 sections and 603 frequency points, including exact band endpoints; winners are re-evaluated at 600 sections. Finite grid only, no global optimum claimed.','',
      '## Celestion uncertainty','',
      'Celestion currently publishes Le=1.73 mH. The 2015 Voice Coil test reports about 0.28–0.30 mH in its inductance-versus-displacement analysis. These values have different measurement contexts and cannot be treated as interchangeable measured broadband impedances. Both are run with an explicitly simplified constant-Le model; 0.29 mH is an uncertainty probe, never a silent correction.','',
      f"With the nominal value the best Celestion geometry has {c['ripple_db']:.2f} dB variation and {c['spl_6500_db']:.2f} dB at 6.5 kHz. See the separate sensitivity row before interpreting its relative ranking. Celestion's current 2.7 mm Xmax definition includes a gap allowance; geometric coil overhang is 1.2 mm. No maximum-SPL claim is made.",'',
      '## Robustness and model limits','',
      'Front chambers of 10/20/40/60 cm³ and rear volumes of 0.5/1/2/4 L were checked for each winning scenario. Factory-sealed drivers ignore the external rear-volume sweep; their identical rows are intentional. Full values are in study.json. Selection uses the declared nominal volumes; sensitivity results are not folded into a hidden score. Numerical section/frequency convergence and real-power conservation passed for every reported scenario. STEP export/re-import checks confirm connected air volumes at millimetre scale.','',
      f"For context, this reduced-order model predicts {s['baseline_diagnostic']['tmm_minus_preserved_fem_at_6500_db']:.2f} dB more output at 6.5 kHz than preserved modal FEM for the original unobstructed 6NMB420 horn with matched assumptions. This is evidence of upper-band model uncertainty, not a correction factor for these new horns.",'',
      '- Plane-wave, lossless acoustic network with a rigid-piston motor and constant inductance. No measured cone breakup, phase-equalising channels, viscothermal losses or nonlinear distortion.',
      '- Central HF housing is a rigid obstruction. HF output, loading, time alignment and crossover summation are absent. A T90A-shaped envelope is retained for clearance only: the T90A itself is specified for ≥7 kHz and is not qualified for this 6.5 kHz target.',
      '- The existing production modal FEM gate rejects annular inlet ports; Docker is unavailable. This run is not production FEM.',
      '- The 15-inch outer horn is unspecified and not simulated. Driver frame diameters, rear chambers, mounting walls and support struts must be included before claiming that the mid assembly fits or preserves the outer horn response.',
      '- The selected STEP starts at the annular inlet and excludes the real cone-following adapter, retention, cable route and supports. A suitable prototype phase plug still needs the actual cone/dustcap contour.',
      '- Driver sensitivity and power ratings are not maximum system output. Manufacturer nominal upper ranges do not validate 6.5 kHz horn operation.','',
      '## Reproduce','',
      'From the repository root, with the same Python 3.12 environment as the earlier annular study:','',
      '```sh','.venv-t90a/bin/python scripts/study_6p5_mid.py --output results/6p5-mid-500-6500',
      '.venv-t90a/bin/python -m pytest scripts/tests/test_t90a_annular.py scripts/tests/test_6p5_mid.py -q','```','',
      '[Dependency versions](requirements.txt). The manifest hashes executed repository sources, driver inputs and all output files. The source Git revision is supplementary; file hashes include working-tree changes.','',
      '## Manufacturer sources','']
    for r in records:
        lines.append(f"- [{r['manufacturer']} {r['model_name']}]({r['parameter_source']}): verified 6.5-inch nominal size; basket envelope {r['physical_envelope_mm']['maximum_width']} × {r['physical_envelope_mm']['depth']} mm.")
    lines += ['- [Celestion independent measurements and alternative inductance context](https://audioxpress.com/article/test-bench-celestion-cf0617m-prosound-midrange-driver)','']
    (out/'README.md').write_text('\n'.join(lines))


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    records,drivers=load_candidates()
    byid={r['driver_id']:r for r in records}
    SEALED_IDS.clear(); SEALED_IDS.update(r['driver_id'] for r in records if r.get('factory_sealed_back'))
    drivers['celestion-cf0617m-le029']=replace(drivers['celestion-cf0617m'],le_h=.00029)
    f=np.unique(np.r_[np.geomspace(300,10000,601),LOW,HIGH])
    rows=[]
    for profile in ('conical','exponential','hyperbolic','os'):
        for deq,mouth,length in itertools.product((.042,.047,.052,.060,.068,.076),(.190,.210,.245,.280),(.125,.150,.200,.250,.300)):
            g=core.Geometry(profile,deq,mouth,length)
            tr=core.transfer(g,f,300)
            for scenario,driver in drivers.items():
                sensitivity=scenario.endswith('-le029');r=byid[driver.driver_id]
                rows.append(dict(scenario_id=scenario,driver_id=driver.driver_id,
                    label=f"{r['manufacturer']} {r['model_name']}"+(' — Le 0.29 mH sensitivity' if sensitivity else ''),
                    sensitivity_only=sensitivity,nominal_band_covers_target=r.get('usable_f_low_hz') is not None and r.get('usable_f_high_hz') is not None and r['usable_f_low_hz']<=LOW and r['usable_f_high_hz']>=HIGH,
                    **g.dimensions(driver),**metrics(f,response(g,f,tr=tr,driver=driver))))
        print(f'{profile}: {len(rows)} driver/geometry scenarios completed',flush=True)
    rows.sort(key=sort_key)
    core.write_csv(out/'search.csv',rows)
    best=[]
    for scenario in drivers:
        row=next(r for r in rows if r['scenario_id']==scenario).copy()
        row.update(metrics(f,response(geometry(row),f,driver=drivers[scenario],segments=600)))
        best.append(row)
    selected=min((r for r in best if not r['sensitivity_only'] and r['nominal_band_covers_target']),key=sort_key)
    sensitivity={};checks={}
    for row in best:
        scenario=row['scenario_id'];g=geometry(row);driver=drivers[scenario]
        checks[scenario]=numerical_checks(g,f,driver)
        sensitivity[scenario]=dict(front_cc={str(v):metrics(f,response(g,f,driver=driver,front_cc=v,segments=600)) for v in (10,20,40,60)},
            rear_litres={str(v):metrics(f,response(g,f,driver=driver,rear_litres=v,segments=600)) for v in (.5,1,2,4)}, factory_sealed_back=driver.driver_id in SEALED_IDS)
        write_response(out/f'{scenario}-response.csv',f,response(g,f,driver=driver,segments=600))
    cad={}
    for label,row in [('selected',selected),('celestion',next(r for r in best if r['scenario_id']=='celestion-cf0617m'))]:
        path=out/label;path.mkdir(exist_ok=True)
        cad[label]=core.export_cad(geometry(row),path)
        core.write_csv(path/'profile-mm.csv',[dict(z_mm=z*1000,outer_radius_mm=geometry(row).radius(z)*1000,inner_radius_mm=32 if z<.1 else 0) for z in np.unique(np.r_[np.linspace(0,row['length_m'],301),.1-1e-9,.1])])
    ref=core.Geometry('os',.076,.24480268,.3,0,0)
    ref_spl=core.response(ref,f,front_cc=0,rear_litres=None,segments=600)['spl']
    preserved=np.genfromtxt(ROOT/'examples/6nmb420-320-5000/response.csv',delimiter=',',names=True)
    summary=dict(status='completed_reduced_order_screen_not_validated_design',selected=selected,best_by_scenario=best,
        passes_nominal_screen=bool(selected['ripple_db']<=6 and selected['mean_spl_db']>=100),
        search=dict(target_hz=[LOW,HIGH],nominal_diameter_inches=6.5,verified_driver_count=len(records),geometry_count=len(rows)//len(drivers),scenario_count=len(rows),front_volume_cc=20,rear_volume_litres=2,body_diameter_mm=64,body_length_mm=100),
        sensitivity=sensitivity,numerical_checks=checks,cad=cad,
        baseline_diagnostic=dict(tmm_minus_preserved_fem_at_6500_db=float(np.interp(HIGH,f,ref_spl)-np.interp(HIGH,preserved['frequency_hz'],preserved['spl_db_2p83vrms_1m'])),scope='Original unobstructed horn only; not a correction to new candidates'),
        compact_reference=min((r for r in rows if not r['sensitivity_only'] and r['nominal_band_covers_target'] and r['mouth_diameter_m']<=.210 and r['length_m']<=.200),key=sort_key))
    (out/'study.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'drivers.json').write_text(json.dumps({'schema_version':2,'drivers':records},indent=2)+'\n')
    (out/'requirements.txt').write_text((ROOT/'examples/6nmb420-t90a-annular/requirements.txt').read_text())
    plot(out,f,best,drivers,selected,sensitivity)
    report(out,summary,records)
    sources={Path(__file__).resolve(),Path(core.__file__).resolve(),core.DRIVER_PATH,*INPUT.glob('*/*.json'),ROOT/'data/catalogue-audit.json',
        ROOT/'examples/6nmb420-320-5000/response.csv',ROOT/'scripts/tests/test_6p5_mid.py',ROOT/'scripts/tests/test_t90a_annular.py'}
    for mod in tuple(sys.modules.values()):
        filename=getattr(mod,'__file__',None)
        if filename:
            p=Path(filename).resolve()
            if p.is_relative_to(ROOT/'packages') and p.is_file():sources.add(p)
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,text=True,capture_output=True,check=True).stdout.strip()
    manifest=dict(source_git_head=head,source_identity_note='Executed file hashes include uncommitted changes.',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)},
        output_sha256={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('SELECTED',json.dumps(selected),flush=True)
    print(f'Results: {out}',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'examples/6p5-mid-500-6500')
    run(parser.parse_args().output)
