"""Reproducible reduced-order 6NMB420 / recessed T90A packaging study.

This is a Webster/plane-wave screen, NOT the production FEM path or a
simulation of cone breakup, phase-equalising channels, or the HF source.
The cylinder terminates inside the horn. Its end is rigid in the mid model.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, asdict, replace
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
for package in ('horn-core', 'horn-drivers', 'horn-analysis'):
    sys.path.insert(0, str(ROOT / 'packages' / package / 'src'))

import numpy as np
from horn_core.profiles import get_radius_func
from horn_core.webster import piston_radiation_impedance, compute_horn_transfer_tmm
from horn_core.acoustics import baffled_piston_on_axis, pressure_level
from horn_drivers.loader import _driver_from_dict
from horn_analysis.transfer_function import compute_driver_operating_point

RHO, C = 1.225, 343.0
DRIVER_PATH = ROOT / 'examples/6nmb420-320-5000/driver.json'
BASE_DRIVER = _driver_from_dict(json.loads(DRIVER_PATH.read_text()))


@dataclass(frozen=True)
class Geometry:
    profile: str
    equivalent_throat_diameter_m: float
    mouth_diameter_m: float
    length_m: float
    body_diameter_m: float = 0.064
    body_length_m: float = 0.100

    @property
    def outer_throat_radius(self):
        return np.hypot(self.equivalent_throat_diameter_m, self.body_diameter_m) / 2

    @property
    def throat_area(self):
        return np.pi * (self.equivalent_throat_diameter_m / 2) ** 2

    def radius(self, z):
        return get_radius_func(self.profile, self.outer_throat_radius,
                               self.mouth_diameter_m / 2, self.length_m)(z)

    def dimensions(self):
        return dict(asdict(self), outer_throat_diameter_mm=2000*self.outer_throat_radius,
                    radial_gap_mm=1000*(self.outer_throat_radius-self.body_diameter_m/2),
                    throat_area_cm2=self.throat_area*1e4,
                    compression_ratio=BASE_DRIVER.sd_m2/self.throat_area)


@lru_cache(maxsize=64)
def mouth_load(frequency_tuple, radius):
    k = 2*np.pi*np.array(frequency_tuple)/C
    return np.array([piston_radiation_impedance(x, radius) for x in k])*RHO*C/(np.pi*radius**2)


def transfer(g, f, segments=300):
    """Pressure/volume-flow matrix; integrate each side of body end separately.

    Continuity of p and U at the area step is a plane-mode approximation.
    End scattering and non-plane modes are deliberately not represented.
    """
    if not (0 <= g.body_length_m < g.length_m and g.body_diameter_m >= 0
            and g.equivalent_throat_diameter_m > 0 and g.mouth_diameter_m > 2*g.outer_throat_radius):
        raise ValueError('Invalid or obstructed annular horn geometry')
    if segments < 2:
        raise ValueError('At least two sections required')
    f = np.asarray(f)
    k = 2*np.pi*f/C
    mat = np.tile(np.eye(2, dtype=complex), (len(f), 1, 1))
    breaks = [0., g.body_length_m, g.length_m] if g.body_length_m > 0 else [0., g.length_m]
    for start, end in zip(breaks[:-1], breaks[1:]):
        n = max(1, int(np.ceil(segments*(end-start)/g.length_m)))
        dz = (end-start)/n
        co, si = np.cos(k*dz), np.sin(k*dz)
        for z in start+(np.arange(n)+.5)*dz:
            inner = g.body_diameter_m/2 if z < g.body_length_m else 0.
            area = np.pi*(g.radius(z)**2-inner**2)
            if area <= 0:
                raise ValueError('Blocked section')
            zc = RHO*C/area
            seg = np.empty_like(mat)
            seg[:, 0, 0] = seg[:, 1, 1] = co
            seg[:, 0, 1] = 1j*zc*si
            seg[:, 1, 0] = 1j*si/zc
            mat = mat @ seg
    load = mouth_load(tuple(f), g.mouth_diameter_m/2)
    pin = mat[:, 0, 0]*load+mat[:, 0, 1]
    uin = mat[:, 1, 0]*load+mat[:, 1, 1]
    return dict(z_acoustic=pin/uin, uout_over_uin=1/uin, load=load, matrix=mat)


def response(g, f, front_cc=20., rear_litres=2., segments=300, tr=None):
    """Ideal rigid-piston driver, shunt front compliance, sealed rear volume.

    Vf is a sensitivity parameter, not a measured/packaged compression chamber.
    The published Mms fallback is retained; added sealed-box stiffness is
    explicit and does not claim a separately measured diaphragm/rear air mass.
    """
    if front_cc < 0 or (rear_litres is not None and rear_litres <= 0):
        raise ValueError('Invalid chamber volume')
    tr = transfer(g, f, segments) if tr is None else tr
    omega = 2*np.pi*f
    chamber_factor = 1+1j*omega*(front_cc*1e-6/(RHO*C*C))*tr['z_acoustic']
    zfront = tr['z_acoustic']/chamber_factor
    driver = BASE_DRIVER
    if rear_litres is not None:
        compliance = 1/(1/driver.cms_m_per_n+RHO*C*C*driver.sd_m2**2/(rear_litres*1e-3))
        driver = replace(driver, cms_m_per_n=compliance)
    point = compute_driver_operating_point(driver, f, (zfront*g.throat_area).real,
                                          (zfront*g.throat_area).imag, g.throat_area, 2.83)
    ucone = point['velocity_rms']*driver.sd_m2
    uin = ucone/chamber_factor
    uout = uin*tr['uout_over_uin']
    pressure = baffled_piston_on_axis(f, uout, np.pi*(g.mouth_diameter_m/2)**2, 1.)
    rad_power = abs(uout)**2*tr['load'].real
    result = dict(pressure=pressure, spl=pressure_level(pressure),
                  z_e=point['electrical_impedance'], x_mm=point['displacement_peak_m']*1000,
                  throat_speed=abs(uin)/g.throat_area,
                  radiation_power=rad_power, horn_power=point['horn_power_w'],
                  input_power=point['input_power_w'])
    return result


def metrics(f, result):
    spl = result['spl']
    band = (f >= 320) & (f <= 7000)
    mid = (f >= 800) & (f <= 2000)
    at7 = float(np.interp(7000, f, spl))
    return dict(ripple_db=float(np.ptp(spl[band])), mean_spl_db=float(np.mean(spl[band])),
                spl_7000_db=at7, level_7000_relative_mid_db=at7-float(np.mean(spl[mid])),
                max_excursion_mm=float(np.max(result['x_mm'][band])),
                max_throat_speed_m_s=float(np.max(result['throat_speed'][band])),
                max_input_power_w=float(np.max(result['input_power'][band])))


def checks(g, f):
    """Independent limiting solution, existing implementation and conservation."""
    fg = np.geomspace(250, 10000, 101)
    ref = Geometry('os', .076, .24480268, .3, 0., 0.)
    custom = transfer(ref, fg, 300)
    existing = compute_horn_transfer_tmm(fg, ref.radius, ref.length_m,
                ref.outer_throat_radius, ref.mouth_diameter_m/2, n_segments=300)
    expected = (existing['z_real']+1j*existing['z_imag'])/ref.throat_area
    agreement = float(np.max(abs(custom['z_acoustic']-expected)/np.maximum(abs(expected), 1e-12)))
    # The no-flare annular limit is obtained to machine precision using a
    # virtually equal mouth radius, body spanning almost all the tube.
    tube = Geometry('conical', .047, np.hypot(.047,.064)+1e-12, .1, .064, .1-1e-12)
    raw = transfer(tube, fg, 200)
    k = 2*np.pi*fg/C
    zc = RHO*C/tube.throat_area
    exact = zc*(raw['load']+1j*zc*np.tan(k*.1))/(zc+1j*raw['load']*np.tan(k*.1))
    tube_error = float(np.max(abs(raw['z_acoustic']-exact)/abs(exact)))
    coarse, fine = response(g, f, segments=300), response(g, f, segments=600)
    grid = np.unique(np.r_[np.geomspace(226.274, 12000, 1201), 320, 7000])
    dense = response(g, grid, segments=600)
    resolution = float(np.max(abs(coarse['spl']-fine['spl'])))
    frequency_error = float(np.max(abs(dense['spl']-np.interp(grid, f, fine['spl']))))
    conservation = float(np.max(abs(fine['radiation_power']-fine['horn_power'])/np.maximum(fine['horn_power'], 1e-20)))
    determinant_error = float(np.max(abs(np.linalg.det(transfer(g, f)['matrix'])-1)))
    result = dict(existing_tmm_relative_impedance_error=agreement, analytic_tube_relative_error=tube_error,
                  segments_300_to_600_max_spl_db=resolution,
                  frequency_grid_doubling_max_interpolation_db=frequency_error,
                  relative_power_conservation_error=conservation, transfer_determinant_error=determinant_error,
                  scope='Numerical reduced-order checks only; no FEM or physical validation')
    assert agreement < 1e-8 and tube_error < 1e-7 and conservation < 1e-8 and determinant_error < 1e-8
    assert resolution < .25 and frequency_error < .25
    return result


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def export_cad(g, out):
    """Millimetre STEP: air passage and separate concept housing, no cone adapter."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.option.setString('Geometry.OCCTargetUnit', 'MM')
    gmsh.model.add('t90a_annular_concept_mm')
    occ = gmsh.model.occ
    wires=[]
    for z in np.linspace(0, g.length_m, 81):
        wires.append(occ.addWire([occ.addCircle(0, 0, z*1000, g.radius(z)*1000)]))
    horn = occ.addThruSections(wires, makeSolid=True)
    body = occ.addCylinder(0,0,0,0,0,g.body_length_m*1000,g.body_diameter_m*500)
    air, _ = occ.cut(horn, [(3,body)])
    occ.synchronize()
    gmsh.write(str(out/'acoustic-air-mm.step'))
    volume = sum(occ.getMass(dim,tag) for dim,tag in air if dim==3)
    gmsh.clear(); gmsh.model.add('concept_housing_mm')
    # 61 mm pocket reserves clearance around 60 mm body; 92 mm includes
    # 87.8 mm nominal envelope/terminals plus assembly allowance. No cable
    # routing, retainers, spokes or cone-following back surface are supplied.
    body = occ.addCylinder(0,0,0,0,0,g.body_length_m*1000,g.body_diameter_m*500)
    pocket = occ.addCylinder(0,0,g.body_length_m*1000-92,0,0,93,30.5)
    occ.cut([(3,body)],[(3,pocket)]); occ.synchronize()
    gmsh.write(str(out/'housing-concept-mm.step'))
    gmsh.clear(); gmsh.model.add('verify_air_mm')
    imported=occ.importShapes(str(out/'acoustic-air-mm.step')); occ.synchronize()
    volumes=[v for v in imported if v[0]==3]
    mass=sum(occ.getMass(*v) for v in volumes)
    bounds=gmsh.model.getBoundingBox(-1,-1)
    # OCC bounding boxes may include spline control points outside the
    # actual surface. Sample the CAD surfaces to verify physical radius.
    sampled=[]
    for _, tag in gmsh.model.getEntities(2):
        lower,upper=gmsh.model.getParametrizationBounds(2,tag)
        uu,vv=np.meshgrid(np.linspace(lower[0],upper[0],81),np.linspace(lower[1],upper[1],81))
        points=np.array(gmsh.model.getValue(2,tag,np.stack([uu.ravel(),vv.ravel()],axis=1).ravel())).reshape(-1,3)
        # For trimmed planar end caps, the rectangular parametrisation also
        # covers points outside the trimmed face; inspect curved faces only.
        if gmsh.model.getType(2,tag) != 'Plane':
            sampled.append(points)
    sampled=np.vstack(sampled)
    max_radius=float(np.max(np.linalg.norm(sampled[:,:2],axis=1)))
    gmsh.finalize()
    assert len(volumes)==1 and abs(mass/volume-1)<1e-6
    assert abs((bounds[5]-bounds[2])-g.length_m*1000)<.01
    assert abs((bounds[3]-bounds[0])-g.mouth_diameter_m*1000)<.01
    assert abs(max_radius-g.mouth_diameter_m*500)<.05
    return dict(units='mm', air_volume_litres=mass/1e6, bbox_mm=list(bounds),
                sampled_surface_max_radius_mm=max_radius,
                single_connected_air_volume=True, reimport_volume_relative_error=abs(mass/volume-1),
                note='Concept packaging only. HF aperture treated as rigid in mid model. No cone adapter or support structure.')


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    f = np.unique(np.r_[np.geomspace(226.274,12000,601),320,7000])
    rows=[]
    for profile in ('os','hyperbolic','exponential','conical'):
        for deq in (.042,.047,.052,.060,.076):
            for length in (.2,.25,.3):
                for mouth in (.18,.21,.245,.28,.32):
                    g=Geometry(profile,deq,mouth,length)
                    m=metrics(f,response(g,f))
                    rows.append(dict(g.dimensions(),**m))
        print(f'Screened {profile}: {len(rows)} geometries',flush=True)
    # Minimum ripple among candidates retaining a broad midrange mean level.
    # This objective is declared, not evidence of a global/physical optimum.
    rows.sort(key=lambda r:(r['mean_spl_db']<100,r['ripple_db']))
    write_csv(out/'search.csv',rows)
    unrestricted=rows[0]
    winner=next(row for row in rows if row['equivalent_throat_diameter_m'] <= .060)
    g=Geometry(**{k:winner[k] for k in Geometry.__dataclass_fields__})
    print('SELECTED',json.dumps(winner),flush=True)
    cases={}
    for cc in (0.,5.,10.,20.,40.,60.):
        cases[f'front_{cc:g}cc']=response(g,f,front_cc=cc,segments=600)
    baseline_g=Geometry('os',.076,.24480268,.3,0.,0.)
    baseline=response(baseline_g,f,front_cc=0,rear_litres=None,segments=600)
    matched_baseline=response(baseline_g,f,front_cc=20,rear_litres=2,segments=600)
    rear_cases={str(v):metrics(f,response(g,f,rear_litres=v,segments=600)) for v in (.5,2.,4.)}
    validation=checks(g,f)
    cad=export_cad(g,out)
    by_throat=[next(row for row in rows if row['equivalent_throat_diameter_m']==d)
               for d in (.042,.047,.052,.060,.076)]
    preserved=np.genfromtxt(ROOT/'examples/6nmb420-320-5000/response.csv',delimiter=',',names=True)
    baseline_bias=float(np.interp(7000,f,baseline['spl'])-
                        np.interp(7000,preserved['frequency_hz'],preserved['spl_db_2p83vrms_1m']))
    summary=dict(selected=winner, unrestricted_reference=unrestricted,best_by_throat=by_throat,
                 front_volume_sensitivity={k:metrics(f,v) for k,v in cases.items()},
                 rear_volume_sensitivity=rear_cases, numerical_checks=validation,cad=cad,
                 baseline_same_front_and_rear=metrics(f,matched_baseline),
                 original_geometry_tmm_minus_preserved_modal_fem_at_7000_db=baseline_bias,
                 search=dict(count=len(rows),front_volume_cc=20,rear_volume_litres=2,
                             target_hz=[320,7000], crossover_target_hz=7000,
                             selection='Among equivalent throats <=60 mm, minimum unfiltered 320–7000 Hz peak-to-peak ripple with mean SPL >=100 dB; otherwise minimum ripple. 76 mm cases retained as references.',
                             body_diameter_mm=64,body_length_mm=100),
                 status='reduced_order_concept_only',
                 meets_6db_unfiltered_ripple_target=winner['ripple_db']<=6,
                 unresolved=['Actual cone/dustcap profile and phase-equalising passages',
                             'T90A response and phase when recessed inside this horn',
                             'Non-plane modes and scattering at centre-body termination',
                             'Cone breakup, distortion and thermal limits',
                             'Buildable chamber volume and excursion clearance',
                             'Separate diaphragm and rear acoustic mass; published Mms retained',
                             'Finite baffle and off-axis summation; no finished crossover'])
    (out/'study.json').write_text(json.dumps(summary,indent=2)+'\n')
    nominal=cases['front_20cc']
    records=[]
    for i,frequency in enumerate(f):
        record=dict(frequency_hz=frequency,mid_pressure_real_pa=nominal['pressure'][i].real,
                    mid_pressure_imag_pa=nominal['pressure'][i].imag,
                    impedance_real_ohm=nominal['z_e'][i].real,impedance_imag_ohm=nominal['z_e'][i].imag,
                    excursion_peak_mm=nominal['x_mm'][i],throat_velocity_rms_m_s=nominal['throat_speed'][i],
                    baseline_same_front_rear_db=matched_baseline['spl'][i],
                    baseline_original_assumptions_tmm_db=baseline['spl'][i])
        record.update({f'spl_{name}_db':value['spl'][i] for name,value in cases.items()})
        records.append(record)
    write_csv(out/'response.csv',records)
    write_csv(out/'profile-mm.csv',[dict(z_mm=z*1000,outer_radius_mm=g.radius(z)*1000,
              inner_radius_mm=g.body_diameter_m*500 if z<g.body_length_m else 0)
              for z in np.unique(np.r_[np.linspace(0,g.length_m,301),g.body_length_m-1e-9,g.body_length_m])])
    plots(g,f,cases,baseline,matched_baseline,out)
    source_paths=[Path(__file__),DRIVER_PATH,ROOT/'packages/horn-core/src/horn_core/webster.py',
                  ROOT/'packages/horn-analysis/src/horn_analysis/transfer_function.py']
    manifest=dict(source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
                  output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!='manifest.json'})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


def plots(g,f,cases,baseline,matched_baseline,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,1,figsize=(10,8),layout='constrained')
    ax=axes[0]
    for name,value in cases.items():
        ax.semilogx(f,value['spl'],label=name.replace('front_','').replace('cc',' cm³'),lw=2 if name=='front_20cc' else 1)
    ax.axvline(7000,color='black',ls=':',label='7 kHz handover target')
    ax.set(xlim=(250,11000),ylim=(70,115),ylabel='Predicted dB SPL, 2.83 V RMS / 1 m',title='Effect of assumed front chamber volume — reduced-order midrange model')
    ax.legend(ncol=4,fontsize=8); ax.grid(alpha=.2)
    ax=axes[1]
    ax.semilogx(f,cases['front_20cc']['spl'],label='Selected annular concept, 20 cm³ front / 2 L rear')
    ax.semilogx(f,matched_baseline['spl'],label='Original 76 mm horn, same 20 cm³ / 2 L assumptions')
    ax.semilogx(f,baseline['spl'],ls='--',alpha=.7,label='Original horn, zero front / no added rear stiffness (TMM)')
    preserved=np.genfromtxt(ROOT/'examples/6nmb420-320-5000/response.csv',delimiter=',',names=True)
    ax.semilogx(preserved['frequency_hz'],preserved['spl_db_2p83vrms_1m'],ls=':',color='grey',label='Preserved original modal FEM; different model and assumptions')
    ax.axvline(7000,color='black',ls=':')
    ax.set(xlim=(250,11000),ylim=(70,115),xlabel='Frequency (Hz)',ylabel='Predicted dB SPL',title='Baseline comparison — no T90A response or crossover summation simulated')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(out/'response-comparison.png',dpi=170);plt.close(fig)
    z=np.linspace(0,g.length_m,401)*1000
    r=np.array([g.radius(x/1000) for x in z])*1000
    fig,ax=plt.subplots(figsize=(11,5),layout='constrained')
    ax.fill_between(z,r,r+5,color='#34495e');ax.fill_between(z,-r,-r-5,color='#34495e')
    ax.fill_between([0,g.body_length_m*1000],-g.body_diameter_m*500,g.body_diameter_m*500,color='#bbb')
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((g.body_length_m*1000-87.8,-30),87.8,60,facecolor='#b78346',edgecolor='black'))
    ax.annotate('T90A envelope: Ø60 × 87.8 mm\nHF face at z = 100 mm',xy=(100,0),xytext=(145,55),arrowprops=dict(arrowstyle='->'))
    ax.annotate(f'Annular inlet: Ø{2*g.outer_throat_radius*1000:.1f} outside / Ø64 inside\n{g.throat_area*1e4:.1f} cm² open; Ø{g.equivalent_throat_diameter_m*1000:.0f} equivalent',xy=(0,g.outer_throat_radius*1000),xytext=(5,105),arrowprops=dict(arrowstyle='->'))
    ax.text(5,-115,'Cone, front chamber, supports and retaining hardware are not drawn.\nHousing pocket: Ø61 × 92 mm. Wall band is illustrative, not exported CAD.',fontsize=9)
    ax.set(xlim=(-10,max(320,g.length_m*1000+10)),ylim=(-140,140),aspect='equal',xlabel='Axial distance from annular inlet (mm)',ylabel='Radius (mm)',title='Recessed T90A packaging concept — not a fabrication drawing')
    ax.grid(alpha=.15);fig.savefig(out/'section.png',dpi=170);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'examples/6nmb420-t90a-annular')
    run(parser.parse_args().output)
