"""Production FEM comparisons with independently known complex tube transfer."""
import numpy as np
import pandas as pd
import pytest
from horn_core.acoustics import C0, RHO0, baffled_piston_on_axis, pressure_level
from horn_solver.solver import run_simulation_from_step
from .conftest import _generate_cylinder_step


@pytest.mark.validation
def test_flanged_tube_impedance_and_output_converge(tmp_path):
    radius, length = .02, .12
    step = tmp_path/'tube.step'
    _generate_cylinder_step(step,radius,length)
    responses=[]
    for h in (.01,.007,.004):
        output=tmp_path/f'tube-{h}.csv'
        run_simulation_from_step(str(step),(200.,1200.),9,{'length':length},str(output),1200.,
                                 mesh_size=h,radiation_model='flanged_piston')
        responses.append(pd.read_csv(output))
    f=responses[-1].frequency.to_numpy(); k=2*np.pi*f/C0
    # Uniform tube impedance transform using an independently assembled closed form.
    from horn_solver.radiation import piston_radiation_impedance
    z_load=np.array([piston_radiation_impedance(ki,radius) for ki in k])
    z_expected=RHO0*C0*(z_load+1j*np.tan(k*length))/(1+1j*z_load*np.tan(k*length))
    df=responses[-1]
    z_actual=df.z_real.to_numpy()+1j*df.z_imag.to_numpy()
    assert np.all(z_actual.real>0), 'Passive termination must not create acoustic energy'
    np.testing.assert_allclose(z_actual,z_expected,rtol=.10,atol=1.)
    pressure_expected=1/(np.cos(k*length)+1j*np.sin(k*length)/z_load)
    mouth_pressure=df.mouth_p_real.to_numpy()+1j*df.mouth_p_imag.to_numpy()
    np.testing.assert_allclose(mouth_pressure,pressure_expected,rtol=.05,atol=.01)
    # Mesh change is checked for the same observable used by optimization.
    observer=[]
    for response in responses:
        u=response.mouth_u_real.to_numpy()+1j*response.mouth_u_imag.to_numpy()
        observer.append(pressure_level(baffled_piston_on_axis(f,u,float(response.mouth_area_m2.iloc[0]))))
    assert np.max(np.abs(observer[-1]-observer[-2])) <= .5
    assert np.max(np.abs(observer[-1]-observer[-2])) < np.max(np.abs(observer[0]-observer[-1]))


@pytest.mark.validation
def test_plane_wave_phase_and_split_sweep(tmp_path):
    radius,length=.02,.12
    step=tmp_path/'tube.step';_generate_cylinder_step(step,radius,length)
    all_frames=[]
    for name,lo,hi,count in [('whole',200.,800.,5),('low',200.,400.,3),('high',400.,800.,3)]:
        out=tmp_path/f'{name}.csv'
        run_simulation_from_step(str(step),(lo,hi),count,{'length':length},str(out),800.,mesh_size=.005,radiation_model='plane_wave')
        all_frames.append(pd.read_csv(out))
    whole=all_frames[0]
    expected=np.exp(-1j*2*np.pi*whole.frequency.to_numpy()*length/C0)
    measured=whole.mouth_p_real.to_numpy()+1j*whole.mouth_p_imag.to_numpy()
    np.testing.assert_allclose(measured,expected,atol=.03)
    joined=pd.concat(all_frames[1:]).drop_duplicates('frequency').sort_values('frequency')
    np.testing.assert_allclose(joined.spl,whole.spl,atol=.5)


@pytest.mark.validation
def test_annular_inlet_area_and_velocity_coupling(tmp_path):
    import gmsh
    from horn_core.parameters import DriverParameters
    from horn_analysis.transfer_function import compute_driver_response
    outer,inner,length=.025,.02,.08
    area=np.pi*(outer**2-inner**2)
    step=tmp_path/'annulus.step'
    gmsh.initialize()
    try:
        gmsh.model.add('annular_air_volume')
        shell=gmsh.model.occ.addCylinder(0,0,0,0,0,length,outer)
        plug=gmsh.model.occ.addCylinder(0,0,0,0,0,length,inner)
        gmsh.model.occ.cut([(3,shell)],[(3,plug)])
        gmsh.model.occ.synchronize();gmsh.write(str(step))
    finally:
        gmsh.finalize()
    run_simulation_from_step(str(step),(200,400),3,{'length':length},str(tmp_path/'a.csv'),400,mesh_size=.006,radiation_model='plane_wave')
    a=pd.read_csv(tmp_path/'a.csv');f=a.frequency.to_numpy()
    np.testing.assert_allclose(a.inlet_area_m2,area,rtol=1e-8)
    np.testing.assert_allclose(a.mouth_area_m2,area,rtol=1e-8)
    np.testing.assert_allclose(a.z_real+1j*a.z_imag,RHO0*C0,rtol=.08)
    driver=DriverParameters('fixture','Test','Motor',200,6,5,.001,.004,.0001,qms=5,qes=.4)
    pressure=compute_driver_response(driver,f,a.z_real.to_numpy(),a.z_imag.to_numpy(),area)
    run_simulation_from_step(str(step),(200,400),3,{'length':length},str(tmp_path/'b.csv'),400,mesh_size=.006,radiation_model='plane_wave',bc_mode='neumann',driver=driver,throat_area=area,z_horn_initial={'frequencies':f,'z_real':a.z_real.to_numpy(),'z_imag':a.z_imag.to_numpy()})
    b=pd.read_csv(tmp_path/'b.csv')
    np.testing.assert_allclose(b.mouth_p_real+1j*b.mouth_p_imag,pressure*(a.mouth_p_real+1j*a.mouth_p_imag),rtol=.04)
