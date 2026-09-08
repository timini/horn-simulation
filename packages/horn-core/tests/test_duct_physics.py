import numpy as np
import pytest
from horn_core.duct import AirProperties, circular_duct_properties, circular_pipe_radiation
from horn_core.webster import compute_horn_transfer_tmm


def test_dissipative_propagation_and_lossless_limit():
    f = np.geomspace(100,10000,50)
    for model in ['kirchhoff','boundary_layer']:
        k, z = circular_duct_properties(f,.007,loss_model=model)
        assert np.all(k.imag < 0)  # exp(-ikx) decays, never amplifies
        assert np.all(z.real > 0)
        k0,z0 = circular_duct_properties(f,100.,loss_model=model)
        np.testing.assert_allclose(k0,2*np.pi*f/343,rtol=1e-5)
        np.testing.assert_allclose(z0,1.225*343,rtol=1e-5)


def test_finite_flange_passivity_and_unflanged_low_frequency_limit():
    ka = np.geomspace(1e-5,1.49,100)
    for width in [0,.002,.007]:
        z = circular_pipe_radiation(ka/.007,.007,width)
        assert np.all(z.real >= 0)
        assert np.all(np.abs((z-1)/(z+1)) <= 1+1e-12)
    z = circular_pipe_radiation(1e-4/.007,.007)
    assert z.imag/1e-4 == pytest.approx(.6133,rel=1e-6)
    assert z.real/1e-8 == pytest.approx(.25,rel=1e-5)


def test_closed_duct_has_zero_flow_and_known_lossless_impedance():
    f = np.array([200.,300.,400.])
    t = compute_horn_transfer_tmm(f,lambda _: .007,.18,.007,.007,radiation_model='closed')
    np.testing.assert_allclose(t['z_real']+1j*t['z_imag'],-1j*1.225*343/np.tan(2*np.pi*f*.18/343),atol=1e-8)
    assert not np.any(t['mouth_volume_velocity_transfer'])


def test_thin_layer_and_air_domains_fail_closed():
    with pytest.raises(ValueError,match='depth / radius'):
        circular_duct_properties([20,30],.0001,loss_model='boundary_layer')
    with pytest.raises(ValueError):
        AirProperties(viscosity=-1)
    with pytest.raises(ValueError,match='ka'):
        circular_pipe_radiation(1.6/.007,.007)
    with pytest.raises(ValueError,match='width / radius'):
        circular_pipe_radiation(100,.007,.05)
