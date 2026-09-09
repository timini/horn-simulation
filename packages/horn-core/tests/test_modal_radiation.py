"""Independent piston limits and structural checks for the aperture operator."""
import numpy as np
import pytest
from scipy.special import struve, j1, j0
from horn_core.modal_radiation import modal_radiation_impedance, modal_baffled_on_axis, radial_roots, mode_spectrum
from horn_core.acoustics import baffled_piston_on_axis


@pytest.mark.parametrize('ka',[.1,.5,1.,3.,8.,20.])
def test_plane_mode_matches_analytic_piston(ka):
    actual=modal_radiation_impedance(ka,4)[0,0]
    expected=1-j1(2*ka)/ka+1j*struve(1,2*ka)/ka
    assert actual==pytest.approx(expected,rel=8e-6,abs=1e-7)


def test_reciprocity_passivity_and_spectral_tail_convergence():
    for ka in (.5,3.,8.):
        matrix=modal_radiation_impedance(ka,8)
        finer=modal_radiation_impedance(ka,8,spectral_limit=4096.)
        np.testing.assert_allclose(matrix,matrix.T,atol=1e-14)
        assert np.linalg.eigvalsh(matrix.real).min()>-1e-13
        assert np.linalg.eigvalsh(matrix.imag).min()>-1e-13
        assert np.max(abs(matrix-finer))<1e-6


def test_bessel_roots_have_finite_correct_limits():
    roots=radial_roots(5)
    for delta in (0.,1e-10,-1e-10):
        actual=np.diag(mode_spectrum(roots+delta,roots))
        np.testing.assert_allclose(actual,j0(roots)/2,atol=1e-12)


def test_uniform_mode_observer_matches_exact_piston_at_near_and_far_points():
    f=np.geomspace(100.,5000.,20);radius=.08
    v=np.full((len(f),1),.3+.1j)
    for distance in (.01,1.,10.):
        expected=baffled_piston_on_axis(f,v[:,0]*np.pi*radius**2,np.pi*radius**2,distance)
        np.testing.assert_allclose(modal_baffled_on_axis(f,v,radius,distance),expected,rtol=1e-11,atol=1e-11)
