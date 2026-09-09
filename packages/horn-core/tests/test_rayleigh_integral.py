"""Independent boundary-integral checks of the uniform infinite-baffle model.

Integrate monopole Green functions rather than reuse the production closed
forms. This does not validate replacing a nonuniform horn mouth by a piston.
Rayleigh formulation: https://euphonics.org/4-3-2-the-rayleigh-integral-and-the-baffled-piston/
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
import pytest
from horn_core.acoustics import baffled_piston_on_axis, C0, RHO0
from horn_core.webster import piston_radiation_impedance


def quadrature(count, low, high):
    nodes, weights = leggauss(count)
    return (nodes+1)*(high-low)/2+low, weights*(high-low)/2


@pytest.mark.parametrize('radius', [.005, .05, .2])
@pytest.mark.parametrize('relative_distance', [.05, 1., 20.])
def test_on_axis_pressure_matches_integrated_monopoles(radius, relative_distance):
    ka = np.geomspace(.01, 20., 31)
    k = ka/radius
    frequency = k*C0/(2*np.pi)
    velocity = .003+.004j  # RMS, including phase
    area = np.pi*radius**2
    distance = relative_distance*radius
    expected = baffled_piston_on_axis(frequency, np.full(len(k),velocity*area), area, distance)
    answers = []
    for count in [64, 128]:
        radial, weights = quadrature(count, 0., radius)
        separation = np.sqrt(distance**2+radial**2)
        # p = i omega rho/(2 pi) integral_S v exp(-ikR)/R dS.
        integral = np.sum(np.exp(-1j*k[:,None]*separation)*radial/separation*weights,axis=1)*2*np.pi
        answers.append(1j*(k*C0)*RHO0*velocity/(2*np.pi)*integral)
    # Check quadrature convergence before comparing the production expression.
    scale = RHO0*C0*abs(velocity)
    np.testing.assert_allclose(answers[0],answers[1],rtol=1e-8,atol=scale*1e-10)
    np.testing.assert_allclose(answers[1],expected,rtol=1e-8,atol=scale*1e-10)


def test_surface_averaged_load_matches_double_rayleigh_integral():
    ka = np.geomspace(.01,20.,31)
    answers=[]
    for count in [128,256]:
        separation, weights = quadrature(count,0.,2.)
        # For unit radius, the convolution of two disk indicators is their
        # overlap area. It reduces the double surface integral to one line;
        # the separation in the Green function cancels the polar Jacobian.
        overlap = 2*np.arccos(separation/2)-separation*np.sqrt(4-separation**2)/2
        integral=np.sum(np.exp(-1j*ka[:,None]*separation)*overlap*weights,axis=1)
        answers.append(1j*ka/np.pi*integral)
    expected=np.array([piston_radiation_impedance(value,1.) for value in ka])
    np.testing.assert_allclose(answers[0],answers[1],rtol=1e-8,atol=1e-10)
    np.testing.assert_allclose(answers[1],expected,rtol=1e-8,atol=1e-10)
