"""Axisymmetric Rayleigh radiation on a circular aperture in an infinite baffle.

Modes are J0(mu_n*r/a)/J0(mu_n), with mu_0=0 and J1(mu_n)=0.
They have equal area norm S. The Hankel spectrum of a mode, divided by
2*pi*a**2, is t*J1(t)/(t**2-mu_n**2). Integrating its outer product over
propagating and evanescent plane waves gives the normalized specific
impedance. Positive quadrature weights preserve reciprocity and passivity.
All phasors use exp(+iwt); pressure and velocity amplitudes are RMS.
"""
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import j0, j1, jn_zeros


@lru_cache(maxsize=32)
def radial_roots(count):
    if not isinstance(count, int) or isinstance(count, bool) or not 1 <= count <= 64:
        raise ValueError('Radial mode count must be an integer from 1 to 64')
    roots=np.r_[0.,jn_zeros(1,count-1)] if count>1 else np.array([0.])
    roots.flags.writeable=False
    return roots


def mode_spectrum(t, roots):
    """Removable Bessel-root singularities are evaluated by their exact limit."""
    t=np.atleast_1d(np.asarray(t,dtype=float))
    roots=np.asarray(roots,dtype=float)
    denominator=t[None,:]**2-roots[:,None]**2
    close=np.abs(t[None,:]-roots[:,None]) < 1e-7
    result=np.divide(t[None,:]*j1(t)[None,:],denominator,
                     out=np.zeros_like(denominator),where=~close)
    return np.where(close,j0(roots)[:,None]/2,result)


@lru_cache(maxsize=8)
def _evanescent_grid(limit, order):
    # A Bessel-product oscillation has period approximately pi in t.
    edges=np.linspace(0.,limit,int(np.ceil(limit/np.pi))+1)
    nodes,weights=leggauss(order)
    midpoint=(edges[:-1]+edges[1:])/2
    halfwidth=np.diff(edges)/2
    return (midpoint[:,None]+halfwidth[:,None]*nodes).ravel(),(halfwidth[:,None]*weights).ravel()


def modal_radiation_impedance(ka, modes=8, *, spectral_limit=2048., quadrature_order=16):
    """Return modal pressure/velocity impedance divided by rho*c.

    The aperture basis has norm S, so its weak impedance is S*rho*c*Z.
    The finite evanescent spectral limit is explicit and must be checked
    for convergence for the intended band and modal count.
    """
    if not np.isfinite(ka) or not 0 < ka <= 30:
        raise ValueError('This quadrature supports 0 < ka <= 30')
    roots=radial_roots(modes)
    if not np.isfinite(spectral_limit) or spectral_limit < max(128.,4*roots[-1],4*ka):
        raise ValueError('Evanescent spectral limit is too small')
    if not isinstance(quadrature_order,int) or not 8 <= quadrature_order <= 64:
        raise ValueError('Quadrature order must be an integer from 8 to 64')
    nodes,weights=leggauss(max(96,int(np.ceil(8*ka))))
    theta=(nodes+1)*np.pi/4
    propagating=mode_spectrum(ka*np.sin(theta),roots)
    real=(propagating*(2*ka**2*np.sin(theta)*weights*np.pi/4))@propagating.T
    s,weights=_evanescent_grid(float(spectral_limit),quadrature_order)
    evanescent=mode_spectrum(np.sqrt(ka**2+s*s),roots)
    imag=(evanescent*(2*ka*weights))@evanescent.T
    return real+1j*imag


def modal_baffled_on_axis(frequencies, velocities, radius, distance, *, rho=1.225, c=343., order=256):
    """Integrate the actual modal aperture velocity at an on-axis point."""
    f=np.asarray(frequencies,dtype=float)
    velocities=np.asarray(velocities,dtype=complex)
    if f.ndim!=1 or velocities.ndim!=2 or velocities.shape[0]!=len(f):
        raise ValueError('One modal velocity row is required per frequency')
    if not np.isfinite(f).all() or np.any(f<=0) or not np.isfinite(velocities).all():
        raise ValueError('Frequency and modal velocity must be finite')
    if not np.isfinite([radius,distance,rho,c]).all() or min(radius,distance,rho,c)<=0:
        raise ValueError('Positive finite geometry and air properties are required')
    roots=radial_roots(velocities.shape[1])
    nodes,weights=leggauss(order)
    radial=(nodes+1)*radius/2
    basis=j0(roots[:,None]*radial/radius)/j0(roots)[:,None]
    velocity=velocities@basis
    travel=np.sqrt(distance**2+radial**2)
    k=2*np.pi*f/c
    integral=np.sum(velocity*np.exp(-1j*k[:,None]*travel)/travel*(radial*weights*radius/2),axis=1)
    return 1j*rho*c*k*integral
