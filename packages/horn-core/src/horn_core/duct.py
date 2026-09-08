"""Circular-duct losses and radiation, using RMS exp(+i omega t).

References: Berggren et al. (2018), doi:10.1016/j.jcp.2018.06.005;
Silva et al. (2009), arXiv:0811.3625, equations 21–22;
Dalmont et al. (2001), doi:10.1006/jsvi.2000.3487.
Finite-flange expression cross-checked against Ernoult's published benchmark.
These are frequency-domain, plane-mode approximations, not exterior-field solvers.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import jve


@dataclass(frozen=True)
class AirProperties:
    c: float = 343.0
    rho: float = 1.225
    gamma: float = 1.4
    viscosity: float = 1.81e-5
    conductivity: float = .0257
    heat_capacity: float = 1006.0

    def __post_init__(self):
        if not np.isfinite(list(vars(self).values())).all() or min(vars(self).values()) <= 0 or self.gamma <= 1:
            raise ValueError('Air properties must be positive and finite; gamma must exceed one')


DEFAULT_AIR = AirProperties()


def boundary_layer_depths(frequency, air=DEFAULT_AIR):
    f = np.asarray(frequency, dtype=float)
    if np.any(~np.isfinite(f)) or np.any(f <= 0):
        raise ValueError('Positive finite frequencies required')
    omega = 2*np.pi*f
    return (np.sqrt(2*air.viscosity/(air.rho*omega)),
            np.sqrt(2*air.conductivity/(air.rho*air.heat_capacity*omega)))


def circular_duct_properties(frequency, radius, *, air=DEFAULT_AIR, loss_model='lossless'):
    """Return complex wavenumber and characteristic specific impedance.

    ``boundary_layer`` is the thin-layer limit of the FEM wall condition;
    ``kirchhoff`` retains exact circular cross-section viscothermal functions.
    A local circular slice of a flaring horn remains a 1D approximation.
    """
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError('Positive finite duct radius required')
    dv, dt = boundary_layer_depths(frequency, air)
    omega = 2*np.pi*np.asarray(frequency)
    if loss_model == 'lossless':
        return omega/air.c + 0j, np.full_like(omega, air.rho*air.c, dtype=complex)
    if loss_model == 'boundary_layer':
        if max(np.max(dv),np.max(dt))/radius > .1:
            raise ValueError('Thin boundary-layer model requires layer depth / radius <= 0.1')
        fv, ft = (1-1j)*dv/radius, (1-1j)*dt/radius
    elif loss_model == 'kirchhoff':
        xv, xt = (1-1j)*radius/dv, (1-1j)*radius/dt
        fv = 2*jve(1,xv)/(xv*jve(0,xv))
        ft = 2*jve(1,xt)/(xt*jve(0,xt))
    else:
        raise ValueError('Unknown loss model')
    density = air.rho/(1-fv)
    compressibility = (1+(air.gamma-1)*ft)/(air.rho*air.c**2)
    return omega*np.sqrt(density*compressibility), np.sqrt(density/compressibility)


def _silva_reflection(ka, flanged):
    if flanged:
        a1,a2,a3,b1,b2,b3,b4,beta,eta = .73,.372,.0231,.244,.723,-.0198,.00366,1.,.8216
    else:
        a1,a2,a3,b1,b2,b3,b4,beta,eta = .8,.266,.0263,.0599,.238,-.0153,.0015,.5,.6133
    x = ka*ka
    modulus = (1+a1*x)/(1+(beta+a1)*x+a2*x*x+a3*x*x*x)
    correction = eta*(1+b1*x)/(1+b2*x+b3*x*x+b4*x*x*x)
    return modulus, correction


class RadiationDomainError(ValueError):
    """A valid geometry lies outside the selected radiation approximation."""


def circular_pipe_radiation(k, radius, flange_width=0.):
    """Normalized specific radiation impedance for a finite circular flange.

    Zero width gives the unflanged Silva approximation. Restricted to ka < 1.5
    for this project's declared plane-mode domain. No measured coefficients fit.
    """
    k = np.asarray(k,dtype=float)
    if not np.isfinite([radius,flange_width]).all() or radius <= 0 or flange_width < 0 or np.any(~np.isfinite(k)) or np.any(k<=0):
        raise ValueError('Invalid circular radiation dimensions or wavenumber')
    ka = k*radius
    if np.any(ka >= 1.5):
        raise RadiationDomainError('Circular pipe radiation requires ka < 1.5')
    if flange_width > radius*(1+1e-10):
        raise RadiationDomainError('Finite-flange approximation restricted to width / radius <= 1')
    r0,l0 = _silva_reflection(ka,False)
    if flange_width == 0:
        reflection = -r0*np.exp(-2j*ka*l0)
    else:
        ri,li = _silva_reflection(ka,True)
        b = radius+flange_width
        ratio = radius/b
        # Complex end corrections combine length correction and reflection loss.
        d0,di = l0+.5j*np.log(r0)/ka, li+.5j*np.log(ri)/ka
        correction = di+ratio*(d0-di)+.057*ratio*(1-ratio**5)
        reflection = -np.exp(-2j*ka*correction)
        reflection -= (.43*ratio*(1-ratio)*np.sin(k*b/(1.85-ratio))**2
                       * np.exp(-1j*k*b*(1+ratio*(2.3-ratio-.3*ka**2))))
    z = (1+reflection)/(1-reflection)
    if np.any(z.real < -1e-10):
        raise RadiationDomainError('Radiation approximation left its passive domain')
    return z
