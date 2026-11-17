"""Pure-Python radius functions for the horn flare profiles.

This module is the single implementation of every profile. The STEP
generator (horn_geometry.generator), the report renderer
(horn_analysis.horn_render) and the TMM/LEM screen all call it, so the
geometry that is screened, meshed and drawn is the same curve.

Only numpy and scipy are required. Each radius function maps axial
position z in [0, L] to the local horn radius r(z), all in metres.

Profiles
--------
conical          r = r_t + (r_m - r_t) z / L
exponential      r = r_t exp(z / x0)                      (hypex T = 1)
hyperbolic       Salmon hypex, default T = 0 (catenoidal)  r = r_t (cosh(z/x0) + T sinh(z/x0))
tractrix         true tractrix, constant tangent length a >= r_m (see tractrix_radius)
os               oblate spheroidal (Geddes) r = sqrt(r_t^2 + z^2 tan^2(theta))
scaled_tractrix  tractrix contour with a = r_m, axially scaled to L (formerly "lecleach")
cd               exponential throat joined C1 to a tangent cone

"lecleach" is a deprecated alias for "scaled_tractrix". The old
implementation was never Le Cleac'h's horn, which has no closed-form
contour (it is computed step by step from a wavefront model; see
B. Kolbrek, "Horn Theory: An Introduction, Part 2", audioXpress 2008).
"""

import math
import warnings
from typing import Callable, Dict, Optional

import numpy as np
from scipy.optimize import brentq


class ProfileGeometryError(ValueError):
    """The requested (r_t, r_m, L) cannot be realised by this profile."""


# -- Simple closed forms ----------------------------------------------------

def conical_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Linear flare: r(z) = r_t + (r_m - r_t) * z / L."""
    return r_t + (r_m - r_t) * z / L


def exponential_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Exponential flare: r(z) = r_t * exp(m * z / L), m = ln(r_m / r_t)."""
    m = math.log(r_m / r_t)
    return r_t * math.exp(m * z / L)


def hypex_flare_rate(r_t: float, r_m: float, L: float, T: float) -> float:
    """Return u = L / x0 so that r_t (cosh u + T sinh u) = r_m.

    With w = exp(u) the condition is (1 + T) w^2 - 2 R w + (1 - T) = 0,
    R = r_m / r_t, whose larger root is the physical one.
    """
    if T < 0:
        raise ValueError("hypex T must be >= 0")
    R = r_m / r_t
    if T == 0:
        return float(np.arccosh(R))
    w = (R + math.sqrt(R * R - (1.0 - T * T))) / (1.0 + T)
    return math.log(w)


def hypex_radius(z: float, r_t: float, r_m: float, L: float, T: float = 0.0) -> float:
    """Salmon hypex flare: r(z) = r_t (cosh(z/x0) + T sinh(z/x0)).

    The area S = pi r^2 follows S_t (cosh(z/x0) + T sinh(z/x0))^2
    (V. Salmon, "A New Family of Horns", JASA 17(3), 1946). x0 is chosen
    so that r(L) = r_m. T = 1 is the exponential horn, T = 0 the
    catenoidal horn, and 0 < T < 1 the hypex horns in between.
    """
    u = hypex_flare_rate(r_t, r_m, L, T)
    s = u * z / L
    if T == 0:
        return float(r_t * np.cosh(s))
    return float(r_t * (math.cosh(s) + T * math.sinh(s)))


DEFAULT_HYPEX_T = 0.0


def hyperbolic_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """The "hyperbolic" profile: hypex with T = 0, i.e. the catenoidal horn.

    r(z) = r_t * cosh(m * z / L), m = acosh(r_m / r_t). Its wall is
    parallel to the axis at the throat. Use hypex_radius for other T.
    """
    return hypex_radius(z, r_t, r_m, L, DEFAULT_HYPEX_T)


def os_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Oblate spheroidal (OS / Geddes) flare: r(z) = sqrt(r_t^2 + (z tan(theta))^2).

    The wall starts parallel to the axis and approaches the cone of
    half-angle theta = atan(sqrt(r_m^2 - r_t^2) / L) as z grows.
    """
    theta = math.atan2(math.sqrt(r_m**2 - r_t**2), L)
    return math.sqrt(r_t**2 + (z * math.tan(theta)) ** 2)


# -- True tractrix ------------------------------------------------------------

def _tractrix_g(t):
    return t - np.tanh(t)


def _tractrix_span(a: float, r_t: float, r_m: float) -> float:
    """Axial distance between radius r_t and radius r_m on a tractrix of constant a."""
    t_t = math.acosh(a / r_t)
    t_m = math.acosh(a / r_m) if a > r_m else 0.0
    return a * (float(_tractrix_g(t_t)) - float(_tractrix_g(t_m)))


def natural_tractrix_length(r_t: float, r_m: float) -> float:
    """Length of the full tractrix horn with a = r_m (wall at 90 deg at the mouth).

    L_nat = a ln((a + sqrt(a^2 - r_t^2)) / r_t) - sqrt(a^2 - r_t^2), a = r_m.
    """
    a = r_m
    root = math.sqrt(a * a - r_t * r_t)
    return a * math.log((a + root) / r_t) - root


_LENGTH_RTOL = 1e-6


def tractrix_constant(r_t: float, r_m: float, L: float) -> float:
    """Tractrix constant a (tangent length) for a horn from r_t to r_m in length L.

    The classical tractrix horn has a = r_m and its natural length
    (natural_tractrix_length). A longer requested length is met by a
    larger a: the same tractrix, truncated at r_m before its wall turns
    perpendicular to the axis (a lower flare cutoff c / (2 pi a)). A
    shorter length has no tractrix through both radii and raises
    ProfileGeometryError.
    """
    if not (0 < r_t < r_m) or L <= 0:
        raise ProfileGeometryError("tractrix needs 0 < r_t < r_m and L > 0")
    L_nat = natural_tractrix_length(r_t, r_m)
    if L < L_nat * (1 - _LENGTH_RTOL):
        raise ProfileGeometryError(
            f"tractrix from r_t={r_t:.6g} m to r_m={r_m:.6g} m cannot be shorter than its "
            f"natural length {L_nat:.6g} m (requested {L:.6g} m)")
    if L <= L_nat * (1 + _LENGTH_RTOL):
        return r_m
    hi = 2.0 * r_m
    while _tractrix_span(hi, r_t, r_m) < L:
        hi *= 2.0
    return float(brentq(lambda a: _tractrix_span(a, r_t, r_m) - L, r_m, hi, xtol=1e-15, rtol=1e-13))


def _tractrix_func(r_t: float, r_m: float, L: float) -> Callable[[float], float]:
    a = tractrix_constant(r_t, r_m, L)
    t_t = math.acosh(a / r_t)
    t_m = math.acosh(a / r_m) if a > r_m else 0.0
    g_t = float(_tractrix_g(t_t))

    def radius(z: float) -> float:
        z = min(max(float(z), 0.0), L)
        target = g_t - z / a
        lo, hi = t_m, t_t
        # g is increasing in t; bisection is robust where dg/dt -> 0 at the mouth.
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if mid - math.tanh(mid) > target:
                hi = mid
            else:
                lo = mid
        return float(a / math.cosh(0.5 * (lo + hi)))

    return radius


def tractrix_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """True tractrix: every wall tangent meets the axis a distance a from the wall.

    With t >= 0 as parameter, r = a / cosh(t) and the axial position
    measured from the a-radius end is a (t - tanh t). The horn runs from
    r_t (z = 0) to r_m (z = L). See tractrix_constant for how L sets a.
    """
    return _tractrix_func(r_t, r_m, L)(z)


# -- Axially scaled tractrix (formerly "lecleach") ---------------------------

def _build_scaled_tractrix_interp(r_t: float, r_m: float, L: float):
    ratio = r_t / r_m
    if ratio == 1:
        return np.array([0., L]), np.array([r_t, r_m])
    t = np.linspace(np.pi - np.arcsin(ratio), np.pi / 2, 500)
    y, x = np.sin(t), np.log(np.tan(t / 2)) + np.cos(t)
    y[0] = ratio
    x -= x[0]
    return x / x[-1] * L, y * r_m


def scaled_tractrix_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Tractrix contour with a = r_m from r_t to r_m, stretched axially to length L.

    The radii are those of the a = r_m tractrix, so the wall is parallel
    to the axis near the throat and perpendicular at the mouth. The axial
    coordinate is scaled linearly to the requested L, so the constant
    tangent-length property holds only when L equals
    natural_tractrix_length(r_t, r_m). This profile was previously
    (mis)named "lecleach".
    """
    x_scaled, y_scaled = _build_scaled_tractrix_interp(r_t, r_m, L)
    return float(np.interp(z, x_scaled, y_scaled))


def lecleach_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Deprecated alias of scaled_tractrix_radius (it never was a Le Cleac'h horn)."""
    warnings.warn("'lecleach' is deprecated; it is the scaled_tractrix profile",
                  DeprecationWarning, stacklevel=2)
    return scaled_tractrix_radius(z, r_t, r_m, L)


# -- CD: exponential throat with a tangent cone -------------------------------

CD_THROAT_FRACTION = 0.3


def cd_flare_rate(r_t: float, r_m: float, L: float, frac: float = CD_THROAT_FRACTION) -> float:
    """Exponential rate m (1/m) that makes the CD join slope-continuous.

    The exponential r_t exp(m z) runs to z_j = frac * L and the cone
    continues along its tangent, so r_t e^{m z_j} (1 + m (L - z_j)) = r_m.
    The left side increases with m, so the root is unique.
    """
    z_j = frac * L

    def mouth_error(m):
        return r_t * math.exp(m * z_j) * (1.0 + m * (L - z_j)) - r_m

    hi = 1.0 / L
    while mouth_error(hi) < 0:
        hi *= 2.0
    return float(brentq(mouth_error, 0.0, hi, xtol=1e-14, rtol=1e-14))


def cd_radius(z: float, r_t: float, r_m: float, L: float) -> float:
    """Exponential throat (first 30% of L) joined with matching slope to a cone.

    The cone is tangent to the exponential at the join, so r and dr/dz are
    both continuous. This is a simple compound flare for screening; it is
    not a Keele constant-directivity horn (no separate horizontal and
    vertical coverage angles or mouth flare).
    """
    m = cd_flare_rate(r_t, r_m, L)
    z_j = CD_THROAT_FRACTION * L
    if z <= z_j:
        return r_t * math.exp(m * z)
    r_j = r_t * math.exp(m * z_j)
    return r_j + r_j * m * (z - z_j)


# -- Dispatch -------------------------------------------------------------------

PROFILE_NAMES = ["conical", "exponential", "hyperbolic", "tractrix", "os", "scaled_tractrix", "cd"]

DEPRECATED_PROFILE_ALIASES: Dict[str, str] = {"lecleach": "scaled_tractrix"}

PROFILE_DESCRIPTIONS = {
    "conical": "Conical",
    "exponential": "Exponential (hypex T=1)",
    "hyperbolic": "Hyperbolic: hypex T=0 (catenoidal)",
    "tractrix": "Tractrix",
    "os": "Oblate spheroidal",
    "scaled_tractrix": "Axially scaled tractrix",
    "cd": "Exponential throat with tangent cone",
}

_PROFILE_DISPATCH = {
    "conical": conical_radius,
    "exponential": exponential_radius,
    "hyperbolic": hyperbolic_radius,
    "tractrix": tractrix_radius,
    "os": os_radius,
    "scaled_tractrix": scaled_tractrix_radius,
    "cd": cd_radius,
}


def canonical_profile(profile: str) -> str:
    """Return the canonical profile name, warning about deprecated aliases."""
    if profile in DEPRECATED_PROFILE_ALIASES:
        new = DEPRECATED_PROFILE_ALIASES[profile]
        warnings.warn(f"Profile '{profile}' is deprecated; use '{new}' (same geometry)",
                      DeprecationWarning, stacklevel=3)
        return new
    if profile not in _PROFILE_DISPATCH:
        raise ValueError(
            f"Unknown profile '{profile}'. Choose from: {PROFILE_NAMES + list(DEPRECATED_PROFILE_ALIASES)}"
        )
    return profile


def profile_geometry_error(profile: str, r_t: float, r_m: float, L: float) -> Optional[str]:
    """Return why (profile, r_t, r_m, L) is not realisable, or None if it is."""
    try:
        name = DEPRECATED_PROFILE_ALIASES.get(profile, profile)
        if name == "tractrix":
            tractrix_constant(r_t, r_m, L)
    except ProfileGeometryError as error:
        return str(error)
    return None


def get_radius_func(
    profile: str, r_t: float, r_m: float, L: float, *, hypex_t: Optional[float] = None,
) -> Callable[[float], float]:
    """Return a closure r(z) for the given profile and geometry.

    Args:
        profile: One of PROFILE_NAMES (or the deprecated alias "lecleach").
        r_t: Throat radius (m).
        r_m: Mouth radius (m).
        L: Horn length (m).
        hypex_t: Hypex T for the "hyperbolic" profile (default 0, catenoidal).

    Returns:
        Callable that maps z in [0, L] to radius (m).

    Raises:
        ProfileGeometryError: the profile cannot pass through both radii in
            length L (a tractrix shorter than its natural length).
    """
    profile = canonical_profile(profile)
    if hypex_t is not None and profile != "hyperbolic":
        raise ValueError("hypex_t applies only to the 'hyperbolic' profile")

    if profile == "tractrix":
        return _tractrix_func(r_t, r_m, L)

    if profile == "scaled_tractrix":
        x_scaled, y_scaled = _build_scaled_tractrix_interp(r_t, r_m, L)
        return lambda z: float(np.interp(z, x_scaled, y_scaled))

    if profile == "hyperbolic":
        T = DEFAULT_HYPEX_T if hypex_t is None else float(hypex_t)
        return lambda z: hypex_radius(z, r_t, r_m, L, T)

    if profile == "cd":
        m = cd_flare_rate(r_t, r_m, L)
        z_j = CD_THROAT_FRACTION * L
        r_j = r_t * math.exp(m * z_j)
        return lambda z: (r_t * math.exp(m * z) if z <= z_j else r_j + r_j * m * (z - z_j))

    raw_func = _PROFILE_DISPATCH[profile]
    return lambda z: raw_func(z, r_t, r_m, L)


def radius_profile(z, r_t: float, r_m: float, L: float, profile: str, **kwargs) -> np.ndarray:
    """Vectorised r(z) for an array of axial positions."""
    func = get_radius_func(profile, r_t, r_m, L, **kwargs)
    return np.array([func(float(zi)) for zi in np.atleast_1d(np.asarray(z, dtype=float))])
