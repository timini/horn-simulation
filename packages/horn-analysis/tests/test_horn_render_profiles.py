"""The report render must draw the same curve that is meshed and screened."""
import numpy as np
import pytest

from horn_analysis.horn_render import _radius_profile
from horn_core.profiles import PROFILE_NAMES, get_radius_func

R_T, R_M, L = 0.0254, 0.182, 0.4


@pytest.mark.parametrize("profile", PROFILE_NAMES)
def test_render_profile_is_horn_core_profile(profile):
    z = np.linspace(0, L, 41)
    func = get_radius_func(profile, R_T, R_M, L)
    np.testing.assert_allclose(_radius_profile(z, R_T, R_M, L, profile),
                               [func(float(v)) for v in z], rtol=1e-12)


def test_rendered_tractrix_flares_from_the_throat():
    """The old render held the tractrix near throat radius for most of its length."""
    r = _radius_profile(np.array([L / 4, L / 2]), R_T, R_M, L, "tractrix")
    assert (r[0] - R_T) / (R_M - R_T) > 0.03
    assert (r[1] - R_T) / (R_M - R_T) > 0.15
