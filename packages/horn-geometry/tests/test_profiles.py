"""Tests for exponential and hyperbolic horn geometry generators."""

import pytest
import numpy as np
from pathlib import Path

try:
    import gmsh
except ImportError:
    gmsh = None

from horn_geometry.generator import (
    create_conical_horn,
    create_exponential_horn,
    create_hyperbolic_horn,
    create_horn,
)


@pytest.mark.skipif(gmsh is None, reason="gmsh not available")
class TestProfileSmoke:
    """Smoke tests: each profile creates a non-empty STEP file."""

    @pytest.fixture(params=["conical", "exponential", "hyperbolic"])
    def profile_name(self, request):
        return request.param

    def test_create_horn_produces_step_file(self, profile_name, tmp_path):
        output = tmp_path / f"{profile_name}_horn.step"
        result = create_horn(
            profile=profile_name,
            throat_radius=0.025,
            mouth_radius=0.1,
            length=0.3,
            output_file=output,
            num_sections=10,
        )
        assert result.exists()
        assert result.stat().st_size > 0


@pytest.mark.skipif(gmsh is None, reason="gmsh not available")
class TestProfileGeometricValidation:
    """Validate geometric properties against analytical formulae."""

    THROAT_R = 0.025
    MOUTH_R = 0.1
    LENGTH = 0.3

    def _analytical_volume_conical(self, r_t, r_m, L):
        """V = (pi/3) * L * (R^2 + R*r + r^2)"""
        return (np.pi / 3) * L * (r_m**2 + r_m * r_t + r_t**2)

    def _analytical_volume_exponential(self, r_t, r_m, L):
        """V = integral_0^L pi * [r_t * exp(m*z/L)]^2 dz
        = pi * r_t^2 * L / (2m) * [exp(2m) - 1]
        where m = ln(r_m / r_t).
        """
        m = np.log(r_m / r_t)
        return np.pi * r_t**2 * L / (2 * m) * (np.exp(2 * m) - 1)

    def _analytical_volume_hyperbolic(self, r_t, r_m, L):
        """V = integral_0^L pi * [r_t * cosh(m*z/L)]^2 dz
        = pi * r_t^2 * L * [1/2 + sinh(2m) / (4m)]
        where m = acosh(r_m / r_t).
        """
        m = np.arccosh(r_m / r_t)
        return np.pi * r_t**2 * L * (0.5 + np.sinh(2 * m) / (4 * m))

    def _get_step_volume(self, step_file: Path) -> float:
        """Load a STEP file with gmsh and compute the volume of all 3D entities."""
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("volume_check")
        gmsh.model.occ.importShapes(str(step_file))
        gmsh.model.occ.synchronize()

        volumes = gmsh.model.occ.getEntities(dim=3)
        total_volume = 0.0
        for dim, tag in volumes:
            mass = gmsh.model.occ.getMass(dim, tag)
            total_volume += mass

        gmsh.finalize()
        return total_volume

    def test_conical_volume(self, tmp_path):
        step = create_conical_horn(
            self.THROAT_R, self.MOUTH_R, self.LENGTH,
            tmp_path / "conical.step", num_sections=2,
        )
        expected = self._analytical_volume_conical(self.THROAT_R, self.MOUTH_R, self.LENGTH)
        actual = self._get_step_volume(step)
        assert np.isclose(actual, expected, rtol=0.01), (
            f"Conical volume mismatch: expected={expected:.6e}, actual={actual:.6e}"
        )

    def test_exponential_volume(self, tmp_path):
        step = create_exponential_horn(
            self.THROAT_R, self.MOUTH_R, self.LENGTH,
            tmp_path / "exponential.step", num_sections=40,
        )
        expected = self._analytical_volume_exponential(self.THROAT_R, self.MOUTH_R, self.LENGTH)
        actual = self._get_step_volume(step)
        assert np.isclose(actual, expected, rtol=0.01), (
            f"Exponential volume mismatch: expected={expected:.6e}, actual={actual:.6e}"
        )

    def test_hyperbolic_volume(self, tmp_path):
        step = create_hyperbolic_horn(
            self.THROAT_R, self.MOUTH_R, self.LENGTH,
            tmp_path / "hyperbolic.step", num_sections=40,
        )
        expected = self._analytical_volume_hyperbolic(self.THROAT_R, self.MOUTH_R, self.LENGTH)
        actual = self._get_step_volume(step)
        assert np.isclose(actual, expected, rtol=0.01), (
            f"Hyperbolic volume mismatch: expected={expected:.6e}, actual={actual:.6e}"
        )

    def test_throat_mouth_radii_match(self, tmp_path):
        """Verify that throat (z=0) and mouth (z=L) surfaces have correct areas."""
        for profile in ["exponential", "hyperbolic"]:
            step = create_horn(
                profile=profile,
                throat_radius=self.THROAT_R,
                mouth_radius=self.MOUTH_R,
                length=self.LENGTH,
                output_file=tmp_path / f"{profile}_area_check.step",
                num_sections=20,
            )
            gmsh.initialize()
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("area_check")
            gmsh.model.occ.importShapes(str(step))
            gmsh.model.occ.synchronize()

            surfaces = gmsh.model.occ.getEntities(dim=2)
            throat_area = None
            mouth_area = None
            for dim, tag in surfaces:
                com = gmsh.model.occ.getCenterOfMass(dim, tag)
                if np.isclose(com[2], 0.0, atol=1e-6):
                    throat_area = gmsh.model.occ.getMass(dim, tag)
                elif np.isclose(com[2], self.LENGTH, atol=1e-6):
                    mouth_area = gmsh.model.occ.getMass(dim, tag)

            gmsh.finalize()

            expected_throat = np.pi * self.THROAT_R**2
            expected_mouth = np.pi * self.MOUTH_R**2
            assert throat_area is not None, f"{profile}: no throat face found at z=0"
            assert mouth_area is not None, f"{profile}: no mouth face found at z=L"
            assert np.isclose(throat_area, expected_throat, rtol=0.01), (
                f"{profile} throat area: expected={expected_throat:.6e}, actual={throat_area:.6e}"
            )
            assert np.isclose(mouth_area, expected_mouth, rtol=0.01), (
                f"{profile} mouth area: expected={expected_mouth:.6e}, actual={mouth_area:.6e}"
            )


@pytest.mark.skipif(gmsh is None, reason="gmsh not available")
class TestCreateHornDispatch:
    """Test the create_horn dispatch function."""

    def test_invalid_profile_raises(self, tmp_path):
        with pytest.raises(ValueError, match="Unknown profile"):
            create_horn(
                profile="parabolic",
                throat_radius=0.025,
                mouth_radius=0.1,
                length=0.3,
                output_file=tmp_path / "bad.step",
            )

    def test_conical_backward_compatible(self, tmp_path):
        """create_horn('conical', ...) should produce the same result as create_conical_horn()."""
        output = tmp_path / "dispatch.step"
        result = create_horn(
            profile="conical",
            throat_radius=0.025,
            mouth_radius=0.1,
            length=0.3,
            output_file=output,
            num_sections=2,
        )
        assert result.exists()
        assert result.stat().st_size > 0
