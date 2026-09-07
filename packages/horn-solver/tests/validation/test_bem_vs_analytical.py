"""Guard against reintroducing a physically mismatched BEM validation claim.

The old tests compared an open-radiating tube against a matched plane-wave load
and assumed a free-standing mouth converges to an infinite baffled piston.
Those are different exterior problems. The current legacy coupling also maps
all boundary faces, including the inlet and walls, rather than only the mouth.
Until a specified exterior domain is implemented and validated, reject that path.
Standalone Bempp sphere/operator tests remain in test_bem_coupling.py.
"""
import pytest
from horn_solver.solver import run_simulation_from_step


@pytest.mark.validation
def test_unvalidated_bem_cannot_be_used_as_a_reference(tmp_path):
    with pytest.raises(NotImplementedError,match="whole-boundary trace"):
        run_simulation_from_step("unneeded.step",(100,200),3,{"length":.1},
                                 str(tmp_path/"invalid.csv"),200,radiation_model="bem")
    assert not (tmp_path/"invalid.csv").exists()
