import numpy as np
import pandas as pd
import pytest
from horn_analysis.single_report import _diagnostics_section, generate_single_report


@pytest.mark.parametrize('column', ['relative_residual', 'residual'])
@pytest.mark.parametrize('bad', [float('inf'), float('nan'), -1, 1e-3])
def test_bad_residual_cannot_disappear_from_report(column, bad):
    result = _diagnostics_section(pd.DataFrame({'frequency':[100,200],column:[1e-14,bad],'converged_reason':[4,4]}))
    assert '1 of 2 frequencies' in result and 'unreliable' in result
    assert 'All reported solver checks passed' not in result


def test_absent_diagnostics_are_not_a_pass():
    assert 'not available' in _diagnostics_section(pd.DataFrame({'frequency':[100,200]}))
    assert 'incomplete' in _diagnostics_section(pd.DataFrame({'residual':[1e-14]}))


def test_negative_convergence_reason_is_counted_once():
    result = _diagnostics_section(pd.DataFrame({'relative_residual':[1e-3,1e-14],'converged_reason':[-3,4]}))
    assert '1 of 2 frequencies' in result


def test_generated_report_includes_actual_solver_diagnostics(tmp_path):
    png = tmp_path/'plot.png';png.write_bytes(b'fixture')
    csv = tmp_path/'result.csv'
    pd.DataFrame({'frequency':[100,200],'relative_residual':[1e-14,np.inf],'converged_reason':[4,4]}).to_csv(csv,index=False)
    result = generate_single_report(.01,.05,.1,'conical',{},str(png),str(png),str(png),str(png),str(png),final_csv=str(csv))
    assert '<h2>Solver checks</h2>' in result and '1 of 2 frequencies' in result
