"""The external comparison must preserve absolute levels and complex phase."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('mmm_check',Path(__file__).resolve().parents[1]/'validate_mmm_radiation.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)


def test_absolute_output_gain_cannot_be_normalized_away():
    result=v.change(np.array([1.,2.])*1j,np.array([2.,4.])*1j)
    assert result['level_db']==pytest.approx(20*np.log10(2.))
    assert result['phase_deg']==0.


def test_opposite_phasor_convention_does_not_silently_match():
    result=v.change(np.array([1j]),np.array([-1j]))
    assert result['phase_deg']==180.


@pytest.mark.parametrize('bad',[0.,np.inf,np.nan])
def test_invalid_reference_cannot_pass(bad):
    with pytest.raises(ValueError,match='Invalid'):
        v.change(np.array([bad]),np.ones(1))


def test_fixed_specification_cannot_change_after_preparation(tmp_path,monkeypatch):
    import json
    monkeypatch.setattr(v,'clean_revision',lambda:'test')
    out=tmp_path/'study';v.prepare(out)
    v.verify(out)
    protocol=json.loads((out/'protocol.json').read_text())
    protocol['geometry']['length']*=2
    (out/'protocol.json').write_text(json.dumps(protocol))
    with pytest.raises(ValueError,match='physical specification'):
        v.verify(out)
