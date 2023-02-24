# tests for params.py 
# Param NML
import pytest
from schimpy import params, namelist

def test_params():
    param_file = 'testdata/param.nml'
    with open(param_file,'r') as fh: param_dict = namelist.parse(fh.read())
    assert 'CORE' in param_dict
    p = params.create_params(param_dict)
    assert p.core.ipre == 0
    p.core.ipre = 1
    assert p.core.ipre == 1
    with pytest.raises(Exception) as ex:
      p.core.ipre=2
    assert p.core.ipre == 1
