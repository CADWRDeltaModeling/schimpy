# tests for params.py
# Param NML
import pathlib
import pytest
from schimpy import param, params, nml
import pandas as pd


@pytest.fixture
def datadir():
    return pathlib.Path(__file__).parent.resolve() / "testdata"


def test_params(datadir):
    param_file = datadir / "param.nml"
    with open(param_file, "r") as fh:
        param_dict = nml.parse(fh.read())
    assert "CORE" in param_dict
    p = params.create_params(param_dict)
    assert p.core.ipre == 0
    p.core.ipre = 1
    assert p.core.ipre == 1
    with pytest.raises(Exception) as ex:
        p.core.ipre = 2
    assert p.core.ipre == 1


def test_param_special(datadir):
    """tests properties that do translation of times and intervals"""
    test_param_file = datadir / "param.nml"
    parms = param.read_params(test_param_file)
    parms["rnday"] = 3000
    parms["dt"] = 90.0
    assert parms["rnday"] == 3000
    parms.run_start = "2010-02-04"
    assert parms.run_start == pd.Timestamp(2010, 2, 4)
    parms.nc_stack = "8H"
    assert parms.nc_stack == pd.tseries.frequencies.to_offset("8H")
    parms.nc_stack = pd.tseries.frequencies.to_offset("1D")
    assert parms["ihfskip"] == 86400 / 90
    parms.hotstart_freq = None
    assert parms["nhot"] == 0
    parms.hotstart_freq = "5D"
    assert parms["nhot"] == 1
    assert parms["nhot_write"] == (5 * 86400 / parms["dt"])

    parms.station_out_freq = "30T"
    assert parms["iout_sta"] == 1
    assert parms["nspool_sta"] == 20
    parms.station_out_freq = None
    assert parms["iout_sta"] == 0
    parms.station_out_freq = "15T"

    print(parms.hotstart_freq)
    print(parms["nhot"])
    parms.hotstart_freq = "10D"
    assert "CORE" in parms.sections()

    dupl = param.read_params(test_param_file)
    differences = not parms.diff(dupl).empty
    assert differences
    dupl2 = param.read_params(test_param_file)
    no_differences = dupl.diff(dupl2).empty
    assert no_differences
    assert parms.validate() is None
