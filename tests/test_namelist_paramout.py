from pathlib import Path
from schimpy import nml


def test_read_param_out():
    dirpath = Path(__file__).parent
    filepath = dirpath / "testdata/param.out.nml"
    with open(filepath, "r") as fin:
        nmldata = nml.parse(fin.read())
    assert nmldata["SCHOUT"]["IOF_HYDRO"]["value"] == "1, 16*0, 2*1, 5*0, 2*1, 14*0,"
