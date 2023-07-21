from schimpy import nml

def test_read_param_out():
    with open("testdata/param.out.nml","r") as fin:
        nmldata = nml.parse(fin.read())
    assert nmldata['SCHOUT']['IOF_HYDRO']['value'] == '1, 16*0, 2*1, 5*0, 2*1, 14*0,'