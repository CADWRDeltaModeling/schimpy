import pytest

from schimpy import nml

def test_namelist_parser():

    file_content = """&NAMELIST1
    ! a full line comment
    VAR1 = 1 ! a line comment
    VAR2 = 2 ! another comment
    / ! a random comment
    &NAMELIST2
    VAR3 = 3
    VAR4 = 4
    /"""
    namelists = nml.parse(file_content)
    print(namelists)

    expected_namelists = {'NAMELIST1': {'comment': "",
            'VAR1': {'value': 1, 'full_line_comment': ' a full line comment', 'inline_comment': ' a line comment'},
            'VAR2': {'value': 2, 'full_line_comment': "", 'inline_comment': ' another comment'}},
        'NAMELIST2': {'comment': "",
            'VAR3': {'value': 3, 'full_line_comment': "", 'inline_comment': ""},
        'VAR4': {'value': 4, 'full_line_comment': "", 'inline_comment': ""}}}

    assert namelists == expected_namelists

def test_write_namelist():
    namelist_data = {
        "first_namelist": {
            "comment": "This is the first namelist",
            "first_value": {
                "value": 1,
                "inline_comment": "This is the first value"
            },
            "second_value": {
                "value": 2
            }
        },
        "second_namelist": {
            "comment": "This is the second namelist",
            "third_value": {
                "value": 3,
                "inline_comment": "This is the third value"
            }
        }
    }
    expected_output = "!This is the first namelist\n&first_namelist\n  first_value = 1 !This is the first value\n  second_value = 2\n/\n!This is the second namelist\n&second_namelist\n  third_value = 3 !This is the third value\n/"
    output = nml.write(namelist_data)
    assert output == expected_output, f"Expected: {expected_output}\nGot: {output}"

    namelist_data = {
        "first_namelist": {
            "first_value": {
                "value": 1
            },
            "second_value": {
                "value": 2
            }
        },
        "second_namelist": {
            "third_value": {
                "value": 3
            }
        }
    }
    expected_output = "&first_namelist\n  first_value = 1\n  second_value = 2\n/\n&second_namelist\n  third_value = 3\n/"
    output = nml.write(namelist_data)
    assert output == expected_output, f"Expected: {expected_output}\nGot: {output}"

@pytest.mark.integration
def test_sample_param_nml(tmp_path):
    url = 'https://raw.githubusercontent.com/schism-dev/schism/master/sample_inputs/param.nml'
    import urllib.request
    with urllib.request.urlopen(url) as response:
        content = response.read().decode("utf-8")
    (tmp_path / 'param-original.nml').write_text(content)
    from schimpy import nml
    params = nml.parse(content)
    first = tmp_path / 'param.nml'
    second = tmp_path / 'param2.nml'
    first.write_text(nml.write(params))
    params2 = nml.parse(first.read_text())
    second.write_text(nml.write(params2))
    import filecmp
    assert filecmp.cmp(first, second), 'Reading and writing from parser to writer is resulting in differences!'


def test_mapping():
    from schimpy import nml
    namelist_dict = {'group1': {'var1': 1, 'var2': 2}, 'group2': {'var3': 3, 'var4': 4}}
    obj = nml.map_to_object(namelist_dict)
    assert obj.group1.var1 == 1
    assert obj.group1.var2 == 2
    assert obj.group2.var3 == 3
    assert obj.group2.var4 == 4
    round_trip_dict = nml.map_to_dict(obj)
    assert namelist_dict == round_trip_dict

