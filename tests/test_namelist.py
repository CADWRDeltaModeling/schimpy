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

    from schimpy import namelist
    namelists = namelist.parse(file_content)
    print(namelists)

    expected_namelists = {'NAMELIST1': {'comment': None,
            'VAR1': {'value': '1', 'full_line_comment': 'a full line comment', 'inline_comment': 'a line comment'},
            'VAR2': {'value': '2', 'full_line_comment': None, 'inline_comment': 'another comment'}},
        'NAMELIST2': {'comment': None,
            'VAR3': {'value': '3', 'full_line_comment': None, 'inline_comment': None},
        'VAR4': {'value': '4', 'full_line_comment': None, 'inline_comment': None}}}

    assert namelists == expected_namelists

def test_write_namelist():
    from schimpy import namelist
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
    expected_output = "&first_namelist\n! This is the first namelist\n  first_value = 1\n  ! This is the first value\n  second_value = 2\n/\n&second_namelist\n! This is the second namelist\n  third_value = 3\n  ! This is the third value\n/"
    output = namelist.write(namelist_data)
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
    output = namelist.write(namelist_data)
    assert output == expected_output, f"Expected: {expected_output}\nGot: {output}"

