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
