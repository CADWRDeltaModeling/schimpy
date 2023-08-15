#!/usr/bin/env python
# -*- coding: utf-8 -*-
import schimpy.schism_linestring as schism_linestring
import os
import pytest


@pytest.fixture
def dirpath_test():
    dir_cur = os.path.dirname(__file__)
    return os.path.join(dir_cur, "testdata/schism_linestring")


def test_read_linestrings_yaml(dirpath_test, tmp_path):
    fpath = os.path.join(dirpath_test, "flowlines.yaml")
    lines = schism_linestring.read_linestrings(fpath)
    assert [l.prop["name"] for l in lines] == ["ocean", "mixed1", "middle"]
    assert [list(l.coords) for l in lines] == [
        [(41.0, 69.0), (41.0, 101.0)],
        [(101.0, 75.0), (69.0, 75.0)],
        [(45.0, 5.0), (45.0, 25.0)],
    ]


def test_read_linestrings_shapefile(dirpath_test):
    fpath = os.path.join(dirpath_test, "test_linestrings.shp")
    lines = schism_linestring.read_linestrings(fpath)
    assert [l.prop["name"] for l in lines] == ["ocean", "mixed1", "middle"]
    assert [list(l.coords) for l in lines] == [
        [(41.0, 69.0), (41.0, 101.0)],
        [(101.0, 75.0), (69.0, 75.0)],
        [(45.0, 5.0), (45.0, 25.0)],
    ]


def test_write_linestrings_yaml(dirpath_test, tmp_path):
    fpath_in = os.path.join(dirpath_test, "flowlines.yaml")
    lines = schism_linestring.read_linestrings(fpath_in)
    fpath_out = tmp_path / "test_linestrings.yaml"
    schism_linestring.write_linestrings(fpath_out, lines)
    lines_readback = schism_linestring.read_linestrings(fpath_out)
    assert [l.prop["name"] for l in lines] == [l.prop["name"] for l in lines_readback]
    assert [list(l.coords) for l in lines] == [list(l.coords) for l in lines_readback]


def test_write_linestrings_shp(dirpath_test, tmp_path):
    fpath_in = os.path.join(dirpath_test, "flowlines.yaml")
    lines = schism_linestring.read_linestrings(fpath_in)
    fpath_out = tmp_path / "test_linestrings.shp"
    schism_linestring.write_linestrings(fpath_out, lines)
