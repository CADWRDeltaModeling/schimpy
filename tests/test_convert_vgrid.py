# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:21:42 2022

@author: zzhang
"""

import pathlib
import pytest
from schimpy.schism_vertical_mesh import convert_vmesh, compare_vmesh, read_vmesh

@pytest.fixture
def datadir():
    return pathlib.Path(__file__).parent.resolve() / 'testdata'


def test_convert_vmesh(datadir):
    old_vgrid_fn = datadir / 'vgrid.in.old'
    new_vgrid_fn = datadir / 'vgrid.in.new'

    output_grid_fn1 = 'vgrid1.in.3d'
    output_grid_fn2 = 'vgrid2.in.3d'
    # converting from old to new vgrid
    convert_vmesh(old_vgrid_fn,
                  output_grid_fn1,
                  input_vgrid_version="5.8",
                  output_vgrid_version="5.10")
    # converting from new to old vgrid.
    convert_vmesh(new_vgrid_fn, output_grid_fn2,
                  input_vgrid_version="5.10",
                  output_vgrid_version="5.8")
    # check if they are still identify with the original grid.
    mesh1 = read_vmesh(old_vgrid_fn, vgrid_version="5.8")
    mesh2 = read_vmesh(output_grid_fn1, vgrid_version="5.10")
    compare_vmesh(mesh1, mesh2)


def test_compare_vmesh(datadir):
    new_vgrid_fn = datadir / 'vgrid.in.new'
    output_grid_fn2 = 'vgrid2.in.3d'
    mesh1 = read_vmesh(new_vgrid_fn, vgrid_version="5.10")
    mesh2 = read_vmesh(output_grid_fn2, vgrid_version="5.8")
    compare_vmesh(mesh1, mesh2)
