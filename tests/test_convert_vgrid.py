# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:21:42 2022

@author: zzhang
"""

from schimpy.schism_vertical_mesh import convert_vmesh, compare_vmesh, read_vmesh

old_vgrid_fn = 'testdata/vgrid.in.old'
new_vgrid_fn = 'testdata/vgrid.in.new'

output_grid_fn1 = 'vgrid1.in.3d'
output_grid_fn2 = 'vgrid2.in.3d'

convert_vmesh(old_vgrid_fn,output_grid_fn1) #converting from old to new vgrid
convert_vmesh(new_vgrid_fn,output_grid_fn2,input_vgrid=5.9,output_vgrid=5.8) #converting from new to old vgrid.

# check if they are still identify with the original grid. 
mesh1 = read_vmesh(old_vgrid_fn,vgrid_version=5.8)
mesh2 = read_vmesh(output_grid_fn1,vgrid_version=5.9)
compare_vmesh(mesh1,mesh2)

mesh1 = read_vmesh(new_vgrid_fn,vgrid_version=5.9)
mesh2 = read_vmesh(output_grid_fn2,vgrid_version=5.8)
compare_vmesh(mesh1,mesh2)

