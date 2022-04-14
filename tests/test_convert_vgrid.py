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
convert_vmesh(new_vgrid_fn,output_grid_fn2,input_vgrid='new',output_vgrid='old') #converting from new to old vgrid.

# check if they are still identify with the original grid. 
mesh1 = read_vmesh(old_vgrid_fn)
mesh2 = read_vmesh(output_grid_fn1,old_vgrid=False)
compare_vmesh(mesh1,mesh2)

mesh1 = read_vmesh(new_vgrid_fn,old_vgrid=False)
mesh2 = read_vmesh(output_grid_fn2)
compare_vmesh(mesh1,mesh2)