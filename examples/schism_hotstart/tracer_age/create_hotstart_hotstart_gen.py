# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 11:12:19 2020

@author: zzhang
"""

import schimpy.schism_hotstart as sh
import matplotlib.pyplot as plt
import xarray as xr

yaml_fn = "./data_in/hotstart_from_prev_hotstart.yaml"
modules = ['TEM','SAL','GEN','AGE']
hotstart_fn = "hotstart2_hotstart_tracer_plg.nc" # output hotstart file
ini_date = '2021-04-20' 

#%% first convert mesh2 into utm coordinate system. 
# from schimpy.schism_mesh import read_mesh, write_mesh
# from schimpy.geo_tools import project_fun
# import numpy as np
# mesh1 =  read_mesh('./data_in/hgrid.gr3',
#                    proj4='EPSG:32610')   # utm coordinates
# mesh2 =  read_mesh('./data_in/hgrid2.gr3',
#                    proj4='EPSG:4326 ')   # lat, lon coordinates
# mesh2_new = sh.project_mesh(mesh2,'EPSG:32610')
# write_mesh(mesh2_new,'./data_in/hgrid2_new.gr3')

#%% create a hotstart file for SCHISM
h = sh.hotstart(yaml_fn,modules=modules,param_nml='./data_in/param.nml.clinic',
                proj4 ='EPSG:32610')
h.create_hotstart()
hnc = h.nc_dataset
  
#%% merge the two hotstart file
nc1 = xr.open_dataset("hotstart.48000.20210609.nc")
var_list = ['tr_nd','tr_nd0','tr_el']
nc1 = nc1.drop_vars(var_list)
hnc_new = xr.merge([nc1, hnc[var_list]])
hnc_new.to_netcdf(hotstart_fn) 

#%% making a 2D surface plot
coll = h.mesh.plot_elems(hnc['tr_el'].values[:,0,-2]) #clim=[0,35])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('GEN_1')
plt.tight_layout(pad=1)

# #%% converting hotstart file to schism output format so that it can be viewd by VisIt
# sh.hotstart_to_outputnc('hotstart2.nc',ini_date,hgrid_fn='./data_in/hgrid2_new.gr3', 
#                          vgrid_fn='./data_in/vgrid2.in.3d',
#                          outname="schout_hotstart2.nc")



