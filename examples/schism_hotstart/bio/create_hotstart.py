# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 11:12:19 2020

@author: zzhang
"""

import schimpy.schism_hotstart as sh
import matplotlib.pyplot as plt
import xarray as xr


yaml_fn = "./hotstart_cosine.yaml"
modules = ['TEM','SAL','COSINE']
hotstart_fn = "hotstart_cosine.nc" # output hotstart file
ini_date = '2015-01-22' 

# create a hotstart file for SCHISM
h = sh.hotstart(yaml_fn,modules=modules,param_nml='data_in/param.nml.clinic',
                crs ='EPSG:26910')
h.create_hotstart()
hnc = h.nc_dataset
hnc.to_netcdf(hotstart_fn)   

#%% making a 2D surface plot for temperature
coll = h.mesh.plot_elems(hnc['tr_el'].values[:,0,0], clim=(0,14)) #clim=[0,35])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional Temperature')
plt.tight_layout(pad=1)

#%% making a 2D surface plot for salinity
coll = h.mesh.plot_elems(hnc['tr_el'].values[:,0,1], clim=(0,35)) #clim=[0,35])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional Salinity')
plt.tight_layout(pad=1)
#%% converting hotstart file to schism output format so that it can be viewd by VisIt
sh.hotstart_to_outputnc(hotstart_fn,ini_date,hgrid_fn='data_in/hgrid.gr3', 
                         vgrid_fn='data_in/vgrid.in.3d',
                         outname="schout_hotstart.nc")





