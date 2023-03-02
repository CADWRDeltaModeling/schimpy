# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 11:12:19 2020

@author: zzhang
"""

import schimpy.schism_hotstart as sh
import matplotlib.pyplot as plt
import xarray as xr
from schimpy.schism_mesh import read_mesh

yaml_fn = "./data_in/hotstart_fromhotstart.yaml"
modules = ['TEM','SAL']
hotstart_fn = "hotstart3.nc" # output hotstart file
ini_date = '2021-04-20' 

#%% first convert mesh2 into utm coordinate system. 
# from schimpy.schism_mesh import read_mesh, write_mesh
# from schimpy.geo_tools import project_fun
# import numpy as np
# mesh1 =  read_mesh('./data_in/hgrid.gr3',
#                    crs='EPSG:26910')   # utm coordinates
# mesh2 =  read_mesh('./data_in/hgrid2.gr3',
#                    crsSG:4326 ')   # lat, lon coordinates
# mesh2_new = sh.project_mesh(mesh2,'EPSG:26910')
# write_mesh(mesh2_new,'./data_in/hgrid2_new.gr3')

#%% create a hotstart file for SCHISM
h = sh.hotstart(yaml_fn,modules=modules,param_nml='./data_in/param.nml.clinic',
                crs ='EPSG:26910')
h.create_hotstart()
hnc = h.nc_dataset
hnc.to_netcdf(hotstart_fn)   

#%% 
hotstart_source = './data_in/hotstart_interpolation/hotstart_it=2688000.nc'
source_data = xr.open_dataset(hotstart_source)
mesh1 = read_mesh('./data_in/hotstart_interpolation/hgrid.gr3',
                  './data_in/hotstart_interpolation/vgrid.in.3d')

mesh2 = read_mesh('./data_in/hgrid.gr3','./data_in/vgrid.in.3d')

#%%
# coll = h.mesh.plot_nodes(h.nc_dataset['elevation'].values,clim=[0,2])
# cb = plt.colorbar(coll)
# #plt.axis('off')
# plt.axis('equal')
# plt.title('Regional Temperature')
# plt.tight_layout(pad=1)

#%% making a 2D surface plot
fig,ax = plt.subplots(1,2,sharex=True)
plt.sca(ax[0])
coll = h.mesh.plot_elems(hnc['tr_el'].values[:,0,-1], clim=[-5,20],ax=ax[0])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional Temperature')
plt.tight_layout(pad=1)

# making a 2D surface plot
plt.sca(ax[1])
coll = mesh1.plot_elems(source_data['tr_el'].values[:,0,-1], clim=[-5,20],ax=ax[1]) #, clim=(14,18)) #clim=[0,35])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional Temperature')
plt.tight_layout(pad=1)

#%% check tke
fig,ax = plt.subplots(1,2,sharex=True)
plt.sca(ax[0])
coll = h.mesh.plot_nodes(hnc['q2'].values[:,0],ax=ax[0],clim=(0,6e-4))
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional q2')
plt.tight_layout(pad=1)

# making a 2D surface plot
plt.sca(ax[1])
coll = mesh1.plot_nodes(source_data['q2'].values[:,0],ax=ax[1],clim=(0,6e-4)) #, clim=(14,18)) #clim=[0,35])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional q2')
plt.tight_layout(pad=1)

#%%
fig,ax = plt.subplots(1,2,sharex=True)
plt.sca(ax[0])
coll = h.mesh.plot_nodes(hnc['xl'].values[:,0],ax=ax[0],clim=(0,6e-2))
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional q2')
plt.tight_layout(pad=1)

# making a 2D surface plot
plt.sca(ax[1])
coll = mesh1.plot_nodes(source_data['xl'].values[:,0],ax=ax[1],clim=(0,6e-2)) #, clim=(14,18)) #clim=[0,35])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional q2')
plt.tight_layout(pad=1)

#%% check u, v and w
fig,ax = plt.subplots(1,2,sharex=True)
plt.sca(ax[0])
coll = h.mesh.plot_edges(hnc['su2'].values[:,-1],ax=ax[0],size=2)
coll.set_clim((-0.01,0.01))
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('u')
plt.tight_layout(pad=1)

# making a 2D surface plot
plt.sca(ax[1])
coll = mesh1.plot_edges(source_data['su2'].values[:,-1], ax=ax[1],size=2) #, clim=(14,18)) #clim=[0,35])
coll.set_clim((-0.01,0.01))
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('u')
plt.tight_layout(pad=1)

#%%
fig,ax = plt.subplots(1,2,sharex=True)
plt.sca(ax[0])
coll = h.mesh.plot_elems(hnc['we'].values[:,1],ax=ax[0],clim=[-0.005,0.005])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('w')
plt.tight_layout(pad=1)

# making a 2D surface plot
plt.sca(ax[1])
coll = mesh1.plot_elems(source_data['we'].values[:,1],ax=ax[1],clim=[-0.005,0.005]) #, clim=(14,18)) #clim=[0,35])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('w')
plt.tight_layout(pad=1)

#%% check water surface
fig,ax = plt.subplots(1,2,sharex=True)
plt.sca(ax[0])
coll = h.mesh.plot_nodes(hnc['eta2'].values,ax=ax[0],clim=[0,2])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional Temperature')
plt.tight_layout(pad=1)

# making a 2D surface plot
plt.sca(ax[1])
coll = mesh1.plot_nodes(source_data['eta2'].values,ax=ax[1],clim=[0,2]) #, clim=(14,18)) #clim=[0,35])
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('water surface')
plt.tight_layout(pad=1)

#%% converting hotstart file to schism output format so that it can be viewd by VisIt
sh.hotstart_to_outputnc('hotstart3.nc',ini_date,hgrid_fn='./data_in/hgrid.gr3', 
                         vgrid_fn='./data_in/vgrid.in.3d',
                         outname="schout_hotstart3.nc")


#%% compare hotstar2.nc and hotstart3.nc
nc2 = xr.open_dataset('hotstart2.nc')
nc3 = xr.open_dataset('hotstart3.nc')

#%%
fig ,ax = plt.subplots(1,2,sharex=True)

coll1 = h.mesh.plot_elems(nc2.tr_el.isel(ntracers=1).values[:,-1],clim=[0,1],ax=ax[0])
coll2 = h.mesh.plot_elems(nc3.tr_el.isel(ntracers=1).values[:,-1],clim=[0,1],ax=ax[1])
cb = plt.colorbar(coll2)



