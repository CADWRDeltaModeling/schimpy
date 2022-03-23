# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22, 2022
Example script to create SCHISM nudging files and 
visualize the results.
@author: zzhang
"""

from schism_nudging import schism_nudging
from schimpy.schism_mesh import read_mesh
import xarray as xr
import matplotlib.pyplot as plt

yaml_fn = 'testdata\schism_nudging\nudge.yaml'

nudging = schism_nudging.nudging(yaml_fn,proj4 ='EPSG:32610')
nudging.read_yaml()
nudging.create_nudging()

#%% plot temperature
nct = xr.open_dataset('TEM_nu.nc')

fig ,ax = plt.subplots(2,2,sharex=True,figsize=(10,8))
plt.sca(ax[0,0])
coll = nudging.plot(nct.map_to_global_node.values-1,nct.tracer_concentration.isel(
    time=-1,nLevels=-1,one=-1).values,clim=[0,22],ax=ax[0,0],edgecolor='face',
    linewidth=1.5)   #-1 is top 0 is bottom.
cb = plt.colorbar(coll)
plt.axis('equal')
plt.axis('off')
plt.title("nu_temp (top)")

wd = read_mesh('temperature_nudge.gr3')
plt.sca(ax[0,1])
v = wd.nodes[nct.map_to_global_node.values-1,2]
coll = nudging.plot(nct.map_to_global_node.values-1,nct.tracer_concentration.isel(
    time=-1,nLevels=-1,one=-1).values*v,clim=[0,0.0005],ax=ax[0,1],
    edgecolor='face', linewidth=1.5)
cb = plt.colorbar(coll)
plt.axis('equal')
plt.axis('off')
plt.title("nu_temp*weight (top)")

plt.sca(ax[1,0])
coll = nudging.plot(nct.map_to_global_node.values-1,nct.tracer_concentration.isel(
    time=-1,nLevels=0,one=-1).values,clim=[0,22],ax=ax[1,0],edgecolor='face',
    linewidth=1.5)   #-1 is top 0 is bottom.
cb = plt.colorbar(coll)
plt.axis('equal')
plt.axis('off')
plt.title("nu_temp (bottom)")

wd = read_mesh('temperature_nudge.gr3')
plt.sca(ax[1,1])
v = wd.nodes[nct.map_to_global_node.values-1,2]
coll = nudging.plot(nct.map_to_global_node.values-1,nct.tracer_concentration.isel(
    time=-1,nLevels=0,one=-1).values*v,clim=[0,0.0005],ax=ax[1,1],
    edgecolor='face', linewidth=1.5)
cb = plt.colorbar(coll)
plt.axis('equal')
plt.axis('off')
plt.title("nu_temp*weight (bottom)")

plt.tight_layout()

#%% plot salinity

ncs = xr.open_dataset('SAL_nu.nc')

fig ,ax = plt.subplots(2,2,sharex=True,figsize=(10,8))
plt.sca(ax[0,0])
coll = nudging.plot(ncs.map_to_global_node.values-1,ncs.tracer_concentration.isel(
    time=-1,nLevels=-1,one=-1).values,clim=[0,35],ax=ax[0,0],edgecolor='face',
    linewidth=1.5)   #-1 is top 0 is bottom.
cb = plt.colorbar(coll)
plt.axis('equal')
plt.axis('off')
plt.title("nu_salt (top)")

wd = read_mesh('temperature_nudge.gr3')
plt.sca(ax[0,1])
v = wd.nodes[ncs.map_to_global_node.values-1,2]
coll = nudging.plot(ncs.map_to_global_node.values-1,ncs.tracer_concentration.isel(
    time=-1,nLevels=-1,one=-1).values*v,clim=[0,0.0001],ax=ax[0,1],edgecolor='face',
    linewidth=1.5)
cb = plt.colorbar(coll)
plt.axis('equal')
plt.axis('off')
plt.title("nu_salt*weight (top)")

plt.sca(ax[1,0])
coll = nudging.plot(ncs.map_to_global_node.values-1,ncs.tracer_concentration.isel(
    time=-1,nLevels=0,one=-1).values,clim=[0,35],ax=ax[1,0],edgecolor='face',
    linewidth=1.5)   #-1 is top 0 is bottom.
cb = plt.colorbar(coll)
plt.axis('equal')
plt.axis('off')
plt.title("nu_salt (bottom)")

wd = read_mesh('temperature_nudge.gr3')
plt.sca(ax[1,1])
v = wd.nodes[ncs.map_to_global_node.values-1,2]
coll = nudging.plot(ncs.map_to_global_node.values-1,ncs.tracer_concentration.isel(
    time=-1,nLevels=0,one=-1).values*v,clim=[0,0.0001],ax=ax[1,1],edgecolor='face',
    linewidth=1.5)
cb = plt.colorbar(coll)
plt.axis('equal')
plt.axis('off')
plt.title("nu_salt*weight (bottom)")

plt.tight_layout()


