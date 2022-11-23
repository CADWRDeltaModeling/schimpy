# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 10:46:35 2020
Example script to generate a hotstart file for a sediment transport run
@author: zzhang
"""

import schimpy.schism_hotstart as sh
import matplotlib.pyplot as plt

yaml_fn = "hotstart.yaml"
modules = ['TEM','SAL','SED']
hotstart_fn = "hotstart_example2.nc" # output hotstart file
ini_date = '2015-11-18' #'2015-10-01'

# create a hotstart file for SCHISM
h = sh.hotstart(yaml_fn,modules=modules,proj4 ='EPSG:26910')
h.create_hotstart()
hnc = h.nc_dataset
hnc.to_netcdf(hotstart_fn)   

#%% making a 2D surface plot
coll = h.mesh.plot_elems(hnc['tr_el'].values[:,-1,0],clim=(15,22))
cb = plt.colorbar(coll)
plt.axis('off')
plt.axis('equal')
plt.title('Regional Temperature')
plt.tight_layout(pad=1)

#%% converting hotstart file to schism output format so that it can be viewd by VisIt
sh.hotstart_to_outputnc(hotstart_fn,ini_date,hgrid_fn='./hgrid.gr3', 
                         vgrid_fn='./vgrid.in.3d',
                         outname="schout_hotstart_sed.nc")

