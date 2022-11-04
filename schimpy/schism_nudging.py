# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 10:19:56 2021

@author: zzhang
To be implemented: other weights function (length scale) from OI in additional 
to gaussian
"""

import copy
import datetime
import yaml
import pandas as pd
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from schimpy.schism_mesh import read_mesh, write_mesh
from schimpy.geo_tools import utm2ll
from shapely.geometry import Polygon
from schimpy  import Interp2D
import time as timer
from vtools.data.vtime import hours,days

class nudging(object):
    """
    A class to create schism nudging
    """
    
    def __init__(self, yaml_fn,crs=None,suffix=None, **kwargs):
        self.yaml_fn = yaml_fn
        self._interpolant = ['nearest','inverse distance']
        self._kernel = ['gaussian']   
        self.crs = crs
        self.output_suffix = suffix
           
    def read_yaml(self):
        """
        read yaml and load mesh grid
        """
        with open(self.yaml_fn) as file:
            info = yaml.load(file, Loader=yaml.FullLoader)
        nudging_info = info['nudging']
        self.info = nudging_info
        self.start_date = nudging_info['start_date']
        self.rnday = nudging_info['rnday']
        self.end_date = self.start_date + datetime.timedelta(self.rnday)
        self.nudge_step = nudging_info['step_nu_tr']
        self.datetime = pd.date_range(self.start_date, self.end_date,
                                      freq=self.nudge_step)
        self.time = pd.to_datetime(self.datetime.values)- \
            pd.to_datetime(self.start_date)
        self.time_seconds = self.time.total_seconds().astype(int)
        self.hgrid_fn = nudging_info['hgrid_input_file']
        self.vgrid_fn = nudging_info['vgrid_input_file']
        self.default_value = nudging_info['default']
        self.vgrid_version = nudging_info['vgrid_version']
        self.mesh = read_mesh(self.hgrid_fn,self.vgrid_fn, self.vgrid_version)
        self.node_x = self.mesh.nodes[:,0]
        self.node_y = self.mesh.nodes[:,1]
        self.node_z = self.mesh.nodes[:,2] 
        self.nnode = self.mesh.n_nodes()
        self.nvrt = self.mesh.n_vert_levels
        self._mesh_gpd = None
        self._z = None
        if self.crs is None:
            if 'crs' in nudging_info['crs']:
                self.crs = nudging_info['crs']
            else:
                self.crs = 'EPSG:26910' # this is required because ROMS only provides lat, lon. 
                print("crs not specified, and assigned to the default crs for UTM10N")
                    
        if self.output_suffix is None:
            if 'output_suffix' in nudging_info:
                self.output_suffix = nudging_info['output_suffix']
            else:
                self.output_suffix = None
        
    
    @property
    def mesh_gpd(self):
        if self._mesh_gpd is None:
            self._mesh_gpd = self.mesh.to_geopandas('point')
        return self._mesh_gpd
    
    @property
    def z(self):
        if self._z is None:
            self._z = self.mesh.build_z()
        return self._z

    def create_nudging(self, create_file=True):
        """
        Parameters
        ----------
        suffix : TYPE, optional
            Options to give a suffix to the output nuding files. 
            The default is 'nu.nc'
        create_file : TYPE, optional
            Options to create nudging files. The default is True.
        Returns
        -------
        None.

        """
        weights_comb, values_comb, imap_comb = self.get_nudging_comb()
        
        weights_var, values_var, imap_var, var_v = self.organize_nudging(
            weights_comb, values_comb, imap_comb) 
        
        self.concatenate_nudge(weights_var, 
                values_var,
                imap_var,
                var_v,
                create_file=True)    
    
    def organize_nudging(self, weights_comb, values_comb, imap_comb):   
        weights_list = [] #dimension [region, var, node]
        values_list = [] # dimensino [region, var, time, map_node, nlevel]
        imap_list = [] # dimension [region, var, map_node] 
        # get a list of variables
        var_list = np.array([])
        for p in self.info['polygons']:
            var_list = np.append(var_list, 
                                 [v['name'] for v in p['interpolant']['variables']])
        var_list = np.unique(var_list)
                
        # eliminate variables that don't exist in records
        weights_list = [] #dimension [region, var, node]
        values_list = [] # dimensino [region, var, time, map_node, nlevel]
        imap_list = [] # dimension [region, var, map_node]
        for pi, p in enumerate(self.info['polygons']):            
            weights_v = []
            values_v = []
            imap_v = []
            var_v = []
            var_info = p['interpolant']['variables']
            var_name = np.array([v['name'] for v in var_info])
             
            for vi, v in enumerate(var_list):
                # reorder according to ['temperature','salinity']
                if v in var_name:
                    pos = np.where(var_name==v)[0][0]
                    if len(imap_comb[pi][pos])>0:
                        weights_v.append(weights_comb[pi][pos])
                        values_v.append(values_comb[pi][pos])
                        imap_v.append(imap_comb[pi][pos])  
                        var_v.append(v)
                    else:
                        weights_v.append([])
                        values_v.append([])
                        imap_v.append([])
                        var_v.append([])
                else:
                    weights_v.append([])
                    values_v.append([])
                    imap_v.append([])
                    var_v.append([])                    
            #if np.any(weights_v):
            weights_list.append(weights_v)
            values_list.append(values_v)
            imap_list.append(imap_v)        
        
        values_var = []  #rearrange to dim [var, region, time, map_node, nlevel]
        imap_var = [] # rearrange to dim [var, region, map_node]
        weights_var = []        
        
        #weights_var = np.array(weights_list).transpose([1,0,2])#rearrange to dim [var, region,node]
        for i, v in enumerate(var_list):
                values_var.append([vl[i] for vl in values_list])
                imap_var.append([im[i] for im in imap_list])
                weights_var.append([wv[i] for wv in weights_list])
        var_v = []
        vi = []
        for i, v in enumerate(var_list):
            if np.any(weights_var[i]):
                vi.append(i)
                var_v.append(v)

        weights_var = [np.array(weights_var[i]) for i in vi]
        values_var = [values_var[i] for i in vi]
        imap_var = [imap_var[i] for i in vi]                           
        return weights_var, values_var, imap_var, var_v

    def get_nudging_comb(self):
        weights_comb = []
        values_comb = []
        imap_comb = []
        for p in self.info['polygons']:
            print("creating nudging for polygon %s"%p['name'])
            weights, values, imap = self.create_region_nudging(p)
            weights_comb.append(weights)
            values_comb.append(values)
            imap_comb.append(imap)
        return weights_comb, values_comb, imap_comb    
    
    def concatenate_nudge(self, weights_var, values_var, imap_var,
                          var_newlist, create_file=True):
            
        for i, v in enumerate(var_newlist):
            # merge all the mapping values
            imap_merged = []
            for l in imap_var[i]:
                imap_merged = np.append(imap_merged, l)
            imap_merged = np.unique(imap_merged).astype(int)
            
            # finding the corresponding indices for each region (map from imap_var to imap_merged)
            idx_list = np.zeros((len(imap_merged),len(self.info['polygons']))
                                )*np.nan
            for j, im in enumerate(imap_merged):
                for k, l in enumerate(imap_var[i]):
                    idx = np.where(l==im)[0]
                    if len(idx)>0:
                        idx_list[j,k] = int(idx[0])
                    else:
                        idx_list[j,k] = np.nan
            
            # check to see if there is any overlap among the regions, if yes, 
            # take the weighted average  
            weights_merged = np.zeros_like(weights_var[0][0])
            values_merged = np.zeros((len(self.time),len(imap_merged),
                                          self.nvrt)) 
            for im in range(len(imap_merged)):
                idx_r = np.where(~np.isnan(idx_list[im,:]))[0].astype(int)
                if len(idx_r)>1: # overlapping regions
                    values_sum = []
                    weights_sum = []  
                    weights_sum2 = []
                    for r in idx_r:
                        w = weights_var[i][r,imap_merged[im]]
                        data_values = values_var[i][r][:,int(idx_list[im,r]),:]
                        if np.isnan(data_values).any():
                            print("warning: nan values are detected in the overlapping region. The nudging values are set to -9999.0 for these regions. ")
                        values_sum.append(w*data_values)
                        weights_sum.append(w)
                        weights_sum2.append(w**2)
                    values_temp = np.sum(values_sum,axis=0)/sum(weights_sum)
                    if np.shape(values_temp)[0]>len(self.time):
                        values_merged[:,im,:] = values_temp[:len(self.time),:]
                    else:
                        values_merged[:,im,:] = values_temp
                    #values_merged[:,im,:] = np.sum(values_sum,axis=0)/sum(weights_sum)
                    weights_merged[imap_merged[im]] = sum(weights_sum2)/sum(weights_sum)
                else:
                    r = idx_r[0]
                    weights_merged[imap_merged[im]]  = weights_var[i][r,imap_merged[im]]
                    values_merged[:,im,:] = values_var[i][r][:,int(idx_list[im,r]),:]
            
            # for the last time, make sure that all nan values are set to -9999.0
            values_merged[np.isnan(values_merged)] = -9999.0
                    
            
            if self.output_suffix is None:
                suffix = ''
            else:
                suffix = self.output_suffix 
            
            if v== 'temperature':
                nudging_fn = "TEM_nu_%s.nc"%suffix     
            elif v=='salinity':
                nudging_fn = "SAL_nu_%s.nc"%suffix
            else:
                raise NotImplementedError(
                    "write %s as nuding nc not implemented"%v)
                           
            rootgrp_T = Dataset(nudging_fn, "w", format="NETCDF4")
            node_dim = rootgrp_T.createDimension("node", len(imap_merged))
            nv_dim = rootgrp_T.createDimension("nLevels", self.nvrt)
            one_dim = rootgrp_T.createDimension("one", 1)
            itime_dim = rootgrp_T.createDimension("time", len(self.time_seconds))
            
            itime_id_T = rootgrp_T.createVariable("time","f",("time",))
            id_map_T = rootgrp_T.createVariable("map_to_global_node","i",
                                                ("node",))
            ivar_T = rootgrp_T.createVariable("tracer_concentration","f",
                                            ("time","node","nLevels","one",))  
                
            id_map_T[:] = np.array(imap_merged) +1  #in schism indices are 1 based   
            itime_id_T[:] = self.time_seconds            
            
           # current var dimension [time, map_node, nlevel]
            ivar_T[:,:,:,:] = values_merged[:,:,:,np.newaxis] 
            rootgrp_T.close()
            print("%s file created"%nudging_fn)  
            
            if self.output_suffix is not None:
                write_mesh(self.mesh, 
                           fpath_mesh="%s_nudge.gr3"%v, 
                           node_attr=weights_merged)
            else:
                write_mesh(self.mesh, 
                           fpath_mesh="%s_nudge_%s.gr3"%(v,self.output_suffix), 
                       node_attr=weights_merged)
        #return values_merged, weights_merged, imap_merged
    
    def create_region_nudging(self, region_info):
        if region_info['type'] == 'roms':
            weights, values_list, imap, time = self.gen_nudge_roms(region_info)
            #values_list =  [vl[:,imap,:] for vl in values_list]
            nvar = len(set([v['name'] for v in 
                                region_info['interpolant']['variables']]))
            weights_list = np.broadcast_to(weights,[nvar,len(weights)])
            imap_list = np.broadcast_to(imap,[nvar,len(imap)])
        elif 'obs' in region_info['type']:  # this is 2D interpolation
            weights_list, values, imap_list = self.gen_nudge_obs(region_info)
            values_list = [np.broadcast_to(v,(self.nvrt, np.shape(v)[0],
                                              np.shape(v)[1])) 
                           if np.any(v) else [] for v in values]
            values_list = [np.transpose(v,[1,2,0]) 
                           if np.any(v) else [] for v in values_list]
        else:
            raise("region type not implemented:%s"%region_info['type'])
        return weights_list, values_list, imap_list
    
    def gen_nudge_roms(self,region_info):
        rjunk = 9998 #!Define junk value for sid; the test is abs()>rjunk
        hr_char = ['03','09','15','21'] #each day has 4 starting hours in ROMS
        small1 = 1.e-2 #used to check area ratios
        weights, obs_df = self.gen_region_weight(region_info['attribute'],
                                                 region_info['vertices'])
        istart = pd.to_datetime(self.start_date)
        istart_year = istart.year
        istart_mon = istart.month
        istart_day = istart.day
        nndays = (self.end_date - self.start_date).days    
        dtout = region_info['interpolant']['dt'] #time step in .nc 
        nt_out = int(pd.Timedelta(self.nudge_step)/pd.Timedelta(dtout)) #output stride
        variable_info = region_info['interpolant']['variables']
        for vi in variable_info:
            if vi['name'] == 'temperature':
                tem_outside = float(vi['none_values'])
                if 'offset' in vi.keys():
                    tem_offset = float(vi['offset'])
            if vi['name'] == 'salinity':
                sal_outside = float(vi['none_values'])        
        ncfile1 = region_info['interpolant']['data']   
        ncfile1 = ncfile1.replace("'","")  # remove the quotation marks if they exist
        irecout = 0
        irecout2 = 0
        #Define limits for variables for sanity checks
        tempmin=0 
        tempmax=30
        saltmin=0
        saltmax=37
        vmag_max=10 #max. |u| or |v|
        ssh_max=8   #max. of |SSH|
        imap = np.nan*np.zeros(self.nnode)
        nodes = np.where(weights>0)[0]      
       
        # # instead of reading lat, lon from hgrid.ll, I can also convert 
        # # hgrid.gr3 to lat, lon.
        # with open('hgrid.ll','r') as reader:
        #     lines = reader.readlines()            
        # xl=[]
        # yl=[]
        # for i in range(self.nnode):
        #     llist = lines[2+i].split()
        #     xl.append(float(llist[1]))
        #     yl.append(float(llist[2]))    
            
        utm_xy = np.array([self.node_x,self.node_y])
        lonlat= utm2ll(utm_xy,self.crs) #roms grid has lat, lon coordinates. 
        xl = lonlat[0]
        yl = lonlat[1]

        # the 'include2' variable is replaced by 'weights'
        # include indicates if spatial interpolation should be applied 
        # include2=0 means skipping interpolation
        # with open('include.gr3','r') as reader:
        #     lines = reader.readlines()            
        # include2 = []
        # for i in range(self.nnode):
        #     llist = lines[2+i].split()
        #     include2.append(float(llist[-1]))
        
        temperature = []
        salinity = []
        output_day = []
        var_name = {'temperature':'temp',
                    'salinity': 'salt',
                    'velocity_u':'uvel',
                    'velocity_v':'vvel',
                    'elevation':'ssh'}
        #var_list = [var_name[v] for v in region_info['interpolant']['variable']]
        
        #start = timer.time()
        for d in range(nndays+1): 
            date = datetime.date(istart_year,istart_mon,istart_day) + \
                datetime.timedelta(d)
            day = (datetime.date(istart_year,istart_mon,istart_day) - 
                self.start_date).days + d
            print('Time out (days)=%d'%day)  
            for hr in hr_char:
                ncfile = "%s%4d%02d%02d%s.nc"%(ncfile1,date.year,date.month,
                                            date.day,hr)                
                ncdata = xr.open_dataset(ncfile)
                #print("time1=%f"%(timer.time() - start))
                
                if d ==0 and hr == hr_char[0]:                    
                    # lon, lat = np.meshgrid(ncdata['lon'].values,
                    #                        ncdata['lat'].values)
                    lat, lon = np.meshgrid(ncdata['lat'].values,
                                           ncdata['lon'].values)
                    hnc = ncdata['depth']
                    lon=lon-360e0 # convert to our long
                    #Compute z-coord. (assuming eta=0)
                    #WARNING: In zm(), 1 is bottom; ilen is surface (SELFE convention)
                    zm = -hnc[-1::-1].values
                    ixlen = np.shape(lat)[0]
                    iylen = np.shape(lat)[1]
                    ilen = len(zm)
                
                if d==0 and hr==hr_char[0]:
                    ilo = 5 #5th record is 00:00 PST UCT-8
                else:
                    ilo = 1 #there is an overlap between files

                time = ncdata['time'].values
                
                if time[-1]< np.datetime64(self.start_date)+ \
                    np.timedelta64(8,'h'):
                    continue
                elif time[0]> np.datetime64(self.end_date) + \
                    np.timedelta64(8,'h'):
                    break                
                
                if time[0] != np.datetime64(date) + np.timedelta64(int(hr),'h'):
                    print("Initial time in %s does not match the filename: %s: modification will be applied "%(ncfile,time[0]))
                    dt_adj = np.datetime64(date) + np.timedelta64(int(hr),'h') - time[0]
                    ncdata['time'] = ncdata['time'] + dt_adj 
                    time = ncdata.time.values

                for t in time[ilo:]:
                    ctime = pd.to_datetime(t)
                    dt = ctime-pd.to_datetime(self.start_date) 
                    dt_in_days = dt.total_seconds()/86400.
                    print(f"Time out at hr stage (days)={ctime}, dt ={dt_in_days} ")                  
                    irecout += 1 

                    #Make sure it2=ilo is output (output at precisely the time interval)
                    # Read T,S,u,v,SSH
                    # WARNING! Make sure the order of vertical indices is 1 
                    # at bottom, ilen at surface; revert if necessary!
                    if (irecout-1)%nt_out !=0:
                        continue

                    irecout2 += 1
                    print("irecount: %d"%irecout)
                           
                    #if  'temp' in var_list:
                    temp = ncdata['temp'].sel(time=t).transpose(
                        'lon','lat','depth').values[:,:,-1::-1]
                    #if 'salt' in var_list:
                    salt = ncdata['salt'].sel(time=t).transpose(
                        'lon','lat','depth').values[:,:,-1::-1]      
                    # if 'uvel' in var_list:
                    #     uvel = ncdata['uvel'].sel(time=t).transpose(
                    #         'lon','lat','depth').values[:,:,-1::-1] 
                    # if 'vvel' in var_list:
                    #     vvel = ncdata['vvel'].sel(time=t).transpose(
                    #         'lon','lat','depth').values[:,:,-1::-1]   
                    # if 'ssh' in var_list:
                    #     ssh = ncdata['ssh'].sel(time=t).transpose(
                    #         'lon','lat')                                     
            
                    kbp=np.ones((ixlen,iylen))  #wet or dry flag
                    
                    dry_xy = np.where( (abs(salt[:,:,-1])>rjunk) | 
                                    (np.isnan(salt[:,:,-1])) )
                    kbp[dry_xy[0],dry_xy[1]] = -1
                    wet_xy = np.where(kbp==1)
                    
                    
                    klev0 = [np.where(~np.isnan(salt[ix,iy,:]))[0][0] for 
                             (ix,iy) in zip(wet_xy[0],wet_xy[1])]
                    
                    #print("time2=%f"%(timer.time() - start))
                    for wi in range(len(klev0)): 
                        salt[wet_xy[0][wi],wet_xy[1][wi],:klev0[wi]] = \
                            salt[wet_xy[0][wi],wet_xy[1][wi],klev0[wi]]
                        temp[wet_xy[0][wi],wet_xy[1][wi],:klev0[wi]] = \
                            temp[wet_xy[0][wi],wet_xy[1][wi],klev0[wi]]   
                    
                    
                    # for i in range(ixlen):
                    #     for j in range(iylen):
                    #         if (abs(salt[i,j,-1]>rjunk)) or np.isnan(salt[i,j,-1]):
                    #             kbp[i,j]=-1 #dry
                    #         if kbp[i,j]==1:
                    #             #extend to the bottom
                    #             klev0 = np.where(~np.isnan(salt[i,j,:]))[0][0]
                    #             salt[i,j,:klev0] = salt[i,j,klev0]
                    #             temp[i,j,:klev0] = temp[i,j,klev0]
                                # if 'uvel' in var_list:
                                #     uvel[i,j,:klev0] = uvel[i,j,klev0]
                                # if 'vvel' in var_list:
                                #     uvel[i,j,:klev0] = uvel[i,j,klev0]                                    
                            
                            # double check to make sure there are no out of bound values
                            # if any(salt[i,j,:]<saltmin) or \
                            #     any(salt[i,j,:]>saltmax) or \
                            #     any(temp[i,j,:]<tempmin) or \
                            #     any(temp[i,j,:]>tempmax):
                            #     # any(abs(uvel[i,j,:]))>vmag_max or \
                            #     # any(abs(vvel[i,j,:]))>vmag_max or \
                            #     # abs(ssh[i,j])>ssh_max
                            #     print("Fatal: no valid values:%s,%d,%d for salt or temp at"
                            #           (ncfile,i,j))                                                      

                    # Compute S,T etc@ invalid pts based on nearest neighbor
                    # Search around neighborhood of a pt                    
                    dryij = np.where(kbp==-1)
                    dryi = dryij[0]
                    dryj = dryij[1]
                    wetij = np.where(kbp==1)
                    weti = wetij[0]
                    wetj = wetij[1] 
                    
                    if wetij[0].size==0:
                        print("no ROMS data available for: %s"%str(ctime))
                        salt = -9999.0*np.ones_like(salt)
                        temp = -9999.0*np.ones_like(temp)
                    else: 
                        for i, j in zip(dryi,dryj):
                            distij = np.abs(weti-i) + np.abs(wetj-j)                        
                            #m = np.where(distij==min(distij))[0][0]
                            m = np.argmin(distij)    
                        # salt[dryi,dryj,:] = salt[weti[m],wetj[m],:]
                        # temp[dryi,dryj,:] = temp[weti[m],wetj[m],:]
                            salt[i,j,:]=salt[weti[m],wetj[m],:]
                            temp[i,j,:]=temp[weti[m],wetj[m],:]    
                            #uvel(i,j,:ilen)=uvel(weti[m],wetj[m],:ilen)
                            #vvel(i,j,:ilen)=vvel(weti[m],wetj[m],:ilen)
                            #ssh(i,j)=ssh(weti[m],wetj[m])
                        
                        # for i in range(ixlen):
                        #     for j in range(iylen):
                        #         if kbp[i,j]==-1:   # invalid pts (dry) 
                        #             #Compute max possible tier # for neighborhood
                        #             mmax=max(i,ixlen-i-1,j,iylen-j-1)
                                    
                        #             m = 0 
                        #             i3 = np.nan
                        #             j3 = np.nan
                        #             #while True: #starting from (i,j) and search on both sides  
                        #             for m in range(0,mmax+1):                                  
                        #                 for ii in range(max(-m,-i),
                        #                                 min(m,ixlen-i-1)):
                        #                     i3 = max(0,min(ixlen-1,i+ii))
                        #                     for jj in range(max(-m,-j),
                        #                                     min(m,iylen-j-1)):
                        #                         j3=max(0,min(iylen-1,j+jj))
                        #                         if kbp[i3,j3]==1:
                        #                             i1=i3
                        #                             j1=j3
                        #                             break 
                        #                     if kbp[i3,j3]==1:
                        #                          break
                        #                 if ~np.isnan(i3) and ~np.isnan(j3): 
                        #                     if kbp[i3,j3]==1:
                        #                          break
                        #                 if m == mmax:
                        #                     print("Max. exhausted:%d,%d,%d"%(
                        #                         i,j,mmax))
                        #                     print('kbp')
                        #                     for ii in range(ixlen):
                        #                         for jj in range(iylen):
                        #                             print(ii,jj,kbp(ii,jj))
                        #                     raise("Max. exhausted!")
                        #                 #m += 1 
                                      
                        #             salt[i,j,:ilen]=salt[i1,j1,:ilen]
                        #             temp[i,j,:ilen]=temp[i1,j1,:ilen]
                                    # double check to make sure there are no out of bound values
                                    # if any(salt[i,j,:]<saltmin) or \
                                    #     any(salt[i,j,:]>saltmax) or \
                                    #     any(temp[i,j,:]<tempmin) or \
                                    #     any(temp[i,j,:]>tempmax):
                                    #     print("Fatal: no valid values:%s,%d,%d for salt or temp at"
                                    #           (ncfile,i,j))  
                        print("dealing with horizontal nan values!") 
                        #print("time3=%f"%(timer.time() - start))
                    
                    if d ==0 and hr == hr_char[0] and t==time[ilo]:
                        ixy=np.zeros((self.nnode,3))*np.nan
                        wild = np.zeros(2)
                        wild2 = np.zeros((4,2))
                        arco = np.zeros((4,self.nnode))

                        x1 = lon[0:-1,0:-1]
                        x2 = lon[1:,0:-1]
                        x3 = lon[1:,1:]
                        x4 = lon[0:-1,1:]
                        y1 = lat[0:-1,0:-1]
                        y2 = lat[1:,0:-1]
                        y3 = lat[1:,1:]
                        y4 = lat[0:-1,1:]   
                        b1=abs(self.signa(x1,x2,x3,y1,y2,y3))
                        b2=abs(self.signa(x1,x3,x4,y1,y3,y4))
                       
                        for i in range(self.nnode): 
                        #for i in nodes:
                            if weights[i]!=0.0:  # the following code applies a spatial interpolation    
                                         
                                a1=abs(self.signa(xl[i],x1,x2,yl[i],y1,y2))
                                a2=abs(self.signa(xl[i],x2,x3,yl[i],y2,y3))
                                a3=abs(self.signa(xl[i],x3,x4,yl[i],y3,y4))
                                a4=abs(self.signa(xl[i],x4,x1,yl[i],y4,y1))
                                rat=abs(a1+a2+a3+a4-b1-b2)/(b1+b2)
                                xy = np.where(rat<small1)
                                ixy[i,0] = xy[0][0]
                                ixy[i,1] = xy[1][0]
                                x1i = x1[int(ixy[i,0]),int(ixy[i,1])]
                                x2i = x2[int(ixy[i,0]),int(ixy[i,1])]
                                x3i = x3[int(ixy[i,0]),int(ixy[i,1])]
                                x4i = x4[int(ixy[i,0]),int(ixy[i,1])]
                                y1i = y1[int(ixy[i,0]),int(ixy[i,1])]
                                y2i = y2[int(ixy[i,0]),int(ixy[i,1])]
                                y3i = y3[int(ixy[i,0]),int(ixy[i,1])]
                                y4i = y4[int(ixy[i,0]),int(ixy[i,1])]
                                a1i = a1[int(ixy[i,0]),int(ixy[i,1])]
                                a2i = a2[int(ixy[i,0]),int(ixy[i,1])]
                                a3i = a3[int(ixy[i,0]),int(ixy[i,1])]
                                a4i = a4[int(ixy[i,0]),int(ixy[i,1])]                                
                                
                                # Find a triangle
                                intri=0  #flag: inside the triangle
                                for l in range(2): #split quad
                                  ap=abs(self.signa(xl[i],x1i,x3i,yl[i],y1i,y3i))
                                  if l==0:  #nodes 1,2,3
                                      bb=abs(self.signa(x1i,x2i,x3i,y1i,y2i,y3i))
                                      wild[l]=abs(a1i+a2i+ap-bb)/bb
                                      if wild[l]<small1*5:
                                          intri=1
                                          arco[0,i]=max(0.,min(1.,a2i/bb))
                                          arco[1,i]=max(0.,min(1.,ap/bb))
                                          break
                                  else: #nodes 1,3,4
                                      bb=abs(self.signa(x1i,x3i,x4i,y1i,y3i,y4i))
                                      wild[l]=abs(a3i+a4i+ap-bb)/bb
                                      if wild[l]<small1*5:
                                          intri=2
                                          arco[0,i]=max(0.,min(1.,a3i/bb))
                                          arco[1,i]=max(0.,min(1.,a4i/bb))
                                          break
      
                                if(intri==0):
                                    raise('Cannot find a triangle:', wild)
          
                                ixy[i,2]=intri
                                arco[2,i]=max(0.,min(1.,1-arco[0,i]-arco[1,i]))     
                            
                                # for ix in range(ixlen-1):   
                                #     for iy in range(iylen-1):  
                                #         x1=lon[ix,iy]
                                #         x2=lon[ix+1,iy]
                                #         x3=lon[ix+1,iy+1]
                                #         x4=lon[ix,iy+1]
                                #         y1=lat[ix,iy]
                                #         y2=lat[ix+1,iy]
                                #         y3=lat[ix+1,iy+1]
                                #         y4=lat[ix,iy+1]
                                #         a1=abs(self.signa(xl[i],x1,x2,yl[i],y1,y2))
                                #         a2=abs(self.signa(xl[i],x2,x3,yl[i],y2,y3))
                                #         a3=abs(self.signa(xl[i],x3,x4,yl[i],y3,y4))
                                #         a4=abs(self.signa(xl[i],x4,x1,yl[i],y4,y1))
                                #         b1=abs(self.signa(x1,x2,x3,y1,y2,y3))
                                #         b2=abs(self.signa(x1,x3,x4,y1,y3,y4))
                                #         rat=abs(a1+a2+a3+a4-b1-b2)/(b1+b2)
                                #         if rat<small1:
                                #           ixy[i,0]=ix
                                #           ixy[i,1]=iy
                                #           # Find a triangle
                                #           intri=0  #flag: inside the triangle
                                #           for l in range(2): #split quad
                                #             ap=abs(self.signa(xl[i],x1,x3,yl[i],y1,y3))
                                #             if l==0:  #nodes 1,2,3
                                #                 bb=abs(self.signa(x1,x2,x3,y1,y2,y3))
                                #                 wild[l]=abs(a1+a2+ap-bb)/bb
                                #                 if wild[l]<small1*5:
                                #                     intri=1
                                #                     arco[0,i]=max(0.,min(1.,a2/bb))
                                #                     arco[1,i]=max(0.,min(1.,ap/bb))
                                #                     break
                                #             else: #nodes 1,3,4
                                #                 bb=abs(self.signa(x1,x3,x4,y1,y3,y4))
                                #                 wild[l]=abs(a3+a4+ap-bb)/bb
                                #                 if wild[l]<small1*5:
                                #                     intri=2
                                #                     arco[0,i]=max(0.,min(1.,a3/bb))
                                #                     arco[1,i]=max(0.,min(1.,a4/bb))
                                #                     break
                
                                #           if(intri==0):
                                #               raise('Cannot find a triangle:', wild)
                    
                                #           ixy[i,2]=intri
                                #           arco[2,i]=max(0.,min(1.,1-arco[0,i]-arco[1,i]))
                        # the following step was calculated at every time step
                        # in the original code but can significanltyly slow down 
                        # the python script.
                        # so it is only calculated once with the exception that
                        # the dry cells are filled at every time step
                        i_out = np.where(np.isnan(ixy[:,0]))[0]
                        i_in = np.where(~np.isnan(ixy[:,0]))[0]
                        npout = len(i_in)
                        imap = i_in
                        ix = ixy[i_in,0].astype(int)
                        iy = ixy[i_in,1].astype(int)
                        intri = ixy[i_in,2].astype(int)
                    
                        lev = np.ones((npout,self.nvrt))*-99
                        vrat = np.zeros((npout,self.nvrt))*-99
                        for il, ii in enumerate(i_in):
                            for k in range(self.nvrt):
                                if(kbp[ix[il],iy[il]]==-1): #ROMS dry cell
                                    lev[il,k] = int(ilen-1-1)  #
                                    vrat[il,k] = 1
                                elif self.z[ii,k]<=zm[int(kbp[ix[il],iy[il]])]:
                                    lev[il,k] = int(kbp[ix[il],iy[il]]-1)
                                    vrat[il,k] = 0
                                elif self.z[ii,k]>=zm[ilen-1]:
                                    lev[il,k] = ilen-1-1
                                    vrat[il,k] = 1
                                else:
                                    lev[il,k] = -99 #flag
                                    for kk in range(ilen-1):
                                        if (self.z[ii,k]>=zm[kk] and 
                                            self.z[ii,k]<=zm[kk+1]):
                                            lev[il,k] = kk
                                            vrat[il,k] = (zm[kk]-self.z[ii,k])/(zm[kk]-zm[kk+1])
                                            break                              
                                if lev[il,k] == -99:
                                    raise('Cannot find a level:', ii,k, 
                                          self.z[ii,k],
                                          zm)
                        kbp = kbp.astype(int)
                        lev = lev.T.astype(int)
                        vrat = vrat.T   
                        
                        wild2 = np.zeros((self.nvrt, npout, 4,2))   #initialize interpolation coefficients 
                        ix = np.broadcast_to(ix,(self.nvrt, npout)) 
                        iy = np.broadcast_to(iy,(self.nvrt, npout))
                        intri1 = np.where(intri==1)[0]
                        intri2 = np.where(intri==2)[0]
                        i_in_1 = i_in[intri1]
                        i_in_2 = i_in[intri2]
                        #print("time=%f"%(timer.time() - start))
                    
                        #tempout=-999*np.ones((self.nvrt, self.nnode))
                        #saltout=-999*np.ones((self.nvrt, self.nnode)) 
                        tempout=np.empty((self.nvrt, self.nnode))
                        saltout=np.empty((self.nvrt, self.nnode)) 
                        tempout[:,i_out] = tem_outside
                        saltout[:,i_out] = sal_outside   
                    
                    id_dry = np.where(kbp[ix,iy]==-1)[0] #ROMS dry cell
                    lev[:,id_dry] = int(ilen-1-1)
                    vrat[:,id_dry] = 1
                    # id_deeper = np.where(self.z[i_in,:]<zm[kbp[ix,iy]]) #deeper than ROMS
                    # lev[id_deeper] = int(kbp[ix,iy]-1)
                    # vrat[id_deeper] = 0
                    # id_shallower = np.where(self.z[i_in,:]>=zm[ilen-1) #shallower than ROMS
                    # lev[id_shallower] = ilen-1-1
                    # vrat[id_shallower] = 1 
                    
                    wild2[:,:,0,0]=temp[ix,iy,lev]*(1-vrat)+temp[ix,iy,lev+1]*vrat
                    wild2[:,:,0,1]=salt[ix,iy,lev]*(1-vrat)+salt[ix,iy,lev+1]*vrat
                    wild2[:,:,1,0]=temp[ix+1,iy,lev]*(1-vrat)+temp[ix+1,iy,lev+1]*vrat
                    wild2[:,:,1,1]=salt[ix+1,iy,lev]*(1-vrat)+salt[ix+1,iy,lev+1]*vrat
                    wild2[:,:,2,0]=temp[ix+1,iy+1,lev]*(1-vrat)+temp[ix+1,iy+1,lev+1]*vrat
                    wild2[:,:,2,1]=salt[ix+1,iy+1,lev]*(1-vrat)+salt[ix+1,iy+1,lev+1]*vrat
                    wild2[:,:,3,0]=temp[ix,iy+1,lev]*(1-vrat)+temp[ix,iy+1,lev+1]*vrat
                    wild2[:,:,3,1]=salt[ix,iy+1,lev]*(1-vrat)+salt[ix,iy+1,lev+1]*vrat
                                        
                    tempout[:,i_in[intri1]] = wild2[:,intri1,0,0]*arco[0,i_in_1] + \
                        wild2[:,intri1,1,0]*arco[1,i_in_1] + \
                        wild2[:,intri1,2,0]*arco[2,i_in_1]
                    
                    saltout[:,i_in[intri1]] = wild2[:,intri1,0,1]*arco[0,i_in_1] + \
                        wild2[:,intri1,1,1]*arco[1,i_in_1] + \
                        wild2[:,intri1,2,1]*arco[2,i_in_1]
                        
                                    #             if intri==1:
                    #                 tempout[k,i] = wild2[0,0]*arco[0,i] + \
                    #                               wild2[1,0]*arco[1,i] + \
                    #                               wild2[2,0]*arco[2,i]
                    #                 saltout[k,i] = wild2[0,1]*arco[0,i] + \
                    #                               wild2[1,1]*arco[1,i] + \
                    #                               wild2[2,1]*arco[2,i]
                    #             else:
                    #                 tempout[k,i] = wild2[0,0]*arco[0,i] + \
                    #                               wild2[2,0]*arco[1,i] + \
                    #                               wild2[3,0]*arco[2,i]
                    #                 saltout[k,i] = wild2[0,1]*arco[0,i] + \
                    #                               wild2[2,1]*arco[1,i] + \
                    #                               wild2[3,1]*arco[2,i]                        
                   
                    tempout[:,i_in[intri2]] = wild2[:,intri2,0,0]*arco[0,i_in_2] + \
                        wild2[:,intri2,2,0]*arco[1,i_in_2] + \
                        wild2[:,intri2,3,0]*arco[2,i_in_2]
                    
                    saltout[:,i_in[intri2]] = wild2[:,intri2,0,1]*arco[0,i_in_2] + \
                        wild2[:,intri2,2,1]*arco[1,i_in_2] + \
                        wild2[:,intri2,3,1]*arco[2,i_in_2]   
                    
                    #print("time4=%f"%(timer.time() - start))
                    
                    #Correct near surface T bias
                    idST = np.where( (kbp[ix,iy]!=-1) & (self.z[i_in,:].T>-10))
                    tempout[idST[0],i_in[idST[1]]]=tempout[idST[0],i_in[idST[1]]]-1
    
                    #if i==23293 and k == self.nvrt-1:
                    # if npout==0 and k == self.nvrt-1 and irecout2==100:
                    #     import pdb
                    #     pdb.set_trace()
                    # Check
                    # if (tempout[k,i]<tempmin or tempout[k,i]>tempmax or \
                    #     saltout[k,i]<saltmin or saltout[k,i]>saltmax):
                    #     raise("Interpolated values invalid:",i,k,
                    #       tempout[k,i],saltout[k,i])  
                    
                    # for i in range(self.nnode):
                    # #for i in nodes:
                    #     if np.isnan(ixy[i,0]) or np.isnan(ixy[i,1]):
                    #         #print("Cannot find a parent element:",i)
                    #         tempout[:,i]=tem_outside
                    #         saltout[:,i]=sal_outside                    
                    #         if weights[i]!=0.:
                    #             print("This node is not covered by the data", i)            
                    #     else:        
                    #         npout+=1
                    #         imap[npout]=i
                    #         ix=int(ixy[i,0])
                    #         iy=int(ixy[i,1])
                    #         intri=int(ixy[i,2])
                    #         #Find vertical level
                    #         tempout[:,i]=-999
                            
                    #         for k in range(self.nvrt):
                    #             if(kbp[ix,iy]==-1): #ROMS dry cell
                    #                 lev=ilen-1-1  #
                    #                 vrat=1
                    #             elif self.z[i,k]<=zm[int(kbp[ix,iy])]:
                    #                 lev=int(kbp[ix,iy])-1
                    #                 vrat=0
                    #             elif self.z[i,k]>=zm[ilen-1]:
                    #                 lev=ilen-1-1
                    #                 vrat=1
                    #             else:
                    #                 lev=-99 #flag
                    #                 for kk in range(ilen-1):
                    #                     if (self.z[i,k]>=zm[kk] and 
                    #                         self.z[i,k]<=zm[kk+1]):
                    #                         lev=kk
                    #                         vrat=(zm[kk]-self.z[i,k])/(zm[kk]-zm[kk+1])
                    #                         break                              
                    #             if lev==-99:
                    #                 raise('Cannot find a level:', i,k, self.z[i,k],
                    #                       zm)
                                    
                    #             wild2 = np.zeros((4,2))
                    #             wild2[0,0]=temp[ix,iy,lev]*(1-vrat)+temp[ix,iy,lev+1]*vrat
                    #             wild2[0,1]=salt[ix,iy,lev]*(1-vrat)+salt[ix,iy,lev+1]*vrat
                    #             wild2[1,0]=temp[ix+1,iy,lev]*(1-vrat)+temp[ix+1,iy,lev+1]*vrat
                    #             wild2[1,1]=salt[ix+1,iy,lev]*(1-vrat)+salt[ix+1,iy,lev+1]*vrat
                    #             wild2[2,0]=temp[ix+1,iy+1,lev]*(1-vrat)+temp[ix+1,iy+1,lev+1]*vrat
                    #             wild2[2,1]=salt[ix+1,iy+1,lev]*(1-vrat)+salt[ix+1,iy+1,lev+1]*vrat
                    #             wild2[3,0]=temp[ix,iy+1,lev]*(1-vrat)+temp[ix,iy+1,lev+1]*vrat
                    #             wild2[3,1]=salt[ix,iy+1,lev]*(1-vrat)+salt[ix,iy+1,lev+1]*vrat
                
                    #             if intri==1:
                    #                 tempout[k,i] = wild2[0,0]*arco[0,i] + \
                    #                               wild2[1,0]*arco[1,i] + \
                    #                               wild2[2,0]*arco[2,i]
                    #                 saltout[k,i] = wild2[0,1]*arco[0,i] + \
                    #                               wild2[1,1]*arco[1,i] + \
                    #                               wild2[2,1]*arco[2,i]
                    #             else:
                    #                 tempout[k,i] = wild2[0,0]*arco[0,i] + \
                    #                               wild2[2,0]*arco[1,i] + \
                    #                               wild2[3,0]*arco[2,i]
                    #                 saltout[k,i] = wild2[0,1]*arco[0,i] + \
                    #                               wild2[2,1]*arco[1,i] + \
                    #                               wild2[3,1]*arco[2,i]
                    #             #if i==23293 and k == self.nvrt-1:
                    #             # if npout==0 and k == self.nvrt-1 and irecout2==100:
                    #             #     import pdb
                    #             #     pdb.set_trace()
                    #             # Check
                    #             # if (tempout[k,i]<tempmin or tempout[k,i]>tempmax or \
                    #             #     saltout[k,i]<saltmin or saltout[k,i]>saltmax):
                    #             #     raise("Interpolated values invalid:",i,k,
                    #             #       tempout[k,i],saltout[k,i])                         
                
                    #             #Correct near surface T bias
                    #             if kbp[ix,iy]!=-1 and self.z[i,k]>-10:
                    #                 tempout[k,i]=tempout[k,i]-1
                
                    #             #Enforce lower bound for temp. for eqstate
                    #             tempout[k,i]=max(0.,tempout[k,i])  
                    print("applying spatial interpolation!")          
                    print("outputting at day ",dt.total_seconds()/86400, npout) 
                    
                    # only save the nodes with valid values. 
                    tempout_in = [temp_t[imap[:npout]] for temp_t in tempout]
                    saltout_in = [salt_t[imap[:npout]] for salt_t in saltout]
                    temperature.append(tempout_in)
                    salinity.append(saltout_in)
                    output_day.append(dt.total_seconds()/86400)
                    # if output_day[-1] == 350.625:
                    #     import pdb
                    #     pdb.set_trace()
                    nudge_step = int(pd.Timedelta(self.nudge_step).total_seconds()/3600)
                    if len(output_day)>1:
                        if output_day[-1]-output_day[-2]>(nudge_step+0.1)/24:
                            print(f"Current time step {output_day[-1]} and previous {output_day[-2]} file{ncfile}")
                            raise ValueError('The difference between current and previous time step is greater than assigned stride')
                    #print("time5=%f"%(timer.time() - start))
        # return weights, output_day, [temperature, salinity], \
        #     imap[:npout].astype(int)+1   #schism is one-based  
        temperature = np.array(temperature)
        salinity = np.array(salinity)
        temperature = np.transpose(temperature,(0,2,1)) # [var,time,node_map, nvrt]
        salinity = np.transpose(salinity,(0,2,1))        
        #Enforce lower bound for temp. for eqstate
        temperature[temperature<0]== 0 
        print("reorganizing matrix!")  
        #print("time=%f"%(timer.time() - start))
        return weights, [temperature, salinity], \
            imap[:npout].astype(int), output_day   #schism is one-based 

    
    @staticmethod
    def signa(x1,x2,x3,y1,y2,y3):
        signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
        return signa
    
    @staticmethod
    def reorder_variables(var_info,var_ordered_list):
        #reorder observed variables according to [temperature,salinity]
        #var_ordered_list = np.array(['temperature','salinity','elevation'])  # this list can be extended if new nudging variables are required.
        idx = [np.where(var_ordered_list == v['name'])[0][0] for v in var_info]
        return [var_info[i] for i in idx]
       
    def gen_nudge_obs(self, region_info):
        cutoff = 1e3  #when weights are less than cutoff, they will be set to zero
        weights, obs_df = self.gen_region_weight(region_info['attribute'],
                                                 region_info['vertices'])
        if len(obs_df)>1:
            obs_sites = obs_df.index.values
        weights[weights<weights.max()/cutoff] = 0
        method = region_info['interpolant']['method']        
        variable_info = region_info['interpolant']['variables'] 
        weights_list = []
        imap_list = []
        values_list = []
        for v in variable_info:
            empty_data = False
            none_values = v['none_values']
            obs = self.read_data(v['data'])
            name = v['name']                                     
            if name in obs.keys(): # this is in case there are multiple variables listed in the datafile
                vdata = obs[name]
            else: # only one variable listed in the file (but can be at multiple locations)
                vdata = obs            
            
            #if np.any(np.isnan(vdata)):  
            if vdata.isnull().any(axis=None):
                if isinstance(vdata, pd.core.series.Series): # pandas array
                    # if none_values == 'interpolate': #other fillna options possible.                        
                    #     vdata = vdata.interpolate()    
                    # elif none_values == 'ambient':  #fillna with -9999. 
                    #     vdata = vdata.fillna(-9999.)
                    # elif none_values == 'ignore': #ignore the entire series when there is any none
                    #     empty_data = True           
                    if none_values == 'error':
                        raise ValueError("nan values detected for %s in %s"%(v,v['data']))
                    elif none_values == 'ambient':
                        vdata = vdata.fillna(-9999.)
                    else:
                        raise NotImplementedError("Fill nan option %s is not implemented"%none_values)
                if isinstance(vdata,pd.core.frame.DataFrame): # multiple obs data,typically multiple stations
                    vdata = vdata.dropna(axis=1,how='all')
                    # if none_values == 'interpolate':
                    #     vdata = vdata.interpolate(axis=0)
                    # elif none_values = 'ambient':
                    #     vdata = vdata.fillna(-9999.)
                    # elif non_values == 'ignore':
                    #     vdata = vdata.dropna(axis=1,how='any') 
                    # if len(vdata.column)==0:
                    #     empty_data = True
                    if none_values == 'error':
                        raise ValueError("nan values detected for %s in %s"%(v,v['data']))
                    elif none_values == 'ambient':
                        vdata = vdata.fillna(-9999.)
                    else:
                        raise NotImplementedError("Fill nan option %s is not implemented"%none_values)
                if isinstance(vdata, xr.core.dataarray.DataArray): #xarray
                    vdata = vdata.dropna(dim='site',how='all')
                    # if none_values == 'interpolate':
                    #     vdata = vdata.interpolate_na(dim='time')
                    if none_values == 'error':
                        raise ValueError("nan values detected for %s in %s"%(v,v['data']))
                    elif none_values == 'ambient':
                        vdata = vdata.fillna(-9999.)
                    # elif none_values == 'ignore':
                    #     vdata = vdata.dropna(dim='site',how='any') 
                    if len(vdata.site)==0:
                        empty_data = True
                elif isinstance(vdata,xr.core.dataarray.Dataset): # the only possibility this happens is that the variable is not in the dataset
                    empty_data = True
            
            if empty_data:
                weights_v = np.zeros(self.nnode)
                imap_v = np.array([]) # empty array
                values_v = np.array([])
                weights_list.append(weights_v)
                values_list.append(values_v)
                imap_list.append(imap_v) 
                continue     
            
            if isinstance(vdata,xr.core.dataarray.DataArray):
                vdata_sites = vdata.site
                # the weights is ordered accodring to obs_sites;
                # and it will need to be reordered according to vdata.site
                w_list = []
                for s in vdata.site.values:
                    idx = np.where(obs_sites==s)[0]
                    w_list.append(weights[idx,:])                        
                weights_v = np.squeeze(np.array(w_list).sum(axis=0)) 
            elif isinstance(vdata,pd.core.frame.DataFrame):
                vdata_sites = vdata.columns
                # the weights is ordered accodring to obs_sites;
                # and it will need to be reordered according to vdata.site
                w_list = []
                for s in vdata.columns:
                    idx = np.where(obs_sites==s)[0]
                    w_list.append(weights[idx,:])
                weights_v = np.squeeze(np.array(w_list).sum(axis=0)) 
            else: # no adjustment for weights for a single obs site
                vdata_sites = ['single_site']
                weights_v = weights
            
            imap_v = np.where(weights_v>0)[0]            
            if (len(obs_df)>1) & (len(vdata_sites)>1): # multiple sites    
                print(vdata_sites)
                print(obs_df)                
                obs_x = obs_df.loc[vdata_sites].x.values.astype(float)
                obs_y = obs_df.loc[vdata_sites].y.values.astype(float)                
                                     
                if method =='nearest':
                    nn_id = []
                    for nx, ny in zip(self.node_x[imap_v],
                                      self.node_y[imap_v]):
                        dist = (nx-obs_x)**2 + (ny-obs_y)**2
                        nn_id.append(np.where(dist==dist.min())[0][0])
                    values_v = vdata.isel(site=nn_id).values    
                elif method == 'inverse_distance':
                    obs_loc = np.array([obs_x,obs_y]).T
                    values_v = []
                    if isinstance(vdata,xr.core.dataarray.DataArray):
                        for t in vdata.indexes['time']:
                            vals = vdata.sel(time=t).values
                            vals[vals==-9999.0] = np.nan
                            obs_loc_t = obs_loc[~np.isnan(vals)] # removing nan points
                            vals = vals[~np.isnan(vals)]
                            if (vals<0).any():
                                raise Exception("negative values detected in %s for %s: "%(v['data'],name)+ 
                                                 vals[vals<0])
                            # removing all the nans.                            
                            invdisttree = Interp2D.Invdisttree(obs_loc_t, vals,
                                           leafsize=10, stat=1)
                            node_xy = np.array([self.node_x[imap_v],
                                                self.node_y[imap_v]]).T
                            values_v.append(invdisttree(node_xy, nnear=4, p=2))

                    else:    
                        for t in vdata.index:
                            vals = vdata.loc[t].values
                            vals[vals==-9999.0] = np.nan
                            obs_loc_t = obs_loc[~np.isnan(vals)] # removing nan points
                            vals = vals[~np.isnan(vals)]     
                            if (vals<0).any():
                                raise Exception("negative values detected in %s for %s: "%(v['data'],name)+ 
                                                 vals[vals<0])
                            invdisttree = Interp2D.Invdisttree(obs_loc_t, vals,
                                           leafsize=10, stat=1)
                            node_xy = np.array([self.node_x[imap_v],
                                                self.node_y[imap_v]]).T
                            values_v.append(invdisttree(node_xy, nnear=4, p=2))
                else:
                     raise NotImplementedError                  
            else: # single site
                imap_v = np.where(weights_v>0)[0]
                values_v = np.broadcast_to(vdata.values,
                           (len(imap_v),
                            len(vdata))).T                 
                
            weights_list.append(weights_v)
            values_list.append(values_v)
            imap_list.append(imap_v)          
                
        # for v in variable_info:
        #     if v in obs.keys():
        #         vdata = obs[v]  
        #         if region_info['type'] == 'single_obs':                               
        #             if np.any(np.isnan(vdata)):  
        #                 if none_values == 'interpolate': #other fillna options possible.                        
        #                     if isinstance(vdata, pd.core.series.Series): # pandas array
        #                         vdata = vdata.interpolate()
        #                     elif isinstance(vdata, xr.core.dataarray.DataArray): #xarray
        #                         vdata = vdata.interpolate_na(dim='time')
        #                     weights_v = weights
        #                     imap_v = np.where(weights_v>0)[0]
        #                     values_v = np.broadcast_to(vdata.values,
        #                                                (len(imap_v),
        #                                                 len(vdata))).T
        #                 elif none_values == 'ignore': #ignore the entire series when there is any none
        #                     weights_v = np.zeros(weights)
        #                     imap_v = np.array([]) # empty array
        #                     values_v = np.array([])
        #             else:
        #                 weights_v = weights
        #                 imap_v = np.where(weights_v>0)[0]
        #                 values_v = np.broadcast_to(vdata.values,
        #                                            (len(imap_v),
        #                                             len(vdata))).T  
        #         elif region_info['type'] == 'multi_obs':  #multiple observational points
        #             # multiple input should be stored in netcdf format only. 
        #             if np.any(np.isnan(vdata)): 
        #                 vdata = vdata.dropna(dim='site',how='all')
        #                 if none_values == 'interpolate':
        #                     vdata = vdata.interpolate_na(dim='time')
        #                 elif none_values == 'ignore':
        #                     vdata = vdata.dropna(dim='site',how='any')
                            
        #                 inc_site = []
        #                 for s in obs.site.values:
        #                     if s in vdata.site.values:
        #                         inc_site.append(True)
        #                     else:
        #                         inc_site.append(False)
        #                 weights_v = weights[inc_site,:].sum(axis=0)
        #             else:
        #                 weights_v = weights
        #             imap_v = np.where(weights_v>0)[0]                
        #             obs_x = obs.x.sel(site=vdata.site).values.astype(float)
        #             obs_y = obs.y.sel(site=vdata.site).values.astype(float)
                                     
        #             if method =='nearest':
        #                 nn_id = []
        #                 for nx, ny in zip(self.node_x[imap_v],
        #                                   self.node_y[imap_v]):
        #                     dist = (nx-obs_x)**2 + (ny-obs_y)**2
        #                     nn_id.append(np.where(dist==dist.min())[0][0])
        #                 values_v = vdata.isel(site=nn_id).values    
        #             elif method == 'inverse_distance':
        #                 obs_loc = np.array([obs_x,obs_y]).T
        #                 values_v = []
        #                 for t in vdata.time:
        #                     vals = vdata.sel(time=t).values
        #                     invdisttree = Interp2D.Invdisttree(obs_loc, vals,
        #                                    leafsize=10, stat=1)
        #                     node_xy = np.array([self.node_x[imap_v],
        #                                         self.node_y[imap_v]]).T
        #                     values_v.append(invdisttree(node_xy, nnear=4, p=2))
        #             else:
        #                  raise NotImplementedError   
        #     else:
        #         weights_v = []
        #         values_v = []
        #         imap_v = []
        #     weights_list.append(weights_v)
        #     values_list.append(values_v)
        #     imap_list.append(imap_v) 
     
        return weights_list, values_list, imap_list          

    def read_data(self, data):
        if data.endswith('csv'):
            obs = pd.read_csv(data,index_col='datetime',parse_dates=['datetime'])
            obs.index.name = 'time'   
            obs = obs[(obs.index>=
                       pd.to_datetime(self.start_date)) &
                      (obs.index<=
                       pd.to_datetime(self.end_date))]
            # time interpolation
            obs = obs.resample(self.nudge_step).nearest()
            
            if len(obs)< len(self.datetime):
                raise ValueError('The input time series may not cover the \
                      entire nudging period: check %s'%data)
            else:
                obs = obs.reindex(self.datetime)
        elif data.endswith('nc'):
            obs = xr.open_dataset(data)
            obs = obs.sel(time=self.datetime)     
        else:
            raise NotImplementedError()
        return obs        

    def gen_region_weight(self, attribute, vertices):
        if isinstance(attribute, str):
            if vertices:
                inpoly = self.in_vertices(vertices)
            weights = np.zeros_like(self.node_x)
            for i, (x, y) in enumerate(zip(self.node_x[inpoly], 
                                           self.node_y[inpoly])):
                weights[inpoly[i]] = eval(attribute)
            obs_df = []
        else:
            if isinstance(attribute['xy'],str): # multiple points
                if attribute['xy'].endswith('nc'):
                    dataset = xr.open_dataset(attribute['xy'])
                    x0 = dataset['x'].values.astype(float)
                    y0 = dataset['y'].values.astype(float)
                    obs_sites = dataset['site'].values
                elif (attribute['xy']).endswith('csv'):
                    dataset = pd.read_csv(attribute['xy'])
                    x0 = dataset['x'].values
                    y0 = dataset['y'].values
                    obs_sites = dataset['site'].values
                obs_df = pd.DataFrame({'site':obs_sites,
                                       'x':x0,
                                       'y':y0})
                obs_df = obs_df.set_index('site')
            else:
                x0,y0 = attribute['xy']
                obs_df = []
                
            if isinstance(x0,np.ndarray):
                weights = np.zeros([len(x0),len(self.node_x)])
                norm_factor = np.ones_like(x0)
                for i, (xc, yc) in enumerate(zip(x0, y0)):
                    weights_c = self.gaussian_weights([xc,yc],[x0,y0],attribute)
                    # perhaps other more sophisticated normalization can be used here instead
                    if sum(weights_c)>1:
                        norm_factor[i] = weights_c[i]/sum(weights_c) # if 1 exceeds
            else:
                weights = np.zeros_like(self.node_x)
                norm_factor = 1
            if vertices != 'None':
                inpoly = self.in_vertices(vertices)
                for i, (x, y) in enumerate(zip(self.node_x[inpoly], 
                                               self.node_y[inpoly])):  
                    weights[inpoly[i]] = self.construct_weights(attribute,x,y,
                                                                x0,y0,norm_factor)
            else:
                for i, (x, y) in enumerate(zip(self.node_x, 
                                               self.node_y)):
                    weights[:,i] = self.construct_weights(attribute,x,y,x0,y0,
                                                          norm_factor)    
        
        return weights, obs_df

    def construct_weights(self, attribute,x,y,x0,y0,norm_factor):
        if isinstance(attribute, str):
            return eval(attribute)
        elif isinstance(attribute, dict):
            if attribute['kernel'] == 'gaussian':
                if isinstance(x0,np.ndarray): # multiple points 
                    weights = self.gaussian_weights([x,y],[x0,y0],attribute)
                    weights = weights*norm_factor/pd.Timedelta(
                        attribute['time_scale']).total_seconds()  
                    return weights # multiple weights
                else: # float or int
                    weight = self.gaussian_weights([x,y],
                                                   [x0,y0],attribute)
                    weight = weight/pd.Timedelta(
                        attribute['time_scale']).total_seconds()
                    return weight # single weight          
            else:
                raise NotImplementedError(
                    "%s kernel not implemented"%attribute['kernel'])
                
    def plot(self,imap,values,**kwargs):
        v = np.zeros(self.nnode)
        v[imap] = values
        col = self.mesh.plot_nodes(v,**kwargs)
        return col
    
    @staticmethod
    def gaussian_weights(xy,xyc,attribute):
        r2 = (xy[0]-xyc[0])**2 + (xy[1]-xyc[1])**2  #can be multiple points
        # calculate weights however normalization is needed
        weight = np.exp(-1*r2/2/attribute['length_scale']**2)
        return weight

    def in_vertices(self, vertices):        
        Polygon(vertices)
        p = Polygon(vertices)        
        inpoly = np.where(self.mesh_gpd.within(p))[0]
        return inpoly  


def create_arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Create hotstart for a schism run")
    parser.add_argument('--yaml_fn',type=str, help='yaml file for nudging',required=True)
    parser.add_argument('--suffix',type=str, help='suffix for generated nudging files', default=None)
    parser.add_argument('--crs',type=str,help='The projection system for the mesh',default=None)
    return parser
    
def main():
    # User inputs override the yaml file inputs.  
    parser = create_arg_parser() 
    args = parser.parse_args()
    
    nudging = nudging(args.yaml_fn,crs=args.crs,
                                     suffix=suffix)
    nudging.read_yaml()
    nudging.create_nudging()
           
if __name__ == "__main__":
    main()