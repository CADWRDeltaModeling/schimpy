# -*- coding: utf-8 -*-
"""
Module for creating a schism hotstart file

Options:
1. Centering:
    * node2D: node (2D surface)
    * node3D: node at whole level
    * edge: edge center at whole level (mainly for horizontal velocty)
    * elem: element center at whole level (mainly for vertical velocity)
    * prism: prism center (for tracers); three different variables tr_el, tr_nd,
    * tr_el need to be generated for the hotstart.nc file.
        *  tr_nd = tr_nd0 (on nodes)
        *  tr_el (on element at middle level)
2. Initializer:
    * text_init: 2D map input from either *.prop or *.ic files; require input text file.
    * simple_trend: can either be a number or an equation that depends on x (lat) and y(lon). e.g,: 0.97 + 1.e-5*x
    * patch_init: regional-based method. Polygon shapefile input required.
    * extrude_casts: 3D nearest neighbourhood interpolation based on Polaris cruise transects.
    * obs_points: interpolated from 1D (time series of) observational points.
    * hotstart_nc: initialize with a previous hotstart output from the same or a different mesh.
    * schout_nc: initialize with a schism output snapshot.

Required data files:
    * param.in
    * hgrid.gr3
    * vgrid.in
    * hotstart.yaml and input files defined in the yaml file.
"""

import yaml
import numpy as np
import xarray as xr
import pandas as pd
from scipy import interpolate
from scipy.spatial import distance
import matplotlib.pyplot as plt
import os
# import from schimpy

from schimpy.schism_mesh import read_mesh, write_mesh, SchismMeshGr3Reader, \
    compare_mesh
from schimpy import geo_tools

class SCHISMHotstart(object):
    """ Mock object """
    pass

def read_hotstart():
    """Mock method to read existing hotstart"""
    pass

class hotstart(object):
    """
    A class to generate hotstart initial condition for SCHISM
    """

    # will change these into yaml file on
    def __init__(self, yaml_fn=None,modules=None,crs=None):
        # read input from yaml files;
        # create elev.ic if it does not exist or not used as an input
        # the modules to turn on.
        self.yaml_fn = yaml_fn
        self.nc_dataset = None
        self.modules = modules
        self.crs = crs
        # if 'HYDRO' in modules:
        #     modules.remove('HYDRO')
        #     modules+=['TEM', 'SAL']
        # self.modules = modules # modules can be none for barotropic runs.
        # if 'param_nml' in kwargs:
        #     self.param_nml = kwargs['param_nml']
        # else:
        #     if any([m in ['SED', 'AGE', 'GEN', 'ECO'] for m in modules]): # these modules require param.nml
        #         raise ValueError('param_nml needs to be defined in %s'%yaml_fn)
        #     else: # for all other modules, param.nml will not be used.
        #         self.param_nml = "param.nml"
        # if crs is None:
        #     self.crs = 'EPSG:26910'  # when crs is not specified, use default UTM10N projection. The projection is not important if all inputs are in the same coordinate system.
        # #    raise ValueError("crs must be specified")
        # else:
        #     print("crs is {}".format(crs))
        #     self.crs = crs

        # if self.modules:  # if modules is None for barotropic run, this step is not required.
        #     self.ntracers, self.ntrs, self.irange_tr, self.tr_mname = \
        #         describe_tracers(self.param_nml, modules=self.modules)
        # else:
        #     self.ntracers = 0
        #     self.tr_mname = []

    def read_yaml(self):
        """
        read yaml and load mesh grid
        """
        with open(self.yaml_fn) as file:
            info = yaml.load(file, Loader=yaml.FullLoader)
        hotstart_info = info['hotstart']
        self.info = hotstart_info
        self.date = hotstart_info['date']
        self.hgrid_fn = hotstart_info['hgrid_input_file']
        self.vgrid_fn = hotstart_info['vgrid_input_file']
        variables = list(hotstart_info.keys())
        # remove these items from the list.
        rfl = ['hgrid_input_file', 'vgrid_input_file', 'date','vgrid_version',
               'crs','param_nml','modules','output_fn']
        variables = [v for v in variables if v not in rfl]
        self.variables = variables
        self.vgrid_version = hotstart_info['vgrid_version']  # 5.8 and below is the old version

        if 'output_fn' in hotstart_info.keys():
            self.output_fn = hotstart_info['output_fn']
        else:
            self.output_fn = 'hotstart.nc'

        if self.modules is None:
            if 'modules' not in hotstart_info.keys():
                modules = 'HYDRO' # hydro is the minimum required module
            else:
                modules = hotstart_info['modules']
            if 'HYDRO' in modules:
                modules.remove('HYDRO')
                modules+=['TEM', 'SAL']
            self.modules = modules # modules can be none for barotropic runs.

        if 'param_nml' in hotstart_info.keys():
            self.param_nml = hotstart_info['param_nml']
        else:
            if self.modules is not None:
                if any([m in ['SED', 'AGE', 'GEN', 'ECO'] for m in self.modules]): # these modules require param.nml
                    raise ValueError('param_nml needs to be defined in %s'%self.yaml_fn)
            self.param_nml = "param.nml" # for all other modules, param.nml will not be used.

        if self.crs is None:
            if 'crs' not in hotstart_info.keys():
                self.crs = 'EPSG:26910'  # when crs is not specified, use default UTM10N projection. The projection is not important if all inputs are in the same coordinate system.
            #    raise ValueError("crs must be specified")
            else:
                crs = hotstart_info['crs']
                print("crs is {}".format(crs))
                self.crs = crs

        if self.modules:  # if modules is None for barotropic run, this step is not required.
            self.ntracers, self.ntrs, self.irange_tr, self.tr_mname = \
                describe_tracers(self.param_nml, modules=self.modules)
        else:
            self.ntracers = 0
            self.tr_mname = []

    def create_hotstart(self):
        self.read_yaml()
        mesh = read_mesh(self.hgrid_fn, self.vgrid_fn, self.vgrid_version)
        self.mesh = mesh
        self.depths = self.mesh.build_z() # this is required to get depths for elevation
        #v = self.variables[1]
        self.hotstart_ini = {}
        for v in self.variables:
            print("creating hotstart for %s" % v)
            initializer = self.info[v]['initializer']

            if ('hotstart_nc' in initializer) and \
                (not self.hotstart_ini) :
                self.hotstart_ini['hotstart_nc_hfn'] = initializer[
                    'hotstart_nc']['source_hgrid']
                self.hotstart_ini['hotstart_nc_vfn'] = initializer[
                    'hotstart_nc']['source_vgrid']
                self.hotstart_ini['source_vgrid_version'] = initializer[
                    'hotstart_nc']['source_vgrid_version']
            elif 'patch_init' in initializer and \
                (not self.hotstart_ini) :
                patch_init =  [list(ini['initializer'].keys())[0]
                               for ini in initializer['patch_init']['regions']]
                if 'hotstart_nc' in patch_init:
                    idx_h = np.where(np.array(patch_init)=='hotstart_nc')[0][0]
                    sub_init = initializer['patch_init'][
                        'regions'][idx_h]['initializer']['hotstart_nc']
                    self.hotstart_ini['hotstart_nc_hfn'] = sub_init['hgrid']
                    self.hotstart_ini['hotstart_nc_vfn'] = sub_init['vgrid']
                    self.hotstart_ini['source_vgrid_version'] = sub_init[
                        'source_vgrid_version']

            if self.hotstart_ini:
                hotstart_mesh = read_mesh(self.hotstart_ini['hotstart_nc_hfn'],
                                          self.hotstart_ini['hotstart_nc_vfn'],
                                          self.hotstart_ini['source_vgrid_version'])
                indices, dist = compare_mesh(hotstart_mesh,self.mesh)
                self.hotstart_ini['hotstart_nc_indices'] = indices
                self.hotstart_ini['hotstart_nc_dist'] = dist

            var = self.generate_3D_field(v)
            if self.nc_dataset is None:
                if 'tke' in self.info.keys():
                    self.initialize_netcdf(default_turbulence=False)
                else:
                    self.initialize_netcdf()
            self.nc_dataset = self.nc_dataset.merge(var)
        #  cumsum_eta is required by the most recent version of schism.
        self.nc_dataset = self.nc_dataset.assign(
            cumsum_eta = self.nc_dataset.elevation)
        if 'COSINE' in self.modules: #COS_mS2=COS_5, COS_mDN=COS_8, COS_mZ1 = COS_6, COS_mZ2 = COS_7
            self.nc_dataset = self.nc_dataset.assign(
                COS_mS2 = self.nc_dataset.COS_5_el)
            self.nc_dataset = self.nc_dataset.assign(
                COS_mDN = self.nc_dataset.COS_8_el)
            self.nc_dataset = self.nc_dataset.assign(
                COS_mZ1 = self.nc_dataset.COS_6_el)
            self.nc_dataset = self.nc_dataset.assign(
                COS_mZ2 = self.nc_dataset.COS_7_el)
            # for COS_sxx, these means the sum values and will be set to zero
            # 'COS_sS2  ','COS_sDN  ','COS_sZ1  ','COS_sZ2  ','COS_nstep'
            self.nc_dataset['COS_sS2']=(
                ('elem','nVert'),np.zeros_like(self.nc_dataset.COS_1_el.values))
            self.nc_dataset['COS_sDN']=(
                ('elem','nVert'),np.zeros_like(self.nc_dataset.COS_1_el.values))
            self.nc_dataset['COS_sZ1']=(
                ('elem','nVert'),np.zeros_like(self.nc_dataset.COS_1_el.values))
            self.nc_dataset['COS_sZ2']=(
                ('elem','nVert'),np.zeros_like(self.nc_dataset.COS_1_el.values))
            self.nc_dataset['COS_nstep']=(
                ('elem','nVert'),np.zeros_like(self.nc_dataset.COS_1_el.values))
        # correct wet/dry cells depending on water level (eta2 value):mostly not needed.
        self.depth = self.mesh.build_z(elev=self.nc_dataset['elevation'].values) # compute depth again for the updated elevation.
        self.wet_dry_check()
        self.map_to_schism()
        return self.nc_dataset

    def generate_3D_field(self, variable):
        v_meta = self.info[variable]
        kwargs_options = {}
        if 'SED' in variable:
            if 'Nbed' not in self.__dict__.keys():
                if os.path.isfile("sediment.nml"):
                    self.sediment_fn = "sediment.nml"
                    params = read_param_nml(self.sediment_fn)
                    self.Nbed = params['Nbed']
                else:
                    raise FileNotFoundError(
                        "sediment.nml is required if SED module is turned on")
            kwargs_options.update({'Nbed': self.Nbed})
        if kwargs_options:
            var = VariableField(v_meta, variable, self.mesh,
                                self.depths, self.date, self.crs,
                                self.tr_mname, **kwargs_options)
        else:
            var = VariableField(v_meta, variable, self.mesh,
                                self.depths, self.date, self.crs,
                                self.tr_mname,self.hotstart_ini)
        return var.GenerateField()

    def initialize_netcdf(self, default_turbulence=True):
        if not self.nc_dataset:  # if the dataset is empty, initialize the nc file
            #self.describe_tracers = [h.info[k]['centering'] for k in h.variables].count('prism_c')
            idry_e = np.zeros(self.mesh.n_elems()).astype(int)
            idry = np.zeros(self.mesh.n_nodes()).astype(int)
            idry_s = np.zeros(self.mesh.n_edges()).astype(int)

            ds = xr.Dataset({'time': ('one', [0.0]),
                             'iths': ('one', [0]),
                             'ifile': ('one', [1]),
                             'idry_e': ('elem', idry_e),
                             'idry_s': ('side', idry_s),
                             'idry': ('node', idry),
                             'nsteps_from_cold':('one',[0])})
            # coords={'one':range(1),
            # 'elem':range(self.mesh.n_elems()),
            # 'node':range(self.mesh.n_nodes()),
            # 'side':range(self.mesh.n_edges())})
            # Now apply default setting for  q2, xl, dfv,dfh,dfq1,dfq2,
            # SED3D_dp, SED3D_rough, SED3D_bed, SED3D_bedfrac
            if default_turbulence:
                # Turbulent kinetic energy
                ar_name = ['q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2']
                for name in ar_name:
                    ds_ar = xr.Dataset({name: (['node', 'nVert'],
                                               np.zeros((self.mesh.n_nodes(), self.mesh.n_vert_levels)))})
                    # coords={'node':range(self.mesh.n_nodes()),
                    # 'nVert':range(self.mesh.n_vert_levels)})
                    ds = ds.merge(ds_ar)
        self.nc_dataset = ds

    def wet_dry_check(self):
        """
        modify idry_e, idry, and idry_s based on eta2 and water depth.
        """
#        if hasattr(self.nc_dataset,'eta2'):
#            z = self.mesh.build_z(elev=self.nc_dataset['eta2'])
#        else:
#            z = self.mesh.build_z()
        self.nc_dataset['z'] = xr.DataArray(
            self.depths, dims=['node', 'nVert'])
        dry_nodes = np.where(self.depths.min(axis=1) >= 0)[0]
        if any(dry_nodes):
            idry = self.nc_dataset['idry'].values
            idry[dry_nodes] = 1  # dry cells have positive z values
            idry_s = np.array([idry[n].mean(axis=0)
                               for n in self.mesh.edges[:, :2]])  # node to side
            # if there is one dry node, the element is dry?
            idry_s[idry_s != 0] = 1
            idry_e = np.array([idry[list(n)].mean(axis=0)
                               for n in self.mesh.elems])  # node to elem
            # if there is one dry node, the side is dry?
            idry_e[idry_e != 0] = 1
            self.nc_dataset['idry_s'] = xr.DataArray(idry_s, dims=['side'])
            self.nc_dataset['idry'] = xr.DataArray(idry, dims=['node'])
            self.nc_dataset['idry_e'] = xr.DataArray(idry_e, dims=['elem'])

    # def compute_depths(self, elevation=None):
    #     """ Calculate points in 3-D grid
    #     """
    #     mesh = self.mesh
    #     vert_grid = mesh.vmesh
    #     nvrt = vert_grid.param['nvrt']
    #     depths = np.full([mesh.n_nodes(), nvrt], np.finfo(np.float32).min)

    #     if vert_grid.ivcor == 1:
    #         if nvrt < 3:
    #             raise ValueError()

    #         for i, node in enumerate(mesh.nodes):
    #             if elevation is None:
    #                 hmod2 = max(0.1, node[2])
    #             else:
    #                 hmod2 = max(0.1, node[2]+elevation[i])
    #             kbps = vert_grid.kbps[i]
    #             depths[i, kbps:] = hmod2 * vert_grid.sigma[i, :nvrt-kbps]
    #             depths[i, :kbps] = depths[i, kbps]
    #         self.depths = depths
    #     elif vert_grid.ivcor == 2:
    #         h_s = vert_grid.param['h_s']
    #         h_c = vert_grid.param['h_c']
    #         kz = vert_grid.param['kz']
    #         c_s = vert_grid.c_s
    #         for i, node in enumerate(mesh.nodes):
    #             # TODO: Maybe add surface elevation?
    #             if elevation is None:
    #                 depth = node[2]
    #             else:
    #                 depth = node[2]+elevation[i]
    #             # S-level
    #             hmod2 = max(0.1, min(depth, h_s))
    #             for k, s in enumerate(vert_grid.sigma):
    #                 k_a = kz + k - 1
    #                 if hmod2 <= h_c:
    #                     depths[i, k_a] = s * hmod2
    #                 else:
    #                     depths[i, k_a] = h_c * s + (hmod2 - h_c) * c_s[k]
    #             # Z-level
    #             for k, d in enumerate(vert_grid.ztot[:-1]):
    #                 depths[i, k] = d
    #         self.depths = depths
    #     else:
    #         raise("Not supported ivcor")

    def map_to_schism(self):
        """
        map the hotstart variable to the variable names required by SCHISM
        hotstart file.
        """
        param_nml = self.param_nml
        schism_hotstart_var = {'eta2': 'elevation',
                               'we': 'velocity_w',
                               'su2': 'velocity_u',
                               'sv2': 'velocity_v'}
        if self.modules:
            ntracers, ntrs, irange_tr, tr_mname = describe_tracers(
                 param_nml, modules=self.modules)
            tr_mname_el = ["%s_el" % n for n in tr_mname]
            tr_mname_nd = ["%s_nd" % n for n in tr_mname]
            tr_mname_nd0 = ["%s_nd0" % n for n in tr_mname]
            schism_hotstart_var['tr_el'] = tr_mname_el
            schism_hotstart_var['tr_nd'] = tr_mname_nd
            schism_hotstart_var['tr_nd0'] = tr_mname_nd0

        nc_dataset = self.nc_dataset
        #v_keys =list(nc_dataset.data_vars)
        if any([v in ['HYDRO','TEM','SAL'] for v in self.modules]):
            mapping_dict = {'temperature_nd': 'TEM_nd',
                            'temperature_nd0': 'TEM_nd0',
                            'temperature_el': 'TEM_el',
                            'salinity_nd': 'SAL_nd',
                            'salinity_nd0': 'SAL_nd0',
                            'salinity_el': 'SAL_el'}

            nc_dataset = nc_dataset.rename(mapping_dict)

        for key, var_values in zip(schism_hotstart_var.keys(),
                                   schism_hotstart_var.values()):
            if isinstance(var_values, str):
                nc_dataset = nc_dataset.rename({var_values: key})
            else:
                var_temp = np.asarray(
                    [nc_dataset[var_values][v].values for v in var_values])
                var_temp = np.transpose(var_temp, (1, 2, 0))
                var_dim = nc_dataset[var_values[0]].dims
                #dim_values = [nc_dataset[var_values[0]][dim].values for dim in var_dim]
                var_dim = var_dim + ('ntracers',)
                # dim_values.append(range(len(var_values)))
                #coords = {}
                # for vd,dv in zip(var_dim,dim_values):
                #coords[vd] = dv
                v1_sub = xr.Dataset({key: (var_dim, var_temp)})
                # coords=coords)
                nc_dataset = nc_dataset.merge(v1_sub)
                nc_dataset = nc_dataset.drop(var_values)

        if self.modules:
            nc_dataset['tracer_list'] = tr_mname
        self.nc_dataset = nc_dataset


class VariableField(object):
    """
    A class initializing each varaible for the hotstart file
    """

    def __init__(self, v_meta, vname, mesh, depths, date, crs, tr_mname,
                 hotstart_ini=None, **kwargs):
        #self.centering = v_meta['centering']
        self.variable_name = vname
        self.centering = self.var_centering() # the centering option is dependent on the varialbe name.
        self.mesh = mesh
        self.depths = depths
        if self.centering == 'edge':
            self.edge_depths = (depths[self.mesh.edges[:,0]] +
                                depths[self.mesh.edges[:,1]])/2
        elif self.centering == 'elem':
            self.elem_depths = [self.mesh.z[e].mean(axis=0) for e
                                in self.mesh.elems]
        self.date = date
        self.crs = crs
        self.n_edges = mesh.n_edges()
        self.n_nodes = mesh.n_nodes()
        self.n_elems = mesh.n_elems()
        self.n_vert_levels = mesh.n_vert_levels
        self.node = mesh.nodes[:, :2]
        self.node_z = mesh.nodes[:,2]
        self.edge = mesh.get_centers_of_sides()[:, :2]
        self.elem = mesh.get_centers_of_elements()
        self.elems = mesh.elems
        self.initializer = self.get_key(v_meta['initializer'])
        self.ini_meta = self.get_value(v_meta['initializer'])
        self.tr_mname = tr_mname
        self.tr_index = None
        self.hotstart_ini = hotstart_ini

        if vname in ['SED3D_bed', 'SED3D_bedfrac']:
            self.Nbed = kwargs['Nbed']
        else:
            self.Nbed = 1
        self.grid = self.get_grid()  # grid can be nodes/elems/edges
        self.n_hgrid = list(self.grid.values())[0][0]
        self.n_vgrid = list(self.grid.values())[1][0]
        self.hgrid = list(self.grid.values())[0][1]
        self.vgrid = list(self.grid.values())[1][1]
        self.hgrid_name = list(self.grid.keys())[0]
        self.vgrid_name = list(self.grid.keys())[1]

        if vname in ['SED3D_bed', 'SED3D_bedfrac']:
            self.n_sdim = list(self.grid.values())[2][0]  # sediment dimension
            self.sdim_name = list(self.grid.keys())[2]

    def var_centering(self):
        # built-in centering options for model input variables.
        tr_key_list = ['tem','sal','age_','gen_','icm','sed_','cos_']
        if any([t in self.variable_name.lower() for t in tr_key_list]):
            variable_type = 'tracer'
        else:
            variable_type = self.variable_name

        vc_mapping = {'elevation':'node2D',
                      'velocity_u':'edge',
                      'velocity_v':'edge',
                      'velocity_w':'elem',
                      'SED3D_dp':'node2D',
                      'SED3D_rough':'node2D',
                      'SED3D_bed':'bed',
                      'SED3D_bedfrac':'bedfrac',
                      'tracer':'prism',
                      'tke':'node3D'
                      }
        if variable_type in vc_mapping.keys():
            return vc_mapping[variable_type]
        else:
            raise ValueError("variable %s does not have a valid input centering option"%self.variable_name)
    def get_grid(self):
        """Getting the number and coordinates of horizontal nodes/elems/edges for
        different centering options and return them as grid.
        """
        options = ['node3D', 'node2D', 'edge',
                   'elem', 'prism', 'bed', 'bedfrac']
        if self.centering not in options:
            raise NotImplementedError('''centering option %s is not implemented
                                      availabel options are: ''' % self.centering +
                                      ",".join(options))
        centering_options = {'node3D': {'node': (self.n_nodes, self.node),
                                        'nVert': (self.n_vert_levels, self.depths)},

                             'prism': {'node': (self.n_nodes, self.node),        # on prism center: mainly for tracers; first average from nodes to centers and then then average to prism center
                                       'nVert': (self.n_vert_levels, self.depths)},
                             'node2D': {'node': (self.n_nodes, self.node),
                                        'surface': (1, self.node_z[:,np.newaxis])},
                             'bed': {'elem': (self.n_elems, self.elem),
                                     'Nbed': (self.Nbed, 0),
                                     'MBEDP': (3, ['layer thickness', 'layer age', 'layer porosity'])},
                             'bedfrac': {'elem': (self.n_elems, self.elem),
                                         'Nbed': (self.Nbed, 0),
                                         'nsed': (3, 0)}}

        if self.centering == 'edge':
            centering_options['edge'] = {'side': (self.n_edges, self.edge),
                                         'nVert': (self.n_vert_levels,
                                                   self.edge_depths)}
        elif self.centering == 'elem':
            centering_options['elem'] = {'elem': (self.n_elems, self.elem),         # on element but at whole level: mainly for w
                                         'nVert': (self.n_vert_levels,
                                                   self.elem_depths)}
        return centering_options[self.centering]

    def define_new_grid(self, mesh):
        """
        Generating hgrid and vgrid based on centering option for a new grid
                Parameters
        """
        options = ['node3D', 'node2D', 'edge',
                   'elem', 'prism', 'bed', 'bedfrac']
        if self.centering not in options:
            raise NotImplementedError('''centering option %s is not implemented
                                      availabel options are: ''' % self.centering +
                                      ",".join(options))
        #mesh = read_mesh(hgrid_fn,vgrid_fn)
        mesh_node = mesh.nodes[:, :2]
        mesh_edge = mesh.get_centers_of_sides()[:, :2]
        mesh_elem = mesh.get_centers_of_elements()
        mesh_depth = mesh.build_z()

        centering_options = {'node3D': {'node': (mesh.n_nodes, mesh_node),
                                        'nVert': (mesh.n_vert_levels, mesh_depth)},
                             'prism': {'node': (mesh.n_nodes, mesh_node),        # on prism center: mainly for tracers; first average from nodes to centers and then then average to prism center
                                       'nVert': (mesh.n_vert_levels, mesh_depth)},
                             'node2D': {'node': (mesh.n_nodes, mesh_node),
                                        'surface': (1, mesh_depth[:,0])},
                             'bed': {'elem': (mesh.n_elems, mesh_elem),
                                     'Nbed': (self.Nbed, 0),
                                     'MBEDP': (3, ['layer thickness', 'layer age', 'layer porosity'])},
                             'bedfrac': {'elem': (mesh.n_elems, mesh_elem),
                                         'Nbed': (self.Nbed, 0),
                                         'nsed': (3, 0)}}
        if self.centering == 'edge':
            edge_depths = (mesh_depth[mesh.edges[:,0]] +
                           mesh_depth[mesh.edges[:,1]])/2
            centering_options['edge'] = {'side': (mesh.n_edges, mesh_edge),
                                         'nVert': (mesh.n_vert_levels, edge_depths)}
        elif self.centering == 'elem':
            elem_depths = [mesh_depth[e].mean(axis=0) for e in mesh.elems]
            centering_options['elem'] = {'elem': (mesh.n_elems, mesh_elem),         # on element but at whole level: mainly for w
                                         'nVert': (mesh.n_vert_levels, elem_depths)}

        return centering_options[self.centering]

    def initializer_options(self):
        """Options for input file initializers
        The input initializer overwrites self.initializer
        """
        initializer = self.initializer
        options = ['text_init', 'simple_trend',
                   'patch_init', 'extrude_casts', 'obs_points',
                   'hotstart_nc','schout_nc']
        if initializer not in options:
            raise NotImplementedError('''initializer %s is not implimented
                                      available options are: ''' % initializer +
                                      ",".join(options))
        initializers = {'text_init': self.text_init,
                        'simple_trend': self.simple_trend,
                        'patch_init': self.patch_init,
                        'extrude_casts': self.extrude_casts,
                        'obs_points': self.obs_points,
                        'hotstart_nc': self.hotstart_nc,
                        'schout_nc': self.schout_nc}
        return initializers[self.initializer]

    def patch_initializer_options(self, initializer):
        """Options for patch initializers
        The input initializer overwrites self.initializer
        """
        options = ['text_init', 'simple_trend', 'extrude_casts', 'obs_points',
                   'hotstart_nc','schout_nc']
        if initializer not in options:
            raise NotImplementedError('''initializer %s is not implimented
                                      available options are: ''' % initializer +
                                      ",".join(options))
        initializers = {'text_init': self.text_init,
                        'simple_trend': self.simple_trend,
                        'extrude_casts': self.extrude_casts,
                        'obs_points': self.obs_points,
                        'hotstart_nc': self.hotstart_nc,
                        'schout_nc': self.schout_nc}
        return initializers[initializer]

    def GenerateField(self):
        #initializer = globals()[self.initializer_options()]()
        initializer = self.initializer_options()
        var = initializer()
        #if np.all(np.isnan(var)):
        #    raise ValueError(
        #        "%s field has nan value everywhere: check input data and variable name" % self.variable_name)
        da = self.create_dataarray(var)
        return da

    @classmethod
    def map_to_3D(cls, v_2d, nz):
        """
        polulate the 2D map values to 3D by applying the uniform values over dpeth.
        """
        v_2d = np.asarray(v_2d)
        return np.tile(v_2d[..., np.newaxis], (1, nz))

    @classmethod
    def get_value(cls, dict_obj):
        return list(dict_obj.values())[0]

    @classmethod
    def get_key(cls, dict_obj):
        return list(dict_obj.keys())[0]

    def elem2node_values(self, vmap):
        """
        equation converting 2D values on element to nodes.
        """
        vmap = np.asarray(vmap)
        vnode = [vmap[list(n)].mean(axis=0) for n in self.mesh.node2elems]
        return vnode

    def node2elem_values(self, vmap):
        vmap = np.asarray(vmap)
        velem = [vmap[el].mean(axis=0) for el in self.elems]
        return velem

    def simple_trend(self, ini_meta=None, inpoly=None):
        """
        Assigning a spatially uniformed value or an equation that's dependent on lat, lon
        """
        if ini_meta:
            value = self.get_value(ini_meta)
        else:
            value = self.get_value(self.ini_meta)  # get variable values

        if inpoly is not None:
            n_hgrid = len(inpoly)
        else:
            n_hgrid = self.n_hgrid

        if isinstance(value, (float, int)):
            # initialize the 3D field
            v_ini = np.full((n_hgrid, self.n_vgrid), value)
            return v_ini
        elif isinstance(value, str):
            if ('max' in value) & ('np.maximumm' not in value):
                value = value.replace('max','np.maximum')
            if ('min' in value) & ('np.minimum' not in value):
                value = value.replace('min','np.minimum')
            if inpoly is not None:
                xy = self.hgrid[inpoly]
                if not isinstance(self.vgrid, (float,int)):
                    z = self.vgrid[inpoly,0]
            else:
                xy = self.hgrid
                if not isinstance(self.vgrid, (float,int)):
                    z = self.vgrid[:,0]
            x = xy[:, 0]
            y = xy[:, 1]
            if self.variable_name in ['SED3D_bed', 'SED3D_bedfrac']:
                value = value.split(',')
                v_ini = np.zeros((n_hgrid, self.n_vgrid, self.n_sdim))
                for i in range(len(value)):
                    v_ini[:, :, i] = eval(value[i])
                return v_ini
            else:
                # x and y based function where x and y are lat and lon.
                vmap = eval(value)
                if isinstance(vmap,(int,float)):
                    vmap = np.ones_like(x)*vmap
                if '2D' in self.centering:
                    return vmap
                else:
                    return self.map_to_3D(vmap, self.n_vgrid)
        else:
            raise("input value not recognised as float, int, or str")

    def text_init(self, ini_meta=None, inpoly=None):
        """
        Reading 2D map from either a .prop (on element) or a .ic (on node) file
        If a prop file is read, the element values will need to be converted to node values.
        """
        if ini_meta:
            text_fn = self.get_value(ini_meta)
        else:
            text_fn = self.get_value(self.ini_meta)
        if text_fn.endswith('.ic') or text_fn.endswith('.gr3'):
            if self.centering == 'node2D':
                reader = SchismMeshGr3Reader()  # the ic file has the same format as the .gr3 file
                reader.read(fpath_mesh=text_fn)
                icmesh = reader.read(fpath_mesh=text_fn)

                if icmesh.n_nodes() == self.n_nodes:
                    vmap = icmesh.nodes[:, 2]
                else:
                    raise Exception(
                        "The node of %s is incompatible with the grid node in hgrid.gr3" % text_fn)

                # if the required input is on element, but values on nodes are provided.
                if len(vmap) != self.n_hgrid:
                    vmap = self.node2elem_values(vmap)
            else:
                raise ValueError('%s gives value on node, but value on %s is required' % (
                    text_fn, self.centering))
        elif text_fn.endswith('.prop'):
            if self.centering in ['elem', 'prism']:
                with open(text_fn) as f:
                    content = f.read().splitlines()
                vmap = [float(x.split()[1]) for x in content]
                if len(vmap) != self.n_elems:
                    raise ValueError(
                        "The element of {} is incompatible with grid element in hgrid.gr3".format(text_fn))
                # if the required input is on nodes, but values on elements are provided
                if len(vmap) != self.n_hgrid:
                    vmap = self.elem2node_values(vmap)
            else:
                raise ValueError('%s gives value on element, but value on %s is required' % (
                    text_fn, self.centering))
        elif text_fn.endswith('.nc'):
            ncdata = xr.open_dataset(text_fn)
            vmap = ncdata[self.variable_name].values
            vdim = np.shape(vmap)
            nh = vdim[0]  # horizontal grid dimension of the input nc file.
            if nh != self.n_hgrid:
                if self.n_hgrid == self.mesh.n_nodes:
                    vmap = self.elem2node_values(vmap)
                elif self.n_hgrid == self.mesh.n_elems:
                    vmap = self.node2elem_values(vmap)
                else:
                    raise ValueError(
                        "The horizontal dim of %f does not match with the mesh grid" % text_fn)
        else:
            raise TypeError(
                "Only .ic, .gr3, .prop, or .nc are allowed as text input")

        vmap = np.squeeze(vmap)
        vdim = np.shape(vmap)
        if len(vdim) == 1:  # 1D array
            if self.centering == 'node2D':
                return vmap[inpoly]
            else:
                v = self.map_to_3D(vmap, self.n_vgrid)
                return v[inpoly, :]
        elif len(vdim) == 2:  # 2D array
            if self.vname not in ['SED3D_bed', 'SED3D_bedfrac']:
                return vmap[inpoly]
            else:
                v = self.map_to_3D(vmap, self.n_vgrid)
                v = np.transpose(v, [0, 2, 1])
                return v[inpoly]
        else:  # if alrady 3D array
            return vmap[inpoly]

    def patch_init(self, poly_fn=None):  # regional based initializer
        try:
            poly_fn = self.ini_meta['regions_filename']
        except KeyError:
            if not poly_fn:
                raise EOFError(
                    "regional polygon file is needed for pathch_ini implementation")
        if poly_fn.endswith('shp') or poly_fn.endswith('ic'):
            # perform contiguity check and return a mapping array if successful.
            mapping = geo_tools.partition_check(self.mesh, poly_fn,
                                                self.centering,
                                                crs=self.crs)
        else:
            raise NotImplementedError("Poly_fn can only be shapefile or ic file")

        if self.ini_meta['smoothing']:
            raise NotImplementedError("Smoothing not implemented yet")

        v_merge = np.zeros((self.n_hgrid, self.n_vgrid))
        for i, r in enumerate(self.ini_meta['regions']):
            ini_meta = r['initializer']
            initializer_key = self.get_key(ini_meta)
            initializer = self.patch_initializer_options(initializer_key)
            ini_meta = self.get_value(ini_meta)
            #inpoly = np.where(mapping == i+1)[0]  # mapping is one-based.
            inpoly = np.where(mapping==r['region'])[0]
            v = initializer(ini_meta, inpoly)
            if np.any(np.isnan(v)):
                raise ValueError("region %f has nan in %s field" %
                                 (r, self.variable_name))
            if self.variable_name == 'tke':
                v_merge[:,inpoly,:] = v
            else:
                if type(v) == xr.core.dataarray.DataArray:
                    v_merge[inpoly, :] = np.array(v)[:,np.newaxis] # for cases where hotstart is returning 2D array
                else:
                    try:
                        v_merge[inpoly, :] = v
                    except ValueError:
                        v_merge[inpoly, :] = v[:,np.newaxis]
            #print(i, r)
        return v_merge

    # USGS cast (2D observational points)
    def extrude_casts(self, ini_meta=None, inpoly=None):
        if not ini_meta:
            ini_meta = self.ini_meta
        station_fn = ini_meta['station']
        data_fn = ini_meta['data']
        variable = ini_meta['variable']
        date = pd.to_datetime(self.date)

        polaris_data = pd.read_csv(data_fn, skiprows=[1])
        # find the nearest date for the cruise
        polaris_data['datetime'] = pd.to_datetime(polaris_data.Date)
        polaris_date = pd.to_datetime(polaris_data.Date.unique())
        cast_date = polaris_date[np.argmin(abs(polaris_date-date))]
        if 'Station' in polaris_data.columns:
            polaris_cast = polaris_data[polaris_data.datetime == cast_date][[
                'Station', variable, 'Depth (m)']]
        else:
            polaris_cast = polaris_data[polaris_data.datetime == cast_date][[
                'Station Number', variable, 'Depth']]
            polaris_cast['Station'] = polaris_cast['Station Number'].astype(
                int)
        # setting 'Station' as the index below does not remove depth in the column
        polaris_cast.set_index('Station', inplace=True)
        if variable.lower() in polaris_cast.columns.str.lower():
            # change column to lower case
            polaris_cast.columns = map(str.lower, polaris_cast.columns)
            indnan = np.where(np.isnan(polaris_cast[variable.lower()]))[0]
            if len(indnan) > 0:
                polaris_cast = polaris_cast[~indnan]
        else:
            raise IndexError("The variable %s is not in %s" %
                             (variable, data_fn))

        # find the corresponding station locations
        stations = pd.read_csv(station_fn,sep=",",header=0)[
            ['Station', 'y', 'x']]
        stations.set_index('Station', inplace=True)
        polaris_cast = polaris_cast.join(stations, on='Station', how='left')

        if inpoly is not None:
            n_hgrid = len(inpoly)
        else:
            n_hgrid = self.n_hgrid
        v = np.zeros((n_hgrid, self.n_vgrid))

        # partition the grid domain by polaris stations based on horizontal distance
        grid_df = pd.DataFrame({})
        g_xy = self.hgrid[inpoly]
        #if self.crs == 'EPSG:26910': # utm 10N
        #    g_xy = geo_tools.utm2ll(g_xy.T).T  # conver to (lon, lat).

        for s in polaris_cast.index.unique():
            s_xy = [stations.loc[s].x, stations.loc[s].y]
            grid_df[str(s)] = distance.cdist([s_xy], g_xy)[0]
        grid_df['nearest'] = grid_df.idxmin(axis=1).astype(float).astype(int)

        if 'depth' not in polaris_cast.columns:
            polaris_cast = polaris_cast.rename(columns={'depth (m)':'depth'})

        # loop through each station and perform vertical interpolations for the nearest grid
        # the depths computed are negative but polaris depths are positive
        grid_depths = self.vgrid[inpoly, :]*-1.0
        for s in polaris_cast.index.unique():
            grid_nearest = np.where(grid_df.nearest == s)[0]
            zc = polaris_cast.loc[s]
            f = interpolate.interp1d(
                zc['depth'], zc[variable.lower()], fill_value='extrapolate')
            for g in grid_nearest:
                depth = grid_depths[g]
                v[g, :] = f(depth)

        if np.any(np.isnan(v)):
            raise ValueError(
                "The interpolated %s field has nan in it" % variable)
        return v

#    def extrude_casts(self,ini_meta=None,inpoly=None): # USGS cast (2D observational points)
#        if not ini_meta:
#            ini_meta = self.ini_meta
#        station_fn = ini_meta['station']
#        data_fn = ini_meta['data']
#        variable = ini_meta['variable']
#        date = pd.to_datetime(self.date)
#        polaris_data = pd.read_csv(data_fn,skiprows=[1])
#        # find the nearest date for the cruise
#        polaris_data['datetime'] = pd.to_datetime(polaris_data.Date)
#        polaris_date = pd.to_datetime(polaris_data.Date.unique())
#        cast_date = polaris_date[np.argmin(abs(polaris_date-date))]
#        polaris_cast = polaris_data[polaris_data.datetime == cast_date][['Station Number',variable,'Depth']]
#        polaris_cast['Station'] = polaris_cast['Station Number'].astype(int)
#        polaris_cast.set_index('Station',inplace=True)
#        if variable.lower() in polaris_cast.columns.str.lower():
#            polaris_cast.columns = map(str.lower,polaris_cast.columns) # change column to lower case
#            indnan = np.where(np.isnan(polaris_cast[variable.lower()]))[0]
#            if len(indnan)>0:
#                polaris_cast = polaris_cast[~indnan]
#        else:
#            raise IndexError("The variable %s is not in %s"%(variable,data_fn))
#
#        # find the corresponding station locations
#        stations = pd.read_csv(station_fn)[['Station','North_Lati','West_Longi']]
#        stations.set_index('Station',inplace=True)
#        polaris_cast = polaris_cast.join(stations,on='Station',how='left')
#        polaris_cast['Point'] = [Point(x,y) for x,y in zip(polaris_cast.West_Longi,polaris_cast.North_Lati)]
#        multi_points = MultiPoint(polaris_cast['Point'].values)
#
#        # loop through all cells
#        if inpoly is not None:
#            n_hgrid = len(inpoly)
#        else:
#            n_hgrid = self.n_hgrid
#        v = np.zeros((n_hgrid,self.n_vgrid))
#        for i, xy in enumerate(self.hgrid[inpoly]):
#            x = xy[0]
#            y = xy[1]
#            p = Point((x,y))
#            zc = polaris_cast.iloc[np.where(polaris_cast.Point ==
#                                       nearest_points(p,multi_points)[1])[0]][[variable.lower(),'depth']]
#            depths = self.vgrid[i,:]*-1.0 # the depths computed are negative but polaris depths are positive
#            # interp- and extrapolate over depths.
#            f = interpolate.interp1d(zc.depth,zc[variable.lower()],fill_value='extrapolate',kind='previous')
#            v[i,:] = f(depths)
#            if np.any(np.isnan(v)):
#                raise ValueError("The interpolated %s field has nan in it"%variable)
#        return v

    def obs_points(self, ini_meta=None, inpoly=None):  # 1D observational points
        """
        obs_file includes Lat, Lon, and values for the variable.
        """
        from . import interp_2d  # now onl 2D interpolation IDW is implemented.
        if ini_meta:
            obs_file = ini_meta['data']
            variable = ini_meta['variable']
            if 'ratio' in ini_meta.keys():
                ratio = ini_meta['ratio']
                if isinstance(ratio, str):
                    ratio = eval(ratio)
            else:
                ratio = 1.0
        else:
            obs_file = self.ini_meta['data']
            variable = self.ini_meta['variable']
            if 'ratio' in self.ini_meta.keys():
                ratio = self.ini_meta['ratio']
                if isinstance(ratio, str):
                    ratio = eval(ratio)
            else:
                ratio = 1.0
        obs_data = pd.read_csv(obs_file)
        obs_data = obs_data.dropna(subset=[variable])
        obs_loc = obs_data[['x', 'y']].values
        #print(obs_loc)
        #print(type(obs_loc))
        #obs_loc = geo_tools.ll2utm([obs_loc[:,0],obs_loc[:,1]],self.crs).T
        #print(obs_loc)
        vals = obs_data[variable].values
        #print(vals)

        if inpoly is not None:
            node_xy = self.hgrid[inpoly]
        else:
            node_xy = self.hgrid
        
        tree = interp_2d.Invdisttree(obs_loc)
        tree.weights(node_xy,p=2)
        vmap = tree.interp(vals) 
        
        if np.any(np.isnan(vmap)):
            raise ValueError("vmap has nan value in it.Check %s" % obs_file)
        v = self.map_to_3D(vmap, self.n_vgrid)
        return v

    def hotstart_nc(self, ini_meta=None, inpoly=None):
        var = self.variable_name
        if var=='temperature':
            var='TEM'
        elif var=='salinity':
            var='SAL'
        if not ini_meta:
            ini_meta = self.ini_meta
        data_source = ini_meta['data_source']  #input hotstart file
        yaml_var = {'elevation':'eta2',
                    'velocity_w': 'we',
                    'velocity_u':'su2',
                    'velocity_v':'sv2',
                    'TEM': 'tr_nd',
                    'SAL':'tr_nd',   #['tr_el','tr_nd','tr_nd0']
                    'tke':['q2', 'xl', 'dfv', 'dfh', 'dfq1',
                           'dfq2'],
                    'SED3D_dp':'SED3D_dp',
                    'SED3D_rough':'SED3D_rough',
                    'SED3D_bed':'SED3D_bed',
                    'SED3D_bedfrac':'SED3D_bedfrac'}
        if (var not in yaml_var) and (var in self.tr_mname):
            yamel_var[var] = 'tr_nd'
        if var in self.tr_mname:
            self.tr_index = np.where(np.array(self.tr_mname)==var)[0][0]

        hotstart_data = xr.open_dataset(data_source)

        if 'source_hgrid' not in ini_meta.keys(): #if the grids are exactly the same
            if self.tr_index is not None:
                v = hotstart_data[yaml_var[var]].sel(ntracers=self.tr_index)
            else:
                v = hotstart_data[yaml_var[var]]
        elif ('source_hgrid' in ini_meta.keys()) and ('source_vgrid' not in ini_meta.keys()):
            if self.tr_index is not None:
                vin = hotstart_data[yaml_var[var]].sel(ntracers=self.tr_index)
            else:
                vin = hotstart_data[yaml_var[var]]
            if 'distance_threshold' in ini_meta.keys():
                v = self.interp_from_mesh(ini_meta['source_hgrid'], vin, inpoly,
                                          dist=ini_meta['distance_threshold'],
                                          method=ini_meta['method'])
            else:
                v = self.interp_from_mesh(ini_meta['source_hgrid'], vin, inpoly)
        elif ('source_hgrid' in ini_meta.keys()) and ('source_vgrid' in ini_meta.keys()):
            if self.tr_index is not None:
                vin = hotstart_data[yaml_var[var]].sel(ntracers=self.tr_index)
            else:
                vin = hotstart_data[yaml_var[var]]
            if 'distance_threshold' in ini_meta.keys():
                 v = self.interp_from_mesh(ini_meta['source_hgrid'], vin,
                                          ini_meta['source_vgrid'],
                                          inpoly=inpoly,
                                          dist_th=ini_meta['distance_threshold'],
                                          method=ini_meta['method'],
                                          vgrid_version=ini_meta['vgrid_version'])
            else:
                v = self.interp_from_mesh(ini_meta['source_hgrid'], vin,
                                          ini_meta['source_vgrid'],
                                          ini_meta['source_vgrid_version'],
                                          inpoly)
        # if type(yaml_var[var]) is list:
        #     ds = xr.Dataset()
        #     for j, name in enumerate(yaml_var[var]):
        #         ds_ar = xr.Dataset({name: (['node', 'nVert'],
        #                                    v[j,:,:])})
        #         # coords={'node':range(self.mesh.n_nodes()),
        #         # 'nVert':range(self.mesh.n_vert_levels)})
        #         ds = ds.merge(ds_ar)
        # else:
        #     ds = v
        return v

    def interp_from_mesh(self, hgrid_fn, vin,vgrid_fn=None, vgrid_version=5.8,
                         inpoly=None, dist_th=None, method=None):
        import rtree.index
        mesh1 = read_mesh(hgrid_fn,vgrid_fn,vgrid_version)

        grid1 = self.define_new_grid(mesh1) # mesh to be interpolated from
        hgrid1 = list(grid1.values())[0][1]
        vgrid1 = list(grid1.values())[1][1]

        if inpoly is None:
            hgrid2 = self.hgrid
            vgrid2 = self.vgrid
        else:
            hgrid2 = self.hgrid[inpoly] # mesh to be interpolated into
            vgrid2 = self.vgrid[inpoly,:]

        compare_mesh_flag = False
        if 'hotstart_nc_hfn' in self.hotstart_ini:
            if (self.hotstart_ini['hotstart_nc_hfn'] == hgrid_fn) and \
            ('node' in self.centering or 'prism' in self.centering):
                compare_mesh_flag = True

        if compare_mesh_flag:
            indices = self.hotstart_ini['hotstart_nc_indices']
            dist = self.hotstart_ini['hotstart_nc_dist']
            if inpoly is not None:
                indices = np.array(indices)[inpoly]
                dist = np.array(dist)[inpoly]
        else:
            mesh1_idx = rtree.index.Rtree()
            # Horizontal: nearest neighbor method
            for i, p in enumerate(hgrid1):
                mesh1_idx.insert(i, np.append(p,p))
            dist = []
            indices = []
            print("interpolating horizontal grid from %s"%(
            hgrid_fn))
            for n2 in hgrid2:
                idx = list(mesh1_idx.nearest(tuple(n2)))[0]
                dist2 = (n2[0]-hgrid1[idx,0])**2 + (n2[1]-hgrid1[idx,1])**2
                dist.append(np.sqrt(dist2))
                indices.append(idx)
            print("horizontal interpolation completed!")
        dist = np.array(dist)
        indices = np.array(indices)

        if dist_th is not None: # same points, apply nearest; different points, apply the method defined
            same_points = np.where(dist<=dist_th)[0]
            diff_points = np.where(dist>dist_th)[0]
            if self.mesh.n_vert_levels != mesh1.n_vert_levels:
                raise("distance threshold should only be defined when the two grids have the same vertical layer numbers")
            if isinstance(vin, xr.core.dataset.Dataset): # vin is array for variable from hotstart file
                vout = np.zeros( (len(vin),np.shape(vgrid2)[0],
                                  np.shape(vgrid2)[1]))
                if method == 'nearest':
                    for j,v in enumerate(list(vin.keys())):
                        #vout[j,:] = vin[v][indices]
                        vout[j,same_points] = vin[v][indices[same_points]]
                    for p in diff_points:
                        z1 = vgrid1[indices[p]]
                        z2 = vgrid2[p]
                        for j,v in enumerate(list(vin.keys())):
                            f = interpolate.interp1d(z1, vin[v][indices[p]], 
                                                    fill_value='extrapolate',
                                                    kind='nearest')
                            vout[j, p, :] = f(z2)
                    print("vertical grid interpolation completed!")
                else: # nearest or equation
                    x = hgrid2[diff_points,0]
                    y = hgrid2[diff_points,1]
                    z = vgrid2[diff_points,0]
                    vout[:,diff_points,:] = eval(method)
                    vout[:,same_points,:] = vin[:,indices[same_points],:]
            else:
                vout = np.zeros_like(vgrid2)
                if method == 'nearest':
                    try:
                       vout[same_points] = vin[indices[same_points]]
                    except ValueError:
                        vout[same_points] = vin[indices[same_points]].to_numpy()[:,np.newaxis]
                    for p in diff_points:
                        z1 = vgrid1[indices[p]]
                        z2 = vgrid2[p]
                        f = interpolate.interp1d(z1, vin[indices[p]], 
                                                fill_value='extrapolate',
                                                kind='nearest')
                        vout[p, :] = f(z2)
                    print("vertical grid interpolation completed!")
                else:
                    x = hgrid2[diff_points,0]
                    y = hgrid2[diff_points,1]
                    z = vgrid2[diff_points,0]
                    vout[diff_points,:] = eval(method)[:,np.newaxis]
                    if np.ndim(vin)<np.ndim(vout):
                        vout[same_points,:] = vin[indices[same_points]].values[:,np.newaxis]
                    else:
                        vout[same_points,:] = vin[indices[same_points],:]
        else:

            if self.centering == 'node2D' or \
                (self.mesh.n_vert_levels == mesh1.n_vert_levels and \
                 ('vinterp' not in self.ini_meta or
                  self.ini_meta['vinterp'] is False)):
                if isinstance(vin, xr.core.dataset.Dataset):
                    vout = np.zeros( ((len(vin.keys())),len(indices),
                                      len(vin.nVert)))
                    for j,v in enumerate(list(vin.keys())):
                        vout[j,:] = vin[v][indices]
                else:
                    vout = vin[indices]
            else:
                print("interpolating vertical grid from %s"%(
                    vgrid_fn))
                if isinstance(vin, xr.core.dataset.Dataset):
                    vout = np.zeros( (len(vin),np.shape(vgrid2)[0],
                                      np.shape(vgrid2)[1]))
                else:
                    vout = np.zeros_like(vgrid2)
                for i, n2 in enumerate(hgrid2):
                    z1 = vgrid1[indices[i]]
                    z2 = vgrid2[i]
                    if isinstance(vin, xr.core.dataset.Dataset):
                        for j,v in enumerate(list(vin.keys())):
                            f = interpolate.interp1d(
                            z1, vin[v][indices[i]], fill_value='extrapolate',
                            kind='nearest')
                            vout[j, i, :] = f(z2)
                    else:
                        f = interpolate.interp1d(
                        z1, vin[indices[i]], fill_value='extrapolate',
                        kind='nearest')
                        vout[i, :] = f(z2)
                print("vertical grid interpolation completed!")
        return vout

    def schout_nc(self, var, ini_meta=None, inpoly=None):
        pass

    def create_dataarray(self, var):
        # for prism (tracer variables), tr_el and tr_nd0 need to be calculated.
        if self.centering == 'prism':
            # convert from node to element values for tracers
            # first average to the node center horizontally
            var_temp = self.node2elem_values(var)
            var_temp = np.asarray(var_temp)
            var_el = np.zeros((self.n_elems, self.n_vert_levels))
            # then average to the element center vertically
            var_el[:, 1:] = (var_temp[:, :-1] + var_temp[:, 1:])/2.0
            var_el[:, 0] = var_el[:, 1]  # add one dummy layer.

            name = self.variable_name

            ds = xr.Dataset({'%s_nd' % name: (['node', 'nVert'], var),
                             '%s_nd0' % name: (['node', 'nVert'], var),
                             '%s_el' % name: (['elem', 'nVert'], var_el)})
            # coords={'elem':range(self.n_elems),
            # 'node':range(self.n_nodes),
            # 'nVert':range(self.n_vert_levels)})
        elif self.centering == 'node2D':
            ds_var = xr.DataArray(np.squeeze(var),  # coords=[range(self.n_hgrid)],
                                  dims=[self.hgrid_name], name=self.variable_name)
            ds = ds_var.to_dataset()
        elif self.centering in ['bed', 'bedfrac']:
            ds_var = xr.DataArray(var,  # coords=[range(self.n_hgrid),range(self.n_vgrid),range(self.n_sdim)],
                                  dims=[self.hgrid_name, self.vgrid_name,
                                        self.sdim_name], name=self.variable_name)
            ds = ds_var.to_dataset()
        elif self.variable_name == 'tke':
            var_name = ['q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2']
            ds = xr.Dataset()
            for j, name in enumerate(var_name):
                ds_ar = xr.Dataset({name: (['node', 'nVert'],
                                            var[j,:,:])})
                # coords={'node':range(self.mesh.n_nodes()),
                # 'nVert':range(self.mesh.n_vert_levels)})
                ds = ds.merge(ds_ar)
        else:
            ds_var = xr.DataArray(var,  # coords=[range(self.n_hgrid),range(self.n_vgrid)],
                                  dims=[self.hgrid_name, self.vgrid_name],
                                  name=self.variable_name)
            ds = ds_var.to_dataset()
        return ds


def read_param_nml(nml_fn):
    """
    read param.in file and generate a dict object with key, value pairs for all the parameters
    """
    param = {}
    with open(nml_fn, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            # only take the portion of the arguement before the comment
            line = line.strip().split('!')[0]
            if (line.startswith('#') or line.startswith('!') or
                line.startswith('&') or line.startswith('/') or
                not line):
                continue
            else:
                try:
                    k, v = line.strip().split('=')
                except ValueError:
                    k, v = line.strip().split('==')
                try:  # try converting to numerics.
                    param[k.strip()] = num(v.strip())
                except ValueError:  # if it's only string
                    param[k.strip()] = v.strip()
    return param


def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


def describe_tracers(param_nml, modules=['HYDRO']):
    """ return the number of tracers, the sequence of the tracers and a list
    of tracer names based on the input list of modules and corresponding
    input files.
    Availabel modules are:
!     1: T (default)
!     2: S (default)
!     3: GEN
!     4: AGE
!     5: SED3D: SED
!     6: EcoSim: ECO
!     7: ICM: ICM and/or ICM_PH
!     8: CoSiNE: COSINE
!     9: Feco: FIB
!    10: TIMOR
!    11: FABM
    The script is mainly translated from schism_ini.F90
    """
    ntrs = []  # number of tracers for each module
    tr_mname = []  # tracer moodule name
    irange_tr_1 = []  # the starting indices for each module
    irange_tr_2 = []  # the ending indices for each module

    tr_mname.append('TEM')
    irange_tr_1.append(0)  # start index
    ntrs.append(1)  # temperature tracer
    irange_tr_2.append(1)  # end index

    tr_mname.append('SAL')
    irange_tr_1.append(sum(ntrs))
    ntrs.append(1)  # T,S counted as separate model
    irange_tr_2.append(sum(ntrs))

    if any([m in ['SED', 'AGE', 'GEN', 'ECO'] for m in modules]):
            param = read_param_nml(param_nml)

    if "GEN" in modules and 'ntracer_gen' in param.keys():
        irange_tr_1.append(sum(ntrs))
        ntrs.append(param['ntracer_gen'])
        for t in range(ntrs[-1]):
            tr_mname.append('GEN_%d' % (t+1))
        irange_tr_2.append(sum(ntrs))

    if "AGE" in modules and 'ntracer_age' in param.keys():
        irange_tr_1.append(sum(ntrs))
        ntrs.append(int(param['ntracer_age']/2))
        for t in range(ntrs[-1]):
            tr_mname.append('AGE_%d' % int(t+1))
        irange_tr_2.append(sum(ntrs))

    if "SED" in modules and 'sed_class' in param.keys():
        irange_tr_1.append(sum(ntrs))
        # sediment concentration for different classes (g/m3??)
        ntrs.append(param['sed_class'])
        for t in range(ntrs[-1]):
            tr_mname.append('SED_%d' % (t+1))
        irange_tr_2.append(sum(ntrs))

    if "ECO" in modules and 'eco_class' in param.keys():
        raise NotImplementedError(
            "Initiating ECO from hotstart.nc may not have been implemented")
        irange_tr_1.append(sum(ntrs))
        ntrs.append(param['eco_class'])
        tr_mname.append('ECO')
        irange_tr_2.append(sum(ntrs))

    if "ICM" in modules:
        irange_tr_1.append(sum(ntrs))
        ICM_name = ['ZB1', 'ZB2', 'PB1', 'PB2', 'PB3', 'RPOC', 'LPOC', 'DOC', 'RPON',
                    'LPON', 'DON', 'NH4', 'NO3', 'RPOP', 'LPOP', 'DOP', 'PO4t', 'SU',
                    'SAT', 'COD', 'DOO', 'TIC', 'ALK', 'CA', 'CACO3']
        if "ICM_PH" in modules:
            ntrs.append(25)
            for i in range(25):
                tr_mname.append('ICM_%d' % (i+1))
        else:
            ntrs.append(21)
            for i in range(21):
                tr_mname.append('ICM_%d' % (i+1))
        # The notation for the names can be found in icm.F90.
        irange_tr_2.append(sum(ntrs))

    if "COSINE" in modules:
        irange_tr_1.append(sum(ntrs))
        ntrs.append(13)
        COS_name = ['NO3', 'SiO4', 'NH4', 'S1', 'S2', 'Z1',
                    'Z2', 'DN', 'DSi', 'PO4', 'DOX', 'CO2', 'ALK']
        for t in range(13):
            tr_mname.append('COS_%d' % int(t+1))
        irange_tr_2.append(sum(ntrs))

    if "FIB" in modules:  # Fecal Indicator Bacteria model
        # It appears that FIB has to be read from gr3 files.
        irange_tr_1.append(sum(ntrs))
        ntrs.append(2)
        FIB_name = ['F_coli', 'E_coli']
        for name in FIB_name:
            tr_mname.append('FIB_%' % name)
        irange_tr_2.append(sum(ntrs))

#    if "TIMOR" in modules:
#        tr_mname.append('TMR')

    if "FABM" in modules:
        raise NotImplementedError("Initializing FABM from hotstart.nc is currently not implemented; \
                                  but FABM can be initialized by 'fabm_schism_hotstart.nc. \
                                  see schism_init.F90 for details")
        irange_tr_1.append(sum(ntrs))
        tr_mname.append('FBM')
        irange_tr_2.append(sum(ntrs))

    irange_tr = np.asarray([irange_tr_1, irange_tr_2])

    # Total # of tracers (including T,S)
    ntracers = sum(ntrs)  # including T,S

    return ntracers, ntrs, irange_tr, tr_mname


def hotstart_to_outputnc(hotstart_fn, init_date, hgrid_fn='hgrid.gr3',
                         vgrid_fn='vgrid.in', vgrid_version = 5.8,
                         outname="hotstart_out.nc"):
    """
    convert hotstart.nc to schism output nc file format that can be read by VisIt.
    """
    # load hotstart file and grid file
    hnc = xr.open_dataset(hotstart_fn)

    # rename from hotstart dim name to output dim name
    newname = {'node': 'nSCHISM_hgrid_node',
               'elem': 'nSCHISM_hgrid_face',
               'side': 'nSCHISM_hgrid_edge',
               'nVert': 'nSCHISM_vgrid_layers',
               'idry': 'wetdry_node',
               'idry_e': 'wetdry_elem'}
    hnc = hnc.rename(newname)

    grid_newname = {'n_hgrid_face': 'nSCHISM_hgrid_face',
                    'n_face_node': 'nMaxSCHISM_hgrid_face_nodes',
                    'n_hgrid_edge': 'nSCHISM_hgrid_edge',
                    'n_hgrid_node': 'nSCHISM_hgrid_node',
                    'hgrid_face_nodes': 'SCHISM_hgrid_face_nodes',
                    'hgrid_edge_nodes': 'SCHISM_hgrid_edge_nodes',
                    'hgrid_node_x': 'SCHISM_hgrid_node_x',
                    'hgrid_node_y': 'SCHISM_hgrid_node_y',
                    'hgrid_face_x': 'SCHISM_hgrid_face_x',
                    'hgrid_face_y': 'SCHISM_hgrid_face_y',
                    'hgrid_edge_x': 'SCHISM_hgrid_edge_x',
                    'hgrid_edge_y': 'SCHISM_hgrid_edge_y'}
    # 'n_vgrid_layers':'nSCHISM_vgrid_layers'}

    if not os.path.exists("hgrid.nc"):
        mesh = read_mesh(hgrid_fn, vgrid_fn, vgrid_version)
        write_mesh(mesh, 'hgrid.nc')
    hgrid_nc = xr.open_dataset("hgrid.nc")
    hgrid_nc = hgrid_nc.rename(grid_newname)

    # VisIt only accepts projected coordinates
    for c in ['edge', 'face', 'node']:
        utm_data = np.asarray([hgrid_nc['SCHISM_hgrid_%s_x' % c].values,
                               hgrid_nc['SCHISM_hgrid_%s_y' % c].values])
        gx = xr.DataArray(
            utm_data[0], dims=hgrid_nc['SCHISM_hgrid_%s_x' % c].dims)
        gx.attrs = hgrid_nc['SCHISM_hgrid_%s_x' % c].attrs
        gx.attrs['standard_name'] = 'projection_x_coordinate'
        gx.attrs['units'] = 'm'

        gy = xr.DataArray(
            utm_data[1], dims=hgrid_nc['SCHISM_hgrid_%s_y' % c].dims)
        gy.attrs = hgrid_nc['SCHISM_hgrid_%s_y' % c].attrs
        gy.attrs['standard_name'] = 'projection_y_coordinate'
        gy.attrs['units'] = 'm'

        # update the values
        hgrid_nc['SCHISM_hgrid_%s_x' % c] = gx
        hgrid_nc['SCHISM_hgrid_%s_y' % c] = gy

    # add time dimension to all variable elements
    time = [np.datetime64(init_date)]

    # if 'time' variable already in hnc, drop it.
    if 'time' in list(hnc.variables):
        hnc = hnc.drop('time')

    # add zcor and time dimension
    # zcor is positive upwards: all wet node values are negative.
    zcor = hnc.z.values[np.newaxis, :]
    hgrid_nc['zcor'] = xr.DataArray(
        zcor, dims=['time', 'nSCHISM_hgrid_node', 'nSCHISM_vgrid_layers'])
    hgrid_nc.zcor.attrs['mesh'] = 'SCHISM_hgrid'
    hgrid_nc.zcor.attrs['data_horizontal_center'] = 'node'
    hgrid_nc.zcor.attrs['data_vertical_center'] = 'full'
    hgrid_nc.zcor.attrs['i23d'] = hgrid_nc.z.i23d
    hgrid_nc.zcor.attrs['ivs'] = hgrid_nc.z.ivs
    hgrid_nc.zcor.attrs['missing_value'] = np.nan
    hgrid_nc = hgrid_nc.drop('z')

    # replace tr_nd by the actual variable name.
    for i, v in enumerate(hnc['tracer_list'].values):
        v_values = hnc.tr_nd.isel(ntracers=i).values[np.newaxis, :]
        ds_v = xr.DataArray(v_values,  # coords=[time,n_face,n_vert],
                            dims=['time', 'nSCHISM_hgrid_node',
                                  'nSCHISM_vgrid_layers'],
                            name=v+'_nd')
        ds_v.attrs['mesh'] = "SCHISM_hgrid"
        ds_v.attrs['data_horizontal_center'] = "node"
        ds_v.attrs['data_vertical_center'] = "full"
        ds_v.attrs['i23d'] = 2  # full level on node
        ds_v.attrs['ivs'] = hgrid_nc.zcor.ivs
        hnc = xr.merge([hnc, ds_v])

    # replace tr_el by the actual variable name.
    for i, v in enumerate(hnc['tracer_list'].values):
        v_values = hnc.tr_el.isel(ntracers=i).values[np.newaxis, :]
        ds_v = xr.DataArray(v_values,  # coords=[time,n_face,n_vert],
                            dims=['time', 'nSCHISM_hgrid_face',
                                  'nSCHISM_vgrid_layers'],
                            name=v+'_el')
        ds_v.attrs['mesh'] = "SCHISM_hgrid"
        ds_v.attrs['data_horizontal_center'] = "elem"
        ds_v.attrs['data_vertical_center'] = "half"
        ds_v.attrs['i23d'] = 6  # 3d half level on element
        ds_v.attrs['ivs'] = hgrid_nc.zcor.ivs
        hnc = xr.merge([hnc, ds_v])
    hnc = hnc.drop(['tr_el', 'tracer_list', 'tr_nd', 'tr_nd0'])

    hvar_t = list(hnc.variables)  # entire list
    hvar = [v for v in hvar_t if 'nSCHISM_hgrid' in hnc[v].dims[0]]

#   i23d: indicates location where outputs are defined. 1:3 - node 2D/3D whole/3D half level
#   4:6 - elem 2D/3D whole/half levels; 7:9 - side 2D/3D whole/half levels
#   the i23d value has to be accurate for VisIt to plot.
    for v in hvar:
        dl = len(hnc[v].dims)
        if dl < 3:  # 1D or 2D array, just expand.
            dl = len(hnc[v].squeeze().dims)
            vn = v
            hnc[vn] = hnc[vn].squeeze().expand_dims(
                {'time': time}).astype(float)
            hnc[vn].attrs['mesh'] = "SCHISM_hgrid"
            dhc = hnc[vn].dims[1].split('_')[-1]
            if dhc == 'face':
                dhc = 'elem'
            hnc[vn].attrs['data_horizontal_center'] = dhc

            if dl == 1:
                # 2D are all on full levels.
                hnc[vn].attrs['data_vertical_center'] = "full"
                if hnc[vn].attrs['data_horizontal_center'] == 'node':
                    hnc[vn].attrs['i23d'] = 1
                elif hnc[vn].attrs['data_horizontal_center'] == 'elem':
                    hnc[vn].attrs['i23d'] = 4
                else:  # on edge
                    hnc[vn].attrs['i23d'] = 7
            else:  # I don't think that all variables are on full, but it's unclear how I can make this distinction.
                hnc[vn].attrs['data_vertical_center'] = "half"
                if hnc[vn].attrs['data_horizontal_center'] == 'node':
                    hnc[vn].attrs['i23d'] = 3
                elif hnc[vn].attrs['data_horizontal_center'] == 'elem':
                    hnc[vn].attrs['i23d'] = 6
                else:  # on edge
                    hnc[vn].attrs['i23d'] = 9
            hnc[vn].attrs['ivs'] = hgrid_nc.zcor.ivs
        else:  # more than 1D, the variable needs to be broken down into a few variables with the name extention "_i"
            c3 = hnc[v].dims[2]  # output the name for the last coordinates.
            if c3 == 'MBEDP':
                hnc[vn].attrs['data_vertical_center'] = "full"
                sl = ['thickness', 'age', 'porosity']  # sediment property list
                for i, n in enumerate(sl):
                    vn = "%s_%s" % (v, n)
                    hnc[vn] = hnc[v].isel(MBEDP=i).squeeze(
                    ).expand_dims('time').astype(float)
                    hnc[vn].attrs['mesh'] = "SCHISM_hgrid"
                    dhc = hnc[vn].dims[1].split('_')[-1]
                    if dhc == 'face':
                        dhc = 'elem'
                    hnc[vn].attrs['data_horizontal_center'] = dhc
                    # so far, everything should be on full
                    hnc[vn].attrs['data_vertical_center'] = "full"
                    hnc[vn].attrs['i23d'] = 4
                    hnc[vn].attrs['ivs'] = hgrid_nc.zcor.ivs
            else:
                hnc[vn].attrs['data_vertical_center'] = "full"
                dl = len(hnc[v].squeeze().dims)
                for i in hnc[c3]:
                    vn = "%s_%d" % (v, i+1)
                    hnc[vn] = hnc[v][:, :, i].squeeze().expand_dims(
                        {'time': time}).astype(float)
                    hnc[vn].attrs['mesh'] = "SCHISM_hgrid"
                    dhc = hnc[vn].dims[1].split('_')[-1]
                    if dhc == 'face':
                        dhc = 'elem'
                    hnc[vn].attrs['data_horizontal_center'] = dhc
                    # so far, everything should be on full
                    hnc[vn].attrs['data_vertical_center'] = "full"
                    if dl == 2:  # each variable is then 1D
                        if hnc[vn].attrs['data_horizontal_center'] == 'node':
                            hnc[vn].attrs['i23d'] = 1
                        elif hnc[vn].attrs['data_horizontal_center'] == 'elem':
                            hnc[vn].attrs['i23d'] = 4
                        elif hnc[vn].attrs['data_horizontal_center'] == 'edge':  # on edge
                            hnc[vn].attrs['i23d'] = 7
                        else:
                            raise AttributeError(
                                "The data_horizontal_center option for %s is unavailable" % v)
                    else:  # I don't think that all variables are on full, but it's unclear how I can make this distinction.
                        hnc[vn].attrs['data_vertical_center'] = "half"
                        if hnc[vn].attrs['data_horizontal_center'] == 'node':
                            hnc[vn].attrs['i23d'] = 3
                            hnc[vn].attrs['i23d'] = 6
                        elif hnc[vn].attrs['data_horizontal_center'] == 'edge':  # on edge
                            hnc[vn].attrs['i23d'] = 9
                        else:
                            raise AttributeError(
                                "The data_horizontal_center option for %s is unavailable" % v)

                    hnc[vn].attrs['ivs'] = hgrid_nc.zcor.ivs
    # hnc = hnc.drop(['SED3D_bed', 'SED3D_bedfrac'])
    output_nc = xr.merge([hgrid_nc, hnc])

    # additional variable attributes needed for VisIt to work: give them some dummy values
    hc = xr.DataArray([0.0], dims=['one'], name='sigma_h_c')
    hc.attrs['long_name'] = 'ocean_s_coordinate h_c constant'
    hc.attrs['units'] = 'meters'
    hc.attrs['positive'] = 'down'

    tb = xr.DataArray([0.0], dims=['one'], name='sigma_theta_b')
    tb.attrs['long_name'] = 'ocean_s_coordinate theta_b constant'
    tf = xr.DataArray([0.0], dims=['one'], name='sigma_theta_f')
    tf.attrs['long_name'] = 'ocean_s_coordinate theta_f constant'

    mdepth = xr.DataArray([0.0], dims=['one'], name='sigma_maxdepth')
    mdepth.attrs['long_name'] = 'ocean_s_coordinate maximum depth cutoff (mixed s over z bound...'
    mdepth.attrs['units'] = 'meters'
    mdepth.attrs['positive'] = 'down'

    cs = xr.DataArray(np.nan*np.ones_like(output_nc.nSCHISM_vgrid_layers),
                      # has to be float here.
                      coords=[np.zeros_like(output_nc.nSCHISM_vgrid_layers)],
                      dims=['sigma'], name='Cs')
    cs.attrs['long_name'] = 'Function C(s) at whole levels'
    cs.attrs['postive'] = 'up'
    # the attributes of sigma is very important to make VisIt work for some reason.
    cs['sigma'].attrs['long_name'] = 'S coordinates at whole levels'
    cs['sigma'].attrs['units'] = '1'
    cs['sigma'].attrs['standard_name'] = 'ocean_s_coordinate'
    cs['sigma'].attrs['positive'] = 'up'
    cs['sigma'].attrs['h_s'] = 0.0
    cs['sigma'].attrs['h_c'] = 0.0
    cs['sigma'].attrs['theta_b'] = 0.0
    cs['sigma'].attrs['theta_f'] = 0.0
    cs['sigma'].attrs['formula_terms'] = 's: sigma eta: elev depth: depth a: sigma_theta_f b: sigma...'
    csf = xr.DataArray([2], dims=['one'], name='coordinate_system_flag')
    df = xr.DataArray([0], dims=['one'], name='dry_value_flag')
    df.attrs['values'] = '0: use last-wet value; 1: use junk'
    md = xr.DataArray([0.01], dims=['one'], name='minimum_depth')

    # the following definition of hgrid is not essential for VisIt, but just good to have for information.
    hgrid = xr.DataArray([-2147483647], dims=['one'], name='SCHISM_hgrid')
    hgrid.attrs['long_name'] = 'Topology data of 2d unstructured mesh'
    hgrid.attrs['topology_dimension'] = 2
    hgrid.attrs['cf_role'] = 'mesh_topology'
    hgrid.attrs['node_coordinates'] = 'SCHISM_hgrid_node_x' 'SCHISM_hgrid_node_y'
    hgrid.attrs['face_node_connectivity'] = 'SCHISM_hgrid_face_nodes'
    hgrid.attrs['edge_coordinates'] = 'SCHISM_hgrid_edge_x' 'SCHISM_hgrid_edge_y'
    hgrid.attrs['face_coordinates'] = 'SCHISM_hgrid_face_x' 'SCHISM_hgrid_face_y'
    hgrid.attrs['edge_node_connectivity'] = 'SCHISM_hgrid_edge_nodes'

    output_nc = xr.merge(
        [hgrid, output_nc, hc, tb, tf, mdepth, cs, csf, df, md])

    # remove a few unnecessary variables that make the VisIt crush.
#    output_nc = output_nc.drop(['hgrid_contour_x','hgrid_contour_y'])
#    output_nc = output_nc.drop(['SED_1','SED_2','SED_0','SED3D_bed_thickness',
#                                'SED3D_bed_age','SED3D_bed_porosity','SED3D_bedfrac_1',
#                                'SED3D_bedfrac_2','SED3D_bedfrac_3',
#                                'SED3D_dp','SED3D_rough'])

    # maybe it doesn't like variables on edge?
    output_nc = output_nc.drop(['q2', 'xl', 'su2', 'sv2'])

    output_nc['time'].attrs['long_name'] = 'Time'
    #output_nc['time'].encoding['units'] = 'days since %s'%init_date
    output_nc['time'].attrs['base_date'] = init_date
    output_nc['time'].attrs['standard_name'] = 'time'

    # add global attributes
    output_nc.attrs['Conventions'] = "CF-1.0, UGRID-1.0"
    output_nc.attrs['title'] = "SCHISM Model output"
    output_nc.attrs['institution'] = "SCHISM Model output"
    output_nc.attrs['source'] = "SCHISM model output version v10"
    output_nc.attrs['references'] = "http://ccrm.vims.edu/schismweb/"
    output_nc.attrs['history'] = "converted from hotstart.nc"
    output_nc.attrs['comment'] = "SCHISM Model output"
    output_nc.attrs['type'] = "SCHISM Model output"
    output_nc.attrs['VisIT_plugin'] = "https://schism.water.ca.gov/library/-/document_library/view/3476283"
    output_nc.attrs['_CoordSysBuilder'] = "ucar.nc2.dataset.conv.CF1Convention"

    output_nc.to_netcdf(outname, format='NETCDF4_CLASSIC')

def project_mesh(mesh,new_crs):  #the original mesh has to have attribute crs.
    # lat, lon: 'EPSG:26910'
    # UTM10N: 'EPSG:26910'
    if mesh.crs==new_crs:
        return mesh
    else:
        projection = geo_tools.project_fun(new_crs)
        new_nodes = projection(mesh.nodes[:,0], mesh.nodes[:,1])
        new_nodes = np.asarray(new_nodes).T
        new_nodes = np.append(new_nodes,mesh.nodes[:,2][:,np.newaxis],axis=1)
        mesh.nodes = new_nodes
        return mesh


def create_arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Create hotstart for a schism run")
    parser.add_argument('--yaml_fn',type=str, help='yaml file for hotstart',required=True)
    parser.add_argument('--modules','--list',nargs='+',help='modules activated in schism', default=None)
    parser.add_argument('--crs',type=str,help='The projection system for the mesh',default=None)
    parser.add_argument('--output_fn',type=str,help='Output hotstart filename.nc',default=None)
    return parser

def main():
    # User inputs override the yaml file inputs.
    parser = create_arg_parser()
    args = parser.parse_args()

    h = hotstart(args.yaml_fn,module=args.modules,crs=args.crs,
                 output_fn=args.output_fn)
    if args.output_fn is not None:
        output_fn = args.output_fn
    else:
        output_fn = h.output_fn
    h.create_hotstart()
    hnc = h.nc_dataset
    hnc.to_netcdf(output_fn)

if __name__ == "__main__":
    main()


# if __name__ == '__main__':
#     h = hotstart('hotstart.yaml', modules=['TEM', 'SAL'],crs='EPSG:26910',param_nml="param.nml")
#     #h = hotstart('hotstart_test.yaml')
#     h.read_yaml()
#     h.create_hotstart()
#     v1 = h.nc_dataset

#     # %% making plot
#     coll = h.mesh.plot_elems(v1['tr_el'].values[:, -1, 1])
#     cb = plt.colorbar(coll)
#     plt.axis('off')
#     plt.title('Regional Temperature')
#     plt.tight_layout(pad=1)