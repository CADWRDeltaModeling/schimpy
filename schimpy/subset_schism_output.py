# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 11:41:00 2021

@author: babban
"""

from osgeo import ogr
from shapely.geometry import shape, Point
import json
import xarray as xr
from math import isnan
import numpy as np
import os
import re
import shutil
from dms_datastore.logging_config import logger
import argparse


def records(file):
    # generator
    reader = ogr.Open(file)
    layer = reader.GetLayer(0)
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        yield json.loads(feature.ExportToJson())

# ------------------------------------------------------------------------------------------------------------
# partition_schout
# ------------------------------------------------------------------------------------------------------------


def partition_schout(in_Schout_Files,
                     partition_shp,
                     combined=True,
                     exclude=[],
                     transform={"salt_depth_ave": ('depth_ave', "salt")}):
    """ Partitions schout binary output files (start with combined) into a set of smaller schout files corresponding to each partition. 

    Initially implementation will be nodal variables but ultimately native centered (edge centered velocity, prism centered tracers). 

    Transform is an opportunity to generate new variables that are well determined with the data from one time step, 
    but not based on time-stenciled operations. 

    Exclude would be an opportunity to prune variables not desired by dropping variables, although certain variables that are required 
    for well-formedness (zcor, wet-dry) and visualization would not be “excludable” so exclude = “all_data” could be used to exclude all 
    but the non-excludable. No output, but side effect is production of files. 
    p0000/schout_1.nc.  """

    # ***********************************************************************
    # 1. Read files

    # 1a. load shapefile
    regions = []
    print('Reading shapefile and first schout file...')
    try:
        mpoly = records(partition_shp)
        for mp in mpoly:
            regions.append(mp)
    except BaseException as err:
        print(f"Error reading partition shape file {err=}, {type(err)=}")

    # 1b. load 1st netcdf file
    try:
        filename = in_Schout_Files[0]
        with xr.open_dataset(filename) as ds:
            data = ds.load()
    except BaseException as err:
        print(f"Error reading schout file list {err=}, {type(err)=}")
        data = None
        raise NotImplementedError(
            "The line above this error will cause an error")
    # 1c. get point list (x and y) values
    x_coord = data['SCHISM_hgrid_face_x']
    y_coord = data['SCHISM_hgrid_face_y']
    hgrid_f = data['nSCHISM_hgrid_face']
    n_hgrid = len(hgrid_f)

    # ***********************************************************************
    # 2. Determine faces with centroids that fall within the shapefile and corresponding
    #   nodes, edges
    #   (should update to use schimpy.mesh directly)

    # 2a. ----Faces---
    facesinROI = []

    print('Determining faces in ROI...')
    for ii in range(n_hgrid):
        point = Point(x_coord[ii], y_coord[ii])
        for mp in regions:
            if point.within(shape(mp['geometry'])):
                facesinROI.append(hgrid_f.item(ii))

    nFacesROI = len(facesinROI)

    # 2b. ---Nodes---
    nodesinROI = []
    hgrid_f_nodes = data['SCHISM_hgrid_face_nodes'].values
    nMax_f_nodes = data['nMaxSCHISM_hgrid_face_nodes'].values

    print('Determining nodes in ROI...')
    for ii in range(nFacesROI):
        facenodes = hgrid_f_nodes[facesinROI[ii]]
        for jj in nMax_f_nodes:
            if (not (isnan(facenodes[jj]))):
                nodesinROI.append(int(facenodes[jj]))

    nodesinROI = list(set(nodesinROI))
    nodesinROI.sort()
    nNodesROI = len(nodesinROI)

    # build mapping list
    node_map_list = {}
    for ii in range(nNodesROI):
        node_map_list[nodesinROI[ii]] = ii+1

    # 2c. ---Edges---
    edgesinROI = []
    hgrid_e_nodes = data['SCHISM_hgrid_edge_nodes'].values
    n_edges = len(hgrid_e_nodes)
    node_np = np.array(nodesinROI)
    nodesinROI = [int(a)-1 for a in nodesinROI]

    print('Determining edges in ROI...')
    for ii in range(n_edges):
        if (hgrid_e_nodes[ii][0] in node_np and hgrid_e_nodes[ii][1] in node_np):
            edgesinROI.append(ii)

    nEdgesROI = len(edgesinROI)

    # ***********************************************************************
    # 3. Build face and edge node indices for subsetted regions.

    # 3a. ---Faces---
    face_to_nodes_local = data['SCHISM_hgrid_face_nodes'].sel(
        nSCHISM_hgrid_face=facesinROI)
    temp_f_nodes_local_vals = face_to_nodes_local.values

    print('Building face node indices for subset...')
    for ii in range(nFacesROI):
        for jj in nMax_f_nodes:
            if (not (isnan(temp_f_nodes_local_vals[ii][jj]))):
                face_to_nodes_local[ii][jj] = node_map_list[temp_f_nodes_local_vals[ii][jj]]

    # 3b. ---Edges---
    edge_to_nodes_local = data['SCHISM_hgrid_edge_nodes'].sel(
        nSCHISM_hgrid_edge=edgesinROI)
    temp_e_nodes_local_vals = edge_to_nodes_local.values

    print('Building edge node indices for subset...')
    for ii in range(nEdgesROI):
        for jj in [0, 1]:
            edge_to_nodes_local[ii][jj] = node_map_list[temp_e_nodes_local_vals[ii][jj]]

    # ***********************************************************************
    # 4. Write out hgrid for subset

    # 4a. Output Folder
    cur_path = os.getcwd()
    subset_path = os.path.join(cur_path, 'subset')
    if (os.path.isdir(subset_path) == False):
        os.mkdir(subset_path)

    filepath = os.path.join(subset_path, 'hgrid.gr3')
    f_grid = open(filepath, "w")

    data_subset = data.sel(nSCHISM_hgrid_face=facesinROI)
    data_subset['SCHISM_hgrid_face_nodes'] = face_to_nodes_local

    data_subset = data_subset.sel(nSCHISM_hgrid_node=nodesinROI)

    data_subset = data_subset.sel(nSCHISM_hgrid_edge=edgesinROI)
    data_subset['SCHISM_hgrid_edge_nodes'] = edge_to_nodes_local

    f_grid.write("Domain_Subset\n")
    f_grid.write("{0} {1} ! # of elements and nodes\n".format(
        nFacesROI, nNodesROI))

    nodex = data_subset['SCHISM_hgrid_node_x'].values
    nodey = data_subset['SCHISM_hgrid_node_y'].values
    elev = data_subset['depth'].values

    for ii in range(nNodesROI):
        f_grid.write("{0:<9} {1:>17.8f} {2:>17.8f} {3:>17.8f}\n".format(
            ii+1, nodex[ii], nodey[ii], elev[ii]))

    f_nodes = data_subset['SCHISM_hgrid_face_nodes'].values

    for ii in range(nFacesROI):
        count = 0

        for jj in nMax_f_nodes:
            if (not (isnan(f_nodes[ii][jj]))):
                count += 1

        f_grid.write("{0} {1}".format(ii+1, count))

        for kk in range(count):
            f_grid.write(" {0}".format(int(f_nodes[ii][kk])))

        f_grid.write("\n")

    f_grid.write("{0} = Number of open boundaries\n".format(0))
    f_grid.write("{0} = Total number of open boundary nodes\n".format(0))

    f_grid.close()

    # if there is a vgrid in the same folder subset vgrid also
    vgrid_file = os.path.join(cur_path, 'vgrid.in')
    vgrid_exist = os.path.isfile(vgrid_file)

    if (vgrid_exist):
        sigma = []
        ivcor = "1"
        nlevels = "1"
        with open(vgrid_file, 'r') as vgrid:
            subset_vgrd_file = os.path.join(subset_path, 'vgrid.in')
            with open(subset_vgrd_file, 'w') as subset_vgrd:
                i = 0
                for line in vgrid:
                    if (i > 1):
                        temp = line.split()
                        sigma.append(temp[1:])
                    elif (i == 0):
                        ivcor = line
                    else:
                        nlevels = line
                    i = i+1
                subset_vgrd.write(ivcor)
                subset_vgrd.write(nlevels)
                for ii in range(nNodesROI):
                    subset_vgrd.write(
                        ("{0:<9}"+" ".join(sigma[ii])+"\n").format(ii+1))

    # ***********************************************************************
    # 5. Select data from each schout file, drop variables, write to file

    # Loop over files and output subset
    for ff in in_Schout_Files:

        # -------------------------------------------------------------------
        # 5. Load file and Select Data
        print("Extracting and saving ROI for {0}...".format(ff))

        if (ff is not in_Schout_Files[0]):
            # del data
            data = None
            with xr.open_dataset(ff) as ds:
                data = ds.load()

        # 5a. ---Faces---
        data_subset = data.sel(nSCHISM_hgrid_face=facesinROI)
        data_subset['SCHISM_hgrid_face_nodes'] = face_to_nodes_local

        # 5b. ---Nodes---
        data_subset = data_subset.sel(nSCHISM_hgrid_node=nodesinROI)

        # 5c. ---Edges---
        data_subset = data_subset.sel(nSCHISM_hgrid_edge=edgesinROI)
        data_subset['SCHISM_hgrid_edge_nodes'] = edge_to_nodes_local

        # 5d. Exclude variables
        data_subset = data_subset.drop_vars(exclude)

        # 5e. Write subset out to file
        filepath = os.path.join(subset_path, ff)
        data_subset.to_netcdf(filepath)
        # data.close()
        del data
        
        
def partition_sribeio(partition_shp,
                     variables=["out2d","zCoordinates","salinity",
                                "horizontalVelX","horizontalVelY"]):
    """ 
        Partitions scribio format output files  into a set of smaller 
        files corresponding to each partition. 
    """

    # ***********************************************************************
    # 1. Read files

    # 1a. load shapefile
    regions = []
    print('Reading shapefile and out2d file...')
    try:
        mpoly = records(partition_shp)
        for mp in mpoly:
            regions.append(mp)
    except BaseException as err:
        print(f"Error reading partition shape file {err=}, {type(err)=}")
    cur_path = os.getcwd()
    out2d_file = None
    out2d_re = re.compile("out2d_[0-9]*.nc")
    
    filenames = next(os.walk(cur_path), (None, None, []))[2] 
    
    for f in filenames:
        if out2d_re.match(f):
            out2d_file = f
            break
    if not(out2d_file):
        logger.error("no out2d file exist, exit")
        print("no out2d file exist, exit")
        exit
        
        
    # 1b. load 1st netcdf file
    try:
        with xr.open_dataset(out2d_file) as ds:
            data_2d = ds.load()
    except BaseException as err:
        print(f"Error reading out2d file list {err=}, {type(err)=}")
        data_2d = None
        raise NotImplementedError(
            "The line above this error will cause an error")
    # 1c. get point list (x and y) values
    x_coord = data_2d['SCHISM_hgrid_node_x']
    y_coord = data_2d['SCHISM_hgrid_node_y']
    hgrid_f = data_2d['nSCHISM_hgrid_face']
    n_hgrid = len(hgrid_f)

    # ***********************************************************************
    # 2. Determine faces with all nodes that fall within the shapefile and corresponding
    #   nodes, edges
    #   

    # 2a. ----Faces---
    # all index stored in faceinROI and nodesinROI starts from 0
    facesinROI = []
    nodesinROI = []
    hgrid_f_nodes = data_2d['SCHISM_hgrid_face_nodes'].values
    nMax_f_nodes = data_2d['nMaxSCHISM_hgrid_face_nodes'].values
    
    print('Determining faces and nodes in ROI...')
    for ii in range(n_hgrid):
        facenodes = hgrid_f_nodes[ii]
        all_nodes_in_region=[True]*nMax_f_nodes
        for jj in nMax_f_nodes:
            not_in_any_region=True
            if (not (isnan(facenodes[jj]))):
                node_id = int(facenodes[jj]-1)
                point = Point(x_coord[node_id], y_coord[node_id])
                #nodesinROI.append(int(facenodes[jj]))
                for mp in regions:
                    if point.within(shape(mp['geometry'])):
                        not_in_any_region=False
                        break
                        
            all_nodes_in_region[jj]=not(not_in_any_region)
        if np.all(all_nodes_in_region):
            facesinROI.append(hgrid_f.item(ii))
            for jj in nMax_f_nodes:
                if (not (isnan(facenodes[jj]))):
                    ## here nodes starts from 1, later will be adjusted from 0 line 354
                    nodesinROI.append(int(facenodes[jj]))
    nFacesROI = len(facesinROI)
    nodesinROI = list(set(nodesinROI))
    nodesinROI.sort()
    nNodesROI = len(nodesinROI)

    # build mapping list
    # this node id map index starts from 1
    node_map_list = {}
    for ii in range(nNodesROI):
        node_map_list[nodesinROI[ii]] = ii+1

    # 2c. ---Edges---
    edgesinROI = []
    hgrid_e_nodes = data_2d['SCHISM_hgrid_edge_nodes'].values
    n_edges = len(hgrid_e_nodes)
    ## node_np starts from 1
    node_np = np.array(nodesinROI)
    nodesinROI = [int(a)-1 for a in nodesinROI]

    print('Determining edges in ROI...')
    for ii in range(n_edges):
        edge_node1 = int(hgrid_e_nodes[ii][0])
        edge_node2 = int(hgrid_e_nodes[ii][1])
        if (edge_node1 in node_np and edge_node2 in node_np):
            edgesinROI.append(ii)

    nEdgesROI = len(edgesinROI)

    # ***********************************************************************
    # 3. Build face and edge node indices for subsetted regions.

    # 3a. ---Faces---
    face_to_nodes_local = data_2d['SCHISM_hgrid_face_nodes'].sel(
        nSCHISM_hgrid_face=facesinROI)
    temp_f_nodes_local_vals = face_to_nodes_local.values

    print('Building face node indices for subset...')
    for ii in range(nFacesROI):
        for jj in nMax_f_nodes:
            if (not (isnan(temp_f_nodes_local_vals[ii][jj]))):
                face_to_nodes_local[ii][jj] =\
                   node_map_list[temp_f_nodes_local_vals[ii][jj]]

    # 3b. ---Edges---
    edge_to_nodes_local = data_2d['SCHISM_hgrid_edge_nodes'].sel(
        nSCHISM_hgrid_edge=edgesinROI)
    temp_e_nodes_local_vals = edge_to_nodes_local.values

    print('Building edge node indices for subset...')
    for ii in range(nEdgesROI):
        for jj in [0, 1]:
            edge_to_nodes_local[ii][jj] =\
                node_map_list[temp_e_nodes_local_vals[ii][jj]]

    # ***********************************************************************
    # 4. Write out hgrid for subset

    # 4a. Output Folder
   
    subset_path = os.path.join(cur_path, 'subset')
    if (os.path.isdir(subset_path) == False):
        os.mkdir(subset_path)

    filepath = os.path.join(subset_path, 'hgrid.gr3')
    f_grid = open(filepath, "w")

    data_subset = data_2d.sel(nSCHISM_hgrid_face=facesinROI)
    data_subset['SCHISM_hgrid_face_nodes'] = face_to_nodes_local

    data_subset = data_subset.sel(nSCHISM_hgrid_node=nodesinROI)

    data_subset = data_subset.sel(nSCHISM_hgrid_edge=edgesinROI)
    data_subset['SCHISM_hgrid_edge_nodes'] = edge_to_nodes_local

    f_grid.write("Domain_Subset\n")
    f_grid.write("{0} {1} ! # of elements and nodes\n".format(
        nFacesROI, nNodesROI))

    nodex = data_subset['SCHISM_hgrid_node_x'].values
    nodey = data_subset['SCHISM_hgrid_node_y'].values
    elev = data_subset['depth'].values

    for ii in range(nNodesROI):
        f_grid.write("{0:<9} {1:>17.8f} {2:>17.8f} {3:>17.8f}\n".format(
            ii+1, nodex[ii], nodey[ii], elev[ii]))

    f_nodes = data_subset['SCHISM_hgrid_face_nodes'].values

    for ii in range(nFacesROI):
        count = 0

        for jj in nMax_f_nodes:
            if (not (isnan(f_nodes[ii][jj]))):
                count += 1

        f_grid.write("{0} {1}".format(ii+1, count))

        for kk in range(count):
            f_grid.write(" {0}".format(int(f_nodes[ii][kk])))

        f_grid.write("\n")

    f_grid.write("{0} = Number of open boundaries\n".format(0))
    f_grid.write("{0} = Total number of open boundary nodes\n".format(0))

    f_grid.close()

    # if there is a vgrid in the same folder subset vgrid also
    vgrid_file = os.path.join(cur_path, 'vgrid.in')
    vgrid_exist = os.path.isfile(vgrid_file)
    subset_vgrid_file = os.path.join(subset_path, 'vgrid.in')
     
    if (vgrid_exist):
        ivcor=1
        with open(vgrid_file,'r') as ff:
            first_line=ff.readline()
            ivcor=int(first_line.split()[0])
            if ivcor == 2:
                shutil.copy(vgrid_file,subset_vgrid_file)
        if ivcor== 1:
            node_sigma = np.loadtxt(vgrid_file, skiprows=3)
            node_bottom = np.loadtxt(vgrid_file, skiprows=2, max_rows=1)
            subset_sigma = node_sigma[:, nodesinROI]
            subset_bottom = node_bottom[:, nodesinROI]
            num_level = subset_sigma.shape[0]
            with open(subset_vgrid_file, 'r') as subset_vgrid:
                subset_vgrid.write("1\n")
                subset_vgrid.write(f"{num_level}\n")
                np.savetxt(subset_bottom, "a")
                np.savetxt(subset_sigma, "a")
            
            
    all_files = next(os.walk(cur_path), (None, None, []))[2]   
    var_res=[re.compile(v+"_([0-9]*).nc") for v in variables]

    # ***********************************************************************
    # 5. Select data from each schout file, drop variables, write to file

    # Loop over files and output subset
    for ff in all_files:

        # -------------------------------------------------------------------
        # 5. Load file and Select Data
         
        is_var = False
        var_id = -1
        for var_re in var_res:
            var_id += 1
            if (var_re.match(ff)):
                is_var = True
                break        
        if (is_var):
            var_name=variables[var_id]
            
            print("Extracting and saving ROI for {0}...".format(ff))
            with xr.open_dataset(ff) as ds:
                    data = ds.load()
            if var_name == "out2d":
                data_subset = data.sel(nSCHISM_hgrid_face=facesinROI)
                data_subset['SCHISM_hgrid_face_nodes'] = face_to_nodes_local
                data_subset = data_subset.sel(nSCHISM_hgrid_node=nodesinROI)
                data_subset = data_subset.sel(nSCHISM_hgrid_edge=edgesinROI)
                data_subset['SCHISM_hgrid_edge_nodes'] = edge_to_nodes_local
            else:  
                data_center = data.variables[var_name].attrs["location"]
                if data_center == "node":
                    data_subset = data.sel(nSCHISM_hgrid_node=nodesinROI)
                elif data_center == "face":
                    data_subset = data.sel(nSCHISM_hgrid_face=facesinROI)
                else:
                    data_subset = data.sel(nSCHISM_hgrid_edge=edgesinROI)

        
            # 5e. Write subset out to file
            filepath = os.path.join(subset_path, ff)
            data_subset.to_netcdf(filepath)
            # data.close()
            del data

#------------------------------------------------------------------------------------------------------------
# subset_schism_output for schout format
#------------------------------------------------------------------------------------------------------------
def subset_schism_output():
    start_file = 2
    end_file = 3
    step = 1
    partition_shp = 'ROI_Suisun.shp'
    in_Schout_Files = []
    exclude =['zcor','Cs']
    
    #files to partition
    for ii in range(start_file,end_file+1,step):
        in_Schout_Files.append("schout_{0}.nc".format(ii))
            
    partition_schout(in_Schout_Files, partition_shp,True,exclude)

#------------------------------------------------------------------------------------------------------------
# subset_schism_output for sribeio format
#------------------------------------------------------------------------------------------------------------
def subset_schism_scribeIO_output(partition_shp):

    partition_sribeio(partition_shp)

# subset_schism_output()


def create_arg_parser():
    """ Create an argument parser
        return: argparse.ArgumentParser
    """

    # Read in the input file
    parser = argparse.ArgumentParser(
        description="""
        
        Extract SCHISM output variables contained in predefined polygons
        and save to a subset folder under current folder
        
        Usage:
        download_hrrr 01/01/2023  g:\temp 15  37.36  39.0 -123.023 -121.16
                     """)
    parser.add_argument('--subset_polygon', default=None, required=True,
                        help='a shape file defines polygons where outputs\
                        are desired')


    return parser

if __name__ == "__main__":
    
    parser = create_arg_parser()
    args = parser.parse_args()
    partition_shp = args.subset_polygon
    subset_schism_scribeIO_output(partition_shp)
# ------------------------------------------------------------------------------------------------------------
# create_partition
# ------------------------------------------------------------------------------------------------------------
def create_partition(mesh, polygons, enforce_exact=False):
    """ Takes a mesh and list of polygons and partitions according to where the centroids (mesh.centroids in schimpy) fall. 
    Parameters
    ------------
    mesh : schimpy.SchismMesh The mesh to partition

    polygons :  Polygons (parsed from shape file or yaml into SCHIMPY

    enforce_exact: bool
    Requires a unique, complete partition. Initially not implemented for True

    Produces meshes and dicitonary(?) maps of local_to_global and global_to_local. (tuple or class)"""


# ------------------------------------------------------------------------------------------------------------
# write_global_local_maps
# ------------------------------------------------------------------------------------------------------------
def write_global_local_maps(dest, global_local, local_global):
    """writes maps out in the schism + visit format. Note that the global owner 
    of a shared node is the lowest rank that contains it."""
