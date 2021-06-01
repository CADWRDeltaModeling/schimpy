# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 15:50:34 2019
@author: zzhang
if "no arguments in initialization list" runtime error is seen for pyproj,
it is because that proj_def.dat file is missing.
Solution: find datadir for pyproj.
for example in C:\\Users\\YourUserName\\AppData\\Local\\Continuum\\anaconda3\\Lib\\site-packages\\pyproj\\datadir
Open the file ‘datadir’.
changed ...\\Anaconda3\\share\proj
to ...\\Anaconda3\\Library\\share (where proj_def.dat is).
"""

import yaml
import geopandas as gpd
import pandas as pd
import numpy as np
import logging
from shapely.geometry import Polygon, MultiPolygon
from pyproj import Proj

def shapely_to_geopandas(features,Proj4=None,shp_fn=None):
    """
    convert shapely features to geopandas and generate shapefiles as needed
    """
    df = pd.DataFrame()
    df['geometry'] = features
    gdf = gpd.GeoDataFrame(df,geometry='geometry')
    if Proj4:
        gdf.crs = Proj4
    else:
        gdf.crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_def"
    if shp_fn:
        gdf.to_file(shp_fn)
        print("%s generated"%shp_fn)
    return gdf

def ic_to_gpd(fn,proj4=None):
    """
    Read ic yaml file and convert the polygons to geopandas format
    """
    with open(fn) as file:
        info = yaml.load(file, Loader=yaml.FullLoader)
    default_value = info['default']

    regions = []
    methods = []
    polys = []
    pid = []
    for i, p in enumerate(info['polygons']):        
        regions.append(p['name'])
        methods.append(p['attribute'])    
        polys.append(Polygon(p['vertices']))
        pid.append(i)
    df = pd.DataFrame({'id':pid, 'region':regions,'method':methods,
                       'geometry':polys})           
    gdf = gpd.GeoDataFrame(df,geometry='geometry')
    gdf.crs = proj4
    return gdf

def contiguity_check(mesh,poly_fn,centering='node',proj4=None):
    """
    Check if the schism mesh division by the polygon features in poly_fn is contiguous.
    The contiguity check is based on either node or element, and the function checks
    1) if there are any orphaned nodes, and
    2) if any nodes/elems were assigned to multiple polygons.
    """
    if poly_fn.endswith('shp'):
        poly_gpd = gpd.read_file(poly_fn)
    else:
        poly_gpd = ic_to_gpd(poly_fn,proj4)
        
    if centering == 'elem':
        mesh_gpd = mesh.to_geopandas(feature_type = 'polygon',proj4=proj4)
        # if the centroid falls in a polygon, the mesh grid is in the polygon
        mesh_gpd['geometry'] = mesh_gpd['geometry'].centroid
    elif centering == 'edge':
        mesh_gpd = mesh.to_geopandas(feature_type = 'edge',proj4=proj4)
    else:
        mesh_gpd = mesh.to_geopandas(feature_type = 'point',proj4=proj4)
        
    poly_gpd = poly_gpd.to_crs(mesh_gpd.crs) # project to the mesh crs. 
    poly_gpd.set_index('id',inplace=True)
    poly_gpd.sort_index(inplace=True)
    for id_n, poly in zip(poly_gpd.index,poly_gpd['geometry']):
        id_name = "id_%s"%str(id_n) # make it one-based.
        mesh_gpd[id_name] = mesh_gpd.within(poly)
    other_id = "id_%s"%str(id_n+1)

    ID_keys = mesh_gpd.keys()[ ['id_' in k for k in mesh_gpd.keys()] ]
    ID_df = mesh_gpd[ID_keys]
    # check if there is at least one id that each cell belongs to
    orphaned_cells = np.where(~np.any(ID_df,axis=1))[0]
 
    if len(orphaned_cells) == len(mesh_gpd):
        raise Exception("coordinate system mismatch: the default for the mesh \
                        is lat lon")
    elif len(orphaned_cells) >= len(mesh_gpd)/2:
        ID_df[other_id] = False
        ID_df[other_id].loc[orphaned_cells] = True
    elif len(orphaned_cells) >=1:
        string = ",".join(orphaned_cells.astype(str))
        raise Exception("Orphaned nodes or cells found at %s"%string)
    # check if there are cells that belong to multiple polygons
    multi_labeled_cells = np.where(np.count_nonzero(ID_df,axis=1)>1)[0]
    if len(multi_labeled_cells) >=1:
        df_multi = ID_df.loc[multi_labeled_cells]
        npoly = len(poly_gpd)
        for i in range(len(df_multi)):
            df = pd.array([False]*npoly)
            id1 = np.where(df_multi.iloc[i])[0][-1]
            df[id1] = True
            df_multi.iloc[i] = df
        ID_df.loc[multi_labeled_cells] = df_multi

        string = ",".join(multi_labeled_cells.astype(str))
        print('''These cells or nodes belong to mulitple polygons %s: 
              categorize as the last True poly'''%string)
    # otherwise the domain is contiguous.
    logging.info("The domain divisino is contiguous!")
    mapping = np.where(ID_df.values)[1]+1 # change to one-based indices
    regions = np.append(poly_gpd.region.values,'other') # this is only used when the number of orphaned cells exceed half of the original cells
    mapping = regions[mapping-1]
    return mapping

def Polylen(poly):
    # find the number of polygons in a polygon feature
    if isinstance(poly,Polygon):
        poly_len = 1
    elif isinstance(poly, MultiPolygon):
        poly_len = len(poly)
    return poly_len

def FindMultiPoly(poly_array):
    # find the indices for multipolygons in a list or array of polygons
    plens =[]
    for poly in poly_array:
        poly_len = Polylen(poly)
        plens.append(poly_len)

    ind = np.where(np.array(plens)>1)[0]
    return ind

def project_fun(proj4=None):
    if proj4:
        projection = Proj(proj4)
    else:
        projection = Proj("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs") # this is utm zone 10.
    return projection

def ll2utm(lonlat,proj4=None):
    """
    lonlat can be numpy arrays. lonlat = np.asarray([lon,lat])
    default proj4 = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"
    """
    projection = project_fun(proj4)
    utm_x, utm_y = projection(lonlat[0], lonlat[1])
    return np.asarray([utm_x,utm_y])

def utm2ll(utm_xy,proj4=None):
    """
    utm_xy can be numpy arrays. utm_xy = np.asarray([utm_x,utm_y])
    default proj4 = "+proj=utm +zone=10, +datum=WGS84 +units=m +no_defs"
    """
    projection = project_fun(proj4)
    lon, lat = projection(utm_xy[0],utm_xy[1],inverse=True)
    return np.asarray([lon,lat])