# -*- coding: utf-8 -*-

import os
import yaml
import json
import numpy as np
import pandas as pd
import logging
import geopandas as gpd
from shapely.geometry import Point, Polygon, LineString, mapping, MultiPolygon
from pyproj import Proj, CRS


def shapely_to_geopandas(features, crs=None, shp_fn=None):
    """Convert shapely features to geopandas and generate shapefiles as needed"""
    df = pd.DataFrame()
    df["geometry"] = features
    gdf = gpd.GeoDataFrame(df, geometry="geometry")
    if crs:
        gdf.crs = crs
    else:
        gdf.crs = "EPSG:26910"
    if shp_fn:
        gdf.to_file(shp_fn)
        print("%s generated" % shp_fn)
    return gdf


def ic_to_gpd(fn, crs=None):
    """Read ic yaml file and convert the polygons to geopandas format"""
    with open(fn) as file:
        info = yaml.load(file, Loader=yaml.FullLoader)
    default_value = info["default"]

    regions = []
    methods = []
    polys = []
    pid = []
    for i, p in enumerate(info["polygons"]):
        regions.append(p["name"])
        methods.append(p["attribute"])
        polys.append(Polygon(p["vertices"]))
        pid.append(i)
    df = pd.DataFrame(
        {"id": pid, "region": regions, "method": methods, "geometry": polys}
    )
    gdf = gpd.GeoDataFrame(df, geometry="geometry")
    gdf.crs = crs
    return gdf


def partition_check(
    mesh,
    poly_fn,
    regions,
    centering="node",
    crs=None,
    allow_overlap=False,
    allow_incomplete=False,
):
    """Check if the schism mesh division by the polygon features in poly_fn is unique and complete.
    The partition check is based on either node or element, and the function checks:

    * if there are any orphaned nodes, and
    * if any nodes/elems were assigned to multiple polygons.

    """
    if poly_fn.endswith("shp"):
        poly_gpd = gpd.read_file(poly_fn)
    else:
        poly_gpd = ic_to_gpd(poly_fn, crs)

    if centering == "elem":
        mesh_gpd = mesh.to_geopandas(feature_type="polygon", crs=crs)
        # if the centroid falls in a polygon, the mesh grid is in the polygon
        mesh_gpd["geometry"] = mesh_gpd["geometry"].centroid
    elif centering == "edge":
        mesh_gpd = mesh.to_geopandas(feature_type="edge", crs=crs)
    else:
        mesh_gpd = mesh.to_geopandas(feature_type="point", crs=crs)

    # only perform projection if crs for both poly_gpd and mesh_gpd are provided
    if (poly_gpd.crs is not None) & (mesh_gpd.crs is not None):
        poly_gpd_crs = CRS.from_user_input(poly_gpd.crs)
        mesh_gpd_crs = CRS.from_user_input(mesh_gpd.crs)

        if not poly_gpd_crs.is_exact_same(mesh_gpd_crs):
            # project to the mesh crs.
            poly_gpd = poly_gpd.to_crs(mesh_gpd.crs)

    region_names = np.array([r["region"] for r in regions])

    for i, poly in enumerate(poly_gpd["geometry"]):
        try:
            id_n = np.argwhere(region_names == poly_gpd.region[i])[
                0, 0
            ]  # re-order based on input yaml file
        except IndexError:
            raise Exception(
                "polygon region %s defined in the shapefile but not used by input yaml file."
                % poly_gpd.region[i]
            )
        id_name = "id_%s" % str(id_n)  # make it one-based.
        mesh_gpd[id_name] = mesh_gpd.within(poly)
    # other_id = "id_%s" % str(i+1)

    ID_keys = mesh_gpd.keys()[["id_" in k for k in mesh_gpd.keys()]]
    ID_df = mesh_gpd[ID_keys]
    ID_df = ID_df.reindex(
        sorted(ID_df.columns), axis=1
    )  # sort according to column names

    # check if there is at least one id that each cell belongs to
    orphaned_cells = np.where(~np.any(ID_df, axis=1))[0]

    if len(orphaned_cells) == len(mesh_gpd):
        raise Exception(
            "coordinate system mismatch: the default for the mesh \
                        is utm-xy. Either specify crs for both grid and the \
                        domain polygons, or convert poly_fn to the same \
                        coordinate system as the model grid"
        )
    # elif len(orphaned_cells) >= len(mesh_gpd)/2:
    #     ID_df[other_id] = False
    #     ID_df[other_id].loc[orphaned_cells] = True
    elif len(orphaned_cells) >= 1:
        if allow_incomplete:
            NotImplementedError(
                "categorizing cells based on nearest distance has not been implemented"
            )
        else:
            string = ",".join(orphaned_cells.astype(str))
            raise Exception("Orphaned nodes or cells found at %s" % string)
    # check if there are cells that belong to multiple polygons
    multi_labeled_cells = np.where(np.count_nonzero(ID_df, axis=1) > 1)[0]
    if len(multi_labeled_cells) >= 1:
        if allow_overlap:
            df_multi = ID_df.loc[multi_labeled_cells]
            npoly = len(poly_gpd)
            for i in range(len(df_multi)):
                df = pd.array([False] * npoly)
                id1 = np.where(df_multi.iloc[i])[0][-1]
                df[id1] = True
                df_multi.iloc[i] = df
            ID_df.loc[multi_labeled_cells] = df_multi
        else:
            raise Exception("overalpping domain detected.")
            string = ",".join(multi_labeled_cells.astype(str))
            print(
                """Multiple cells or nodes belong to more than one polygons %s: 
                  categorized as the last True poly"""
            )
            logging.info(
                "The following cells or nodes belong to more than one \
                         region: %s"
                % string
            )
    else:
        # otherwise the domain is contiguous.
        logging.info("The domain division is contiguous!")
    mapping = np.where(ID_df.values)[1] + 1  # change to one-based indices
    # this is only used when the number of orphaned cells exceed half of the original cells
    # regions = np.append(poly_gpd.region.values, 'other')
    mapping = region_names[mapping - 1]
    return mapping


def Polylen(poly):
    # find the number of polygons in a polygon feature
    if isinstance(poly, Polygon):
        poly_len = 1
    elif isinstance(poly, MultiPolygon):
        poly_len = len(poly)
    return poly_len


def FindMultiPoly(poly_array):
    # find the indices for multipolygons in a list or array of polygons
    plens = []
    for poly in poly_array:
        poly_len = Polylen(poly)
        plens.append(poly_len)

    ind = np.where(np.array(plens) > 1)[0]
    return ind


def project_fun(crs=None):
    if crs:
        projection = Proj(crs)
    else:
        # this is utm zone 10.
        projection = Proj("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")
    return projection


def ll2utm(lonlat, crs=None):
    """lonlat can be numpy arrays. lonlat = np.asarray([lon,lat])
    default crs = epsg:26910
    """
    projection = project_fun(crs)
    utm_x, utm_y = projection(lonlat[0], lonlat[1])
    return np.asarray([utm_x, utm_y])


def utm2ll(utm_xy, crs=None):
    """
    utm_xy can be numpy arrays. utm_xy = np.asarray([utm_x,utm_y])
    default crs = epsg:26910
    """
    projection = project_fun(crs)
    lon, lat = projection(utm_xy[0], utm_xy[1], inverse=True)
    return np.asarray([lon, lat])


def geometry2coords(geo_obj):
    coords = list(mapping(geo_obj["geometry"])["coordinates"])
    coords = json.loads(json.dumps(coords))  # nestted tuper to list
    return coords


def geometry2coords_points(geo_obj):
    return list(mapping(geo_obj["geometry"])["coordinates"][0:2])


def shp2yaml(shp_fn, yaml_fn=None, crs=None):
    """Convert a shapefile to yaml file

    Parameters
    ----------
    shp_fn : str
        Input shape filename
    yaml_fn : str
        Output yaml filename
    crs : str, optional
        Output projection. The default is None.

    Returns
    -------
    None.

    """
    gdf = gpd.read_file(shp_fn)
    gdf = gdf.dropna()
    if crs:
        gdf = gdf.to_crs(crs)
    gdf.reset_index(drop=True, inplace=True)
    stype = gdf.geometry.type[0]
    if stype == "LineString":
        stype = "linestrings"
        gdf["coordinates"] = gdf.apply(geometry2coords, axis=1)
        df = gdf.drop(columns="geometry")
        df_yaml = df.to_dict("records")
        df_yaml = {stype: df_yaml}
    elif stype == "Polygon":
        stype = "polygons"
        gdf["vertices"] = gdf.apply(geometry2coords, axis=1)
        df = gdf.drop(columns="geometry")
        df_yaml = df.to_dict("records")
        df_yaml = {stype: df_yaml}
    else:  # if points, no stype needs to be specified in the yaml file.
        stype = shp_fns[sf]
        gdf["coordinates"] = gdf.apply(geometry2coords, axis=1)
        df = gdf.drop(columns="geometry")
        df.set_index("name", inplace=True)
        df_yaml = df.T.to_dict("records")
        df_yaml = {stype: df_yaml}
    if not yaml_fn:
        yaml_fn = "%s.yaml" % os.path.splitext(shp_fn)[0]
    with open(yaml_fn, "w") as file:
        yaml_data = yaml.safe_dump(df_yaml, file)


def yaml2shp(fn, shp_fn=None, crs=None):
    with open(fn, "r") as f:
        yaml_data = yaml.load(f)
    stype = list(
        set(["polygons", "linestrings", "points"]).intersection(yaml_data.keys())
    )[0]
    yaml_df = pd.DataFrame(yaml_data[stype])
    if stype == "linestrings":
        features = [LineString(yc) for yc in yaml_df["coordinates"]]
        yaml_df["geometry"] = features
        gdf = gpd.GeoDataFrame(yaml_df, geometry="geometry")
        gdf = gdf.drop("coordinates", axis=1)
    elif stype == "polygons":
        features = [Polygon(np.squeeze(yc)) for yc in yaml_df["vertices"]]
        yaml_df["geometry"] = features
        gdf = gpd.GeoDataFrame(yaml_df, geometry="geometry")
        gdf = gdf.drop("vertices", axis=1)
    elif stype == "points":  # this only applies dicu
        yaml_df = yaml_df.T
        features = [Point(yc) for yc in yaml_df.values]
        yaml_df["geometry"] = features
        gdf = gpd.GeoDataFrame(yaml_df, geometry="geometry")
        gdf = gdf.drop([0, 1], axis=1)
        gdf["name"] = gdf.index
    if not crs:
        gdf.crs = "EPSG:26910"
    else:
        gdf.crs = crs
    gdf = gdf.rename(columns={"name": "region"})
    if not shp_fn:
        # try to infer shapefile name from the yaml file
        shp_fn = "%s.shp" % os.path.splitext(fn)[0]
    gdf.to_file(shp_fn)
