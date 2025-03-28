# %%

import pandas as pd
import os
from schimpy.schism_mesh import read_mesh
from schimpy.schism_linestring import read_linestrings, write_linestrings
from schimpy.schism_polygon import read_polygons, write_polygons

from shapely.geometry import Point, LineString
from shapely.geometry.polygon import Polygon
import geopandas as gpd
# %%
# Path to original preprocessing directory
raw_prepro_dir = "../raw_preprocessing/"

# User-specified buffer distance beyond the mesh extent
buffer = 50

# List of input files to process
linestring_yaml = ["open_boundary.yaml", "depth_enforcement_linestrings.yaml", "depth_enforcement_ss_linestrings.yaml"]
polygon_yaml = ["depth_enforcement_polygons.yaml", "depth_enforcement_ss_polygons.yaml"]
shapefiles = ["minmaxlayer_slr_0_mod105.shp"]


# Obtain mesh extent
meshfile = "../mesh_input/test_sdg_mesh.2dm"
mesh = read_mesh(meshfile)
nodes = mesh.nodes

# Use zip to group the elements at each index and find the maximum
max_coordinates = [max(values) for values in zip(*nodes)]
min_coordinates = [min(values) for values in zip(*nodes)]

# Extract domain coverage
domain_polygon = Polygon(
    [
        (min_coordinates[0] - buffer, max_coordinates[1] + buffer),
        (max_coordinates[0] + buffer, max_coordinates[1] + buffer),
        (max_coordinates[0] + buffer, min_coordinates[1] - buffer),
        (min_coordinates[0] - buffer, min_coordinates[1] - buffer),
    ]
)
#

# Filter linestrings
filtered_linestrings = {}
removed_linestrings = {}
for yaml_file in linestring_yaml:
    print("--" * 50)
    print(f"Processing {yaml_file}...")
    # Read the linestrings from the YAML file
    linestring = read_linestrings(os.path.join(raw_prepro_dir, yaml_file))

    # Filter out linestrings that are outside the domain polygon
    filtered_linestrings[yaml_file] = []
    removed_linestrings[yaml_file] = []
    for line in linestring:
        x = line.xy[0].tolist()
        y = line.xy[1].tolist()

        seg = LineString(zip(x, y))
        if domain_polygon.contains(seg):
            print(f"Kept:  {line.prop['name']}")
            filtered_linestrings[yaml_file].append(line)
        else:
            print(f"Removed: {line.prop['name']}")
            removed_linestrings[yaml_file].append(line)

    # Writing filtered linestrings to file
    print(f"Writing filtered linestrings to {yaml_file}...")
    write_linestrings(yaml_file, filtered_linestrings[yaml_file])

# %%

# Filter polygons
# %%
filtered_polygons = {}
removed_polygons = {}
for yaml_file in polygon_yaml:
    print("--" * 50)
    print(f"Processing {yaml_file}...")
    # Read the polygons from the YAML file
    polygon = read_polygons(os.path.join(raw_prepro_dir, yaml_file))

    # Filter out linestrings that are outside the domain polygon
    filtered_polygons[yaml_file] = []
    removed_polygons[yaml_file] = []
    for poly in polygon:

        x = poly.boundary.xy[0].tolist()
        y = poly.boundary.xy[1].tolist()

        seg = LineString(zip(x, y))
        if domain_polygon.contains(seg):
            print(f"Kept:  {poly.prop['name']}")
            filtered_polygons[yaml_file].append(poly)
        else:
            print(f"Removed: {poly.prop['name']}")
            removed_polygons[yaml_file].append(poly)

    # Writing filtered polygons to file
    print(f"Writing filtered polygons to {yaml_file}...")
    write_polygons(yaml_file, filtered_polygons[yaml_file])


# %%

#%%
# Process shape files

# Load the shapefile into a GeoDataFrame
filtered_shapes = {}
removed_shapes = {}
for shapefile in shapefiles:

    filtered_shapes[shapefile] = []
    removed_shapes[shapefile] = []

    gdf = gpd.read_file(os.path.join(raw_prepro_dir, shapefile))

    index = 0
    for shape in gdf.geometry:
        if shape.intersects(domain_polygon):
            print(f"Kept:  index {index}")
            filtered_shapes[shapefile].append(shape)
        else:
            print(f"Removed: index {index}")
            removed_shapes[shapefile].append(shape)
            gdf.drop(index, inplace=True)
        index += 1

    # Save the filtered shapes to a new shapefile
    print(f"Writing filtered shapes to {shapefile}...")
    gdf.to_file(shapefile)
# %%

