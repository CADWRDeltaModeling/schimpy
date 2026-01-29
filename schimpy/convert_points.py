#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line tool to convert SCHISM points (source and sink) in YAML to Shapefile"""
import click
import geopandas as gpd
import pandas as pd
import warnings
import warnings
from shapely.geometry import Point
from schimpy.schism_sources_sinks import yaml2df
import yaml


def points_to_shp(fpath, points):
    """Convert a DataFrame to a Shapefile using geopandas."""
    gdf = gpd.GeoDataFrame(
        points,
        geometry=[Point(row.x, row.y) for _, row in points.iterrows()],
        crs="EPSG:26910",
    )
    # Only keep relevant columns
    gdf = gdf[["geometry", "sites", "stype"]]
    gdf = gdf.rename(columns={"sites": "site"})
    gdf.to_file(fpath, driver="ESRI Shapefile")

def points_to_yaml(fpath, points):
    # Ensure columns exist
    if "sites" not in points.columns or "stype" not in points.columns:
        raise ValueError("Shapefile must contain 'sites' and 'stype' fields.")
    # Ensure columns exist
    if "sites" not in points.columns or "stype" not in points.columns:
        raise ValueError("Shapefile must contain 'sites' and 'stype' fields.")

    dict_file = {}
    for _, row in points.iterrows():
        dict_file[str(row.sites)] = [round(float(row.x), 2), round(float(row.y), 2)]
    with open(fpath, "w") as file:
        yaml.dump(dict_file, file, default_flow_style=False, allow_unicode=True)
        
def points_to_bp(fpath, points):
    """Convert a DataFrame to a SCHISM build points (.bp) file."""
    import numpy as np
    
    # Calculate cumulative distance along the transect
    points = points.reset_index(drop=True)
    
    # Write to file
    with open(fpath, 'w') as f:
        # Write header
        f.write("point x y z name\n")
        f.write(f"{len(points)}  ! Number of points\n")
        
        # Write each point
        for i, (_, row) in enumerate(points.iterrows(), start=1):
            x = float(row['x'])
            y = float(row['y'])
            z = float(row['z'])
            name = row['sites'] if 'sites' in row else f"point_{i}"
            f.write(f"{i} {x} {y} {z} {name}\n")
            
    print(f"\nBuild points file written to {fpath}.")
    
    
def shp_to_df(fpath):
    """
    Read a point shapefile using geopandas and return a DataFrame
    with 'sites', 'stype', 'x', and 'y' columns.
    """
    gdf = gpd.read_file(fpath)
    
    # Extract coordinates
    gdf["x"] = gdf.geometry.x
    gdf["y"] = gdf.geometry.y
    
    
    # Select relevant columns
    cols = ["x", "y"]
    if "sites" in gdf.columns:
        cols.insert(0, "sites")
    if "stype" in gdf.columns:
        cols.insert(1 if "sites" in gdf.columns else 0, "stype")
    if "z" in gdf.columns:
        cols.insert(1 if "z" in gdf.columns else 0, "z")
    
    if "sites" not in gdf.columns or "stype" not in gdf.columns:
        warnings.warn("Shapefile must contain 'sites' and 'stype' fields for conversion to yaml")
    
    df = gdf[cols].copy()
    return df


def read_points(fpath):
    """Read points from a YAML file."""
    if fpath.endswith(".yaml"):
        return yaml2df(fpath)
    elif fpath.endswith(".shp"):
        return shp_to_df(fpath)
    else:
        raise ValueError("Not supported input file type. Only YAML (.yaml) and Shapefiles (.shp) are supported.")


def write_points(fpath, points):
    """Write points to a Shapefile."""
    if fpath.endswith(".shp"):
        return points_to_shp(fpath, points)
    elif fpath.endswith(".yaml"):
        return points_to_yaml(fpath, points)
    elif fpath.endswith(".bp"):
        return points_to_bp(fpath, points)
    else:
        raise ValueError("Not supported output file type. Only Shapefile (.shp) is supported.")


@click.command(
    help="Convert SCHISM points (source and sink) between YAML and Shapefile formats."
)
@click.option(
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="Input file (.yml, .yaml, .in, or .shp).",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Output file (.yml, .yaml, .in, .bp, or .shp).",
    help="Output file (.yml, .yaml, .in, .bp, or .shp).",
)
@click.help_option("-h", "--help")
def convert_points_cli(input, output):
    """CLI wrapper for converting SCHISM points."""
    main(input, output)


def main(input, output):
    """Function for converting SCHISM points."""
    points = read_points(input)
    write_points(output, points)


if __name__ == "__main__":
    convert_points_cli()

