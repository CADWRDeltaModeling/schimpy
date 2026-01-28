#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line tool to convert SCHISM points (source and sink) in YAML to Shapefile"""
import click
import geopandas as gpd
import pandas as pd
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

    dict_file = {}
    for _, row in points.iterrows():
        dict_file[str(row.sites)] = [round(float(row.x), 2), round(float(row.y), 2)]
    with open(fpath, "w") as file:
        yaml.dump(dict_file, file, default_flow_style=False, allow_unicode=True)
    
def shp_to_df(fpath):
    """
    Read a point shapefile using geopandas and return a DataFrame
    with 'sites', 'stype', 'x', and 'y' columns.
    """
    gdf = gpd.read_file(fpath)
    # Ensure columns exist
    if "sites" not in gdf.columns or "stype" not in gdf.columns:
        raise ValueError("Shapefile must contain 'sites' and 'stype' fields.")
    # Extract coordinates
    gdf["x"] = gdf.geometry.x
    gdf["y"] = gdf.geometry.y
    # Select relevant columns
    df = gdf[["sites", "stype", "x", "y"]].copy()
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
    else:
        raise ValueError("Not supported output file type. Only Shapefile (.shp) is supported.")


@click.command(
    help="Convert SCHISM points (source and sink) between YAML and Shapefile formats."
)
@click.help_option("-h", "--help")
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
    help="Output file (.yml, .yaml, .in, or .shp).",
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
