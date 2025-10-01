#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line tool to convert SCHISM points (source and sink) in YAML to Shapefile"""
import click
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point


def df_to_shp(fpath, df):
    """Convert a DataFrame to a Shapefile using geopandas."""
    gdf = gpd.GeoDataFrame(
        df,
        geometry=[Point(row.x, row.y) for _, row in df.iterrows()],
        crs="EPSG:26910",
    )
    # Only keep relevant columns
    gdf = gdf[["geometry", "sites", "stype"]]
    gdf = gdf.rename(columns={"sites": "site"})
    gdf.to_file(fpath, driver="ESRI Shapefile")


def read_points(fpath):
    """Read points from a YAML file."""
    if fpath.endswith(".yaml"):
        return yaml2df(fpath)
    else:
        raise ValueError("Not supported file type")


def write_points(fpath, df):
    """Write points to a Shapefile."""
    if fpath.endswith(".shp"):
        return df_to_shp(fpath, df)
    else:
        raise ValueError("Not supported file type")


@click.command(
    help="Convert SCHISM points (source and sink) between YAML and Shapefile formats."
)
@click.help_option("-h", "--help")
@click.option(
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="Input file (YAML).",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Output file (Shapefile).",
)
@click.help_option("-h", "--help")
def convert_points_cli(input, output):
    """CLI wrapper for converting SCHISM points."""
    main(input, output)


def main(input, output):
    """Function for converting SCHISM points."""
    if input.endswith(".yaml"):
        points = read_points(input)
        if output.endswith(".shp"):
            write_points(output, points)
        else:
            raise ValueError(
                "Unsupported output file type. Only Shapefile (.shp) is supported for YAML input."
            )
    else:
        raise ValueError("Unsupported input file type. Only YAML (.yaml) is supported.")


if __name__ == "__main__":
    convert_points_cli()
