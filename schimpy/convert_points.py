#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line tool to convert SCHISM points (source and sink) in YAML to Shapefile"""
import click
import pandas as pd
import fiona
from schimpy.schism_sources_sinks import yaml2df


def df_to_shp(fpath, df):
    """Convert a DataFrame to a Shapefile."""
    # Open a fiona object
    pointShp = fiona.open(
        fpath,
        mode="w",
        driver="ESRI Shapefile",
        schema={"geometry": "Point", "properties": [("site", "str"), ("stype", "str")]},
        crs="EPSG:26910",
    )

    for index, row in df.iterrows():
        rowDict = {
            "geometry": {"type": "Point", "coordinates": (row.x, row.y)},
            "properties": {"site": row.sites, "stype": row.stype},
        }
        pointShp.write(rowDict)

    pointShp.close()


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
