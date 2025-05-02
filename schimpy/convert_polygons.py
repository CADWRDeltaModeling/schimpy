#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line tool to convert SCHISM polygons in YAML
to Shapefile and vice versa
"""
from schimpy.schism_polygon import read_polygons, write_polygons
import click


@click.command(help="Convert SCHISM polygons between YAML and Shapefile formats.")
@click.help_option("-h", "--help")
@click.option(
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="Input file (YAML or Shapefile).",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Output file (YAML or Shapefile).",
)
@click.help_option("-h", "--help")
def convert_polygons_cli(input, output):
    """CLI wrapper for converting polygon files."""
    convert_polys(input, output)


def convert_polys(input, output):
    if input.endswith(".yaml"):
        polygons = read_polygons(input)
        if output.endswith(".shp"):
            write_polygons(output, polygons)
        else:
            raise ValueError("Not supported output file type")
    elif input.endswith(".shp"):
        polygons = read_polygons(input)
        write_polygons(output, polygons)
    else:
        raise ValueError(
            "Unsupported input file type. Only YAML (.yaml) or Shapefile (.shp) are supported."
        )


if __name__ == "__main__":
    convert_polygons_cli()
