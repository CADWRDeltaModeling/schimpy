#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line tool to convert SCHISM polygons in YAML
to Shapefile and vice versa
"""
from schimpy.schism_polygon import read_polygons, write_polygons
import click


@click.command(help="Convert SCHISM polygons between YAML and Shapefile formats.")
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
@click.option(
    "--yaml-root",
    default=None,
    type=str,
    help="Dot-separated path to the parent of the polygons key in a nested YAML file, "
    "e.g. --yaml-root=mesh.depth_enforcement. Only used with YAML input.",
)
@click.help_option("-h", "--help")
def convert_polygons_cli(input, output, yaml_root):
    """CLI wrapper for converting polygon files."""
    convert_polys(input, output, yaml_root=yaml_root)


def convert_polys(input, output, yaml_root=None, envvar=None):
    if input.endswith(".yaml"):
        polygons = read_polygons(input, yaml_root=yaml_root, envvar=envvar)
        if output.endswith(".shp"):
            write_polygons(output, polygons)
        else:
            raise ValueError("Not supported output file type")
    elif input.endswith(".shp"):
        if yaml_root is not None:
            raise ValueError("--yaml-root is only supported with YAML input files.")
        polygons = read_polygons(input)
        write_polygons(output, polygons)
    else:
        raise ValueError(
            "Unsupported input file type. Only YAML (.yaml) or Shapefile (.shp) are supported."
        )


if __name__ == "__main__":
    convert_polygons_cli()
