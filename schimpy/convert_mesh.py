#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Mesh converter"""
from schimpy.schism_mesh import read_mesh, write_mesh
import click


@click.command(
    help="Convert a mesh from one format to another. The format is decided by the extensions automatically."
)
@click.option(
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="Input mesh file.",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Output mesh file.",
)
@click.option(
    "--crs",
    required=False,
    type=str,
    help="CRS string for the projection.",
)
@click.help_option("-h", "--help")
def convert_mesh_cli(input, output, crs):
    """CLI wrapper for converting mesh files."""
    convert_mesh(input, output, crs)


def convert_mesh(input, output, crs):
    """Converts mesh files to a new format (e.g., shapefile, gr3, etc.)."""
    mesh = read_mesh(input)
    write_mesh(mesh, output, crs=crs)


if __name__ == "__main__":
    convert_mesh_cli()
