#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line tool to convert SCHISM line strings in YAML
to Shapefile and vice versa
"""
from schimpy.schism_linestring import read_linestrings, write_linestrings
import click


@click.command(help="Convert SCHISM line strings between YAML and Shapefile formats.")
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
def convert_linestrings_cli(input, output):
    """CLI wrapper for converting line string files."""
    main(input, output)


def main(input, output):
    """Function for converting line string files."""
    if input.endswith(".yaml"):
        linestrings = read_linestrings(input)
        if output.endswith(".shp"):
            write_linestrings(output, linestrings)
        else:
            raise ValueError(
                "Unsupported output file type. Only Shapefile (.shp) is supported for YAML input."
            )
    elif input.endswith(".shp"):
        linestrings = read_linestrings(input)
        if output.endswith(".yaml"):
            write_linestrings(output, linestrings)
        else:
            raise ValueError(
                "Unsupported output file type. Only YAML (.yaml) is supported for Shapefile input."
            )
    else:
        raise ValueError(
            "Unsupported input file type. Only YAML (.yaml) or Shapefile (.shp) are supported."
        )


if __name__ == "__main__":
    convert_linestrings_cli()
