#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line tool to convert SCHISM line strings in YAML
to Shapefile and vice versa
"""
from schimpy.schism_linestring import read_linestrings, write_linestrings
from schimpy.prepare_schism import get_structures_from_yaml
from schimpy.schism_setup import check_and_suggest
from schimpy.util.yaml_load import yaml_from_file
import fiona
from shapely.geometry import LineString, mapping
import click

schema = {
    "geometry": "LineString",
    "properties": {"name": "str", "type": "str", "configuration": "str"},
}


def convert_structures_main(input_yaml, output_fn):
    """Function for converting line string files from a structures yaml"""
    if input_yaml.endswith(".yaml"):
        inputs = yaml_from_file(input_yaml)
        structures = inputs.get("structures")
        structures = get_structures_from_yaml(structures)
        with fiona.open(
            output_fn, "w", driver="ESRI Shapefile", schema=schema, crs="EPSG:26910"
        ) as shp:
            for struct in structures:
                line = LineString(struct["end_points"])
                shp.write(
                    {
                        "geometry": mapping(line),
                        "properties": {
                            "name": struct["name"],
                            "type": struct["type"],
                            "configuration": str(struct["configuration"]),
                        },
                    }
                )
    else:
        raise ValueError(
            "Unsupported input file type. Only YAML (.yaml) supported for structures."
        )


def convert_linestrings_main(input_fn, output_fn):
    """Function for converting line string files."""
    if input_fn.endswith(".yaml"):
        linestrings = read_linestrings(input_fn)
        if output_fn.endswith(".shp"):
            write_linestrings(output_fn, linestrings)
        else:
            raise ValueError(
                "Unsupported output file type. Only Shapefile (.shp) is supported for YAML input."
            )
    elif input_fn.endswith(".shp"):
        linestrings = read_linestrings(input_fn)
        if output_fn.endswith(".yaml"):
            write_linestrings(output_fn, linestrings)
        else:
            raise ValueError(
                "Unsupported output file type. Only YAML (.yaml) is supported for Shapefile input_fn."
            )
    else:
        raise ValueError(
            "Unsupported input_fn file type. Only YAML (.yaml) or Shapefile (.shp) are supported."
        )


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
@click.option(
    "-s",
    "--structures",
    is_flag=True,
    default=False,
    help="Flag to indicate conversion of structures YAML.",
)
def convert_linestrings_cli(input, output, structures):
    """CLI wrapper for converting line string files."""
    if structures:
        convert_structures_main(input, output)
    else:
        convert_linestrings_main(input, output)


if __name__ == "__main__":
    convert_linestrings_cli()
