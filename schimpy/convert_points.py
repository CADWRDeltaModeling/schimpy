#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Command line tool to convert SCHISM points (source and sink) in YAML to Shapefile
"""
import argparse
import pandas as pd
import fiona
from schimpy.schism_sources_sinks import yaml2df


def df_to_shp(fpath, df):

    # open a fiona object
    pointShp = fiona.open(
        fpath,
        mode="w",
        driver="ESRI Shapefile",
        schema={"geometry": "Point", "properties": [("site", "str"),("stype", "str")]},
        crs="EPSG:26910",
    )

    for index, row in df.iterrows():
        rowDict = {
            "geometry": {"type": "Point", "coordinates": (row.x, row.y)},
            "properties": {"site": row.sites, "stype": row.stype},
        }
        pointShp.write(rowDict)

    pointShp.close()

    return


def read_points(fpath):
    if fpath.endswith(".yaml"):
        return yaml2df(fpath)
    else:
        raise ValueError("Not supported file type")


def write_points(fpath, df):
    if fpath.endswith(".shp"):
        return df_to_shp(fpath, df)
    else:
        raise ValueError("Not supported file type")


def create_arg_parser():
    """Create an argument parser"""
    parser = argparse.ArgumentParser(
        description="convert SCHISM points (source and sink) in YAML to .shp (GIS)"
    )
    parser.add_argument("--input", help="input file")
    parser.add_argument("--output", help="output file")
    return parser


def main():
    """A main function to convert polygon files"""
    parser = create_arg_parser()
    args = parser.parse_args()
    if args.input.endswith(".yaml"):
        points = read_points(args.input)
        if args.output.endswith(".shp"):
            write_points(args.output, points)
        else:
            raise ValueError("Not supported output file type")
    else:
        raise ValueError("Not supported input file type")


if __name__ == "__main__":
    main()
