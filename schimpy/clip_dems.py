#!/usr/bin/env python
import math
import click

try:
    from osgeo import gdal
    from osgeo.gdalconst import *

    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal
    from gdalconst import *

import subprocess
import numpy as np
import sys
import os
import gc

DEFAULT_NA_FILL = 6.0


def bounding_coords(image):
    ds = gdal.Open(image, GA_ReadOnly)
    gt = ds.GetGeoTransform()
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    xlo = (gt[0], gt[3])
    xhi = (xlo[0] + ds.RasterXSize * gt[1], xlo[1] + ds.RasterYSize * gt[5])
    print("Upper left and lower right coords of requested region: %s %s" % (xlo, xhi))
    ds = None
    return xlo, xhi


def clip_dem(
    xlo,
    xhi,
    demlist="dem.txt",
    outformat="AAIGrid",
    hshift=False,
    prefix="clipped",
    verbose=False,
):
    filelist = [
        x.strip()
        for x in open(demlist, "r").readlines()
        if x and len(x) > 1 and not x.startswith("#")
    ]
    iout = 0
    if outformat == "AAIGrid":
        extension = "asc"
    elif outformat == "JPEG":
        extension = "jpg"
    elif outformat == "PNG":
        extension = "png"
    elif outformat == "GTiff":
        extension = "tif"
    elif outformat == "AIG":
        extension = "adf"
    else:
        raise ValueError(
            "format not supported? (this may be easily fixed if you add the extension you want and its gdal code)"
        )

    for demfile in filelist:
        if demfile.startswith("-"):
            demfile = demfile[1:].strip()
        print("Checking: %s" % demfile)
        ds = gdal.Open(demfile, GA_ReadOnly)
        gt = ds.GetGeoTransform()
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        ds_xlo = (gt[0], gt[3])
        ds_dx = gt[1]
        ds_dy = gt[5]
        ds_xhi = (
            ds_xlo[0] + ds.RasterXSize * ds_dx,
            ds_xlo[1] + ds.RasterYSize * ds_dy,
        )

        if (
            ds_xhi[0] > xlo[0]
            and ds_xhi[1] < xlo[1]
            and ds_xlo[0] < xhi[0]
            and ds_xlo[1] > xhi[1]
        ):
            complete = (
                ds_xlo[0] < xlo[0]
                and ds_xlo[1] > xlo[1]
                and ds_xhi[0] > xhi[0]
                and ds_xhi[1] < xhi[1]
            )

            win_xlo = (max(xlo[0], ds_xlo[0]), min(xlo[1], ds_xlo[1]))
            win_xlo = (
                math.floor((win_xlo[0] - ds_xlo[0]) / ds_dx) * ds_dx + gt[0],
                math.ceil((win_xlo[1] - ds_xlo[1]) / ds_dy) * ds_dy + gt[3],
            )
            win_xhi = (min(xhi[0], ds_xhi[0]), max(xhi[1], ds_xhi[1]))
            win_xhi = (
                math.ceil((win_xhi[0] - ds_xlo[0]) / ds_dx) * ds_dx + gt[0],
                math.floor((win_xhi[1] - ds_xlo[1]) / ds_dy) * ds_dy + gt[3],
            )

            print(
                "\n********\nDataset %s intersects region and is complete: %s"
                % (demfile, complete)
            )
            if verbose:
                print(
                    "Upper left and lower right coords of requested region: %s %s"
                    % (xlo, xhi)
                )
                print("ds_dx %s ds_dy %s" % (ds_dx, ds_dy))
                print("ds origin %s %s" % ds_xlo)
                print("Data set upper left %s, lower right: %s" % (ds_xlo, ds_xhi))
                print("Final upper left %s, lower right: %s" % (win_xlo, win_xhi))
                approx_size = (
                    (win_xhi[0] - win_xlo[0])
                    * (win_xhi[1] - win_xlo[1])
                    / ds_dx
                    * ds_dy
                )
                print("Approx number of raster cells: %s " % int(approx_size))
            outname = "%s_%s.%s" % (prefix, iout, extension)
            if hshift:
                shift = "-a_ullr %s %s %s %s " % (
                    win_xlo[0] + ds_dx / 2.0,
                    win_xlo[1] - ds_dy / 2.0,
                    win_xhi[0] + ds_dx / 2.0,
                    win_xhi[1] - ds_dy / 2.0,
                )
            else:
                shift = ""
            iout += 1
            quiet = "" if verbose else "--quiet"
            command = "gdal_translate %s -projwin %s %s %s %s %s -of %s %s %s" % (
                quiet,
                win_xlo[0],
                win_xlo[1],
                win_xhi[0],
                win_xhi[1],
                shift,
                outformat,
                demfile,
                outname,
            )
            if verbose:
                print("Calling gdal_translate with command:\n %s" % command)
            p = subprocess.Popen(command.split(), shell=True)
            err = p.wait()
            if err:
                raise Exception("Command failed:\n %s" % command)
            print("Output file: %s" % outname)


@click.command(
    help="Trim each DEM on a prioritized list. The coordinates used for clipping is supplied either directly as an upper left and lower right coordinate or indirectly using the bounding coordinates of a sample image. In practice this script is usually used with images saved from SMS.\n\n"
    "Arguments:\n"
    " DEMLIST file containing prioritized (high to low) list of dems."
)
@click.argument(
    "demlist",
    type=click.Path(exists=True),
)
@click.option(
    "--coords",
    type=(float, float, float, float),
    default=None,
    help="Bounding coordinates to which DEMs will be clipped (upper left x, y, lower right x, y).",
)
@click.option(
    "--image",
    "infile",
    type=click.Path(exists=True),
    default=None,
    help="Image or DEM used to infer bounding coordinates for clipping. Use jpeg for the image. This argument is mutually exclusive with --coords. If a sample is provided its upper left and lower right corner will be used.",
)
@click.option(
    "--prefix",
    default="clipped",
    help="Prefix used for output file names.",
)
@click.option(
    "--outformat",
    default="AAIGrid",
    help="Output format, default is AAIGrid (ArcInfo ASCII).",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Enable verbose output.",
)
@click.option(
    "--hshift",
    is_flag=True,
    default=False,
    help="(Deprecated) Shift DEM by half cell for applications that incorrectly interpret the location of the origin and data centering of a DEM. This is a bug fix for SMS < 11.1",
)
@click.help_option("-h", "--help")
def clip_dems_cli(coords, infile, prefix, outformat, verbose, hshift, demlist):
    """
    Trim each DEM on a prioritized list. The coordinates used for clipping are supplied either directly as an upper left and lower right coordinate or indirectly using the bounding coordinates of a sample image. In practice, this script is usually used with images saved from SMS.
    """
    if not (coords or infile):
        raise click.UsageError(
            "Either --coords or --image argument is required. See --help for usage."
        )
    if coords and infile:
        raise click.UsageError(
            "Arguments --coords and --image cannot both be supplied. See --help for usage."
        )

    if coords:
        x0 = (coords[0], coords[1])
        x1 = (coords[2], coords[3])
    else:
        x0, x1 = bounding_coords(infile)

    clip_dem(x0, x1, demlist, outformat, hshift, prefix, verbose)


if __name__ == "__main__":
    clip_dems_cli()
