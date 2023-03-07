#!/usr/bin/env python
import math
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
    ds = gdal.Open( image, GA_ReadOnly )
    gt = ds.GetGeoTransform()
    cols=ds.RasterXSize
    rows=ds.RasterYSize
    xlo = (gt[0],gt[3])
    xhi = (xlo[0]+ds.RasterXSize*gt[1],xlo[1]+ds.RasterYSize*gt[5])
    print(f"Upper left and lower right coords of requested region: {xlo} {xhi}")
    ds = None
    return xlo,xhi

def clip_dem(xlo,xhi,demlist="dem.txt",outformat="AAIGrid",hshift=False,prefix="clipped",verbose=False):
    filelist = [x.strip() for x in open(demlist,"r").readlines() if x and len(x) > 1 and not x.startswith("#")]
    iout = 0
    if outformat == "AAIGrid":
        extension = 'asc'
    elif outformat == "JPEG":
        extension = "jpg"
    elif outformat == "PNG":
        extension = "png"        
    elif outformat == "GTiff":
        extension = "tif"
    elif outformat == "AIG":
        extension = "adf"
    else:
        raise ValueError("format not supported? (this may be easily fixed if you add the extension you want and its gdal code)")

    for demfile in filelist:
        if demfile.startswith("-"): demfile = demfile[1:].strip()
        print(f"Checking: {demfile}")
        ds = gdal.Open(demfile, GA_ReadOnly )
        gt = ds.GetGeoTransform()
        cols=ds.RasterXSize
        rows=ds.RasterYSize
        ds_xlo = (gt[0],gt[3])
        ds_dx = gt[1]
        ds_dy = gt[5]
        ds_xhi = (ds_xlo[0]+ds.RasterXSize*ds_dx,ds_xlo[1]+ds.RasterYSize*ds_dy)

        if (ds_xhi[0] > xlo[0] and ds_xhi[1] < xlo[1] and \
            ds_xlo[0] < xhi[0] and ds_xlo[1] > xhi[1]):
            complete = ds_xlo[0] < xlo[0] and ds_xlo[1] > xlo[1] and \
                   ds_xhi[0] > xhi[0] and ds_xhi[1] < xhi[1]


            win_xlo = (max(xlo[0], ds_xlo[0]),min(xlo[1], ds_xlo[1]))
            win_xlo = (math.floor((win_xlo[0] - ds_xlo[0])/ds_dx)*ds_dx + gt[0], \
                       math.ceil((win_xlo[1] - ds_xlo[1])/ds_dy)*ds_dy + gt[3])
            win_xhi = (min(xhi[0], ds_xhi[0]),max(xhi[1], ds_xhi[1]))
            win_xhi = (math.ceil((win_xhi[0] - ds_xlo[0])/ds_dx)*ds_dx + gt[0], \
                       math.floor((win_xhi[1] - ds_xlo[1])/ds_dy)*ds_dy + gt[3])

            print("\n********\nDataset %s intersects region and is complete: %s" % (demfile, complete))
            if verbose:
                print(f"Upper left and lower right coords of requested region: {xlo} {xhi}")
                print(f"ds_dx {ds_dx} ds_dy {ds_dy}")
                print("ds origin %s %s" % ds_xlo)
                print(f"Data set upper left {ds_xlo}, lower right: {ds_xhi}")
                print(f"Final upper left {win_xlo}, lower right: {win_xhi}")
                approx_size = (win_xhi[0] - win_xlo[0])*(win_xhi[1] - win_xlo[1])/ds_dx*ds_dy
                print(f"Approx number of raster cells: {int(approx_size)} ")
            outname = f"{prefix}_{iout}.{extension}"
            if hshift:
                shift = f"-a_ullr {win_xlo[0] + ds_dx / 2.0} {win_xlo[1] - ds_dy / 2.0} {win_xhi[0] + ds_dx / 2.0} {win_xhi[1] - ds_dy / 2.0} "
            else:
                shift = ""
            iout += 1
            quiet = "" if verbose else "-quiet"
            command = f"gdal_translate {quiet} -projwin {win_xlo[0]} {win_xlo[1]} {win_xhi[0]} {win_xhi[1]} {shift} -of {outformat} {demfile} {outname}"
            if verbose: print("Calling gdal_translate with command:\n %s" % command)
            p = subprocess.Popen(command.split(), shell=True)
            if err := p.wait():
                raise Exception("Command failed:\n %s" % command)
            print(f"Output file: {outname}")

def create_arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Trim each DEM on a prioritized list. The coordinates used for clipping is supplied either directly as an upper left and lower right coordinate or indirectly using the bounding coordinates of a sample image. In practice this script is usually used with images saved from SMS")
    parser.add_argument('--coords', type=float, nargs=4, metavar=('ul_x','ul_y','lr_x','lr_y'), default = (None, None, None, None),
                        help='bounding coordinates to which DEMs will be clipped (upper left, lower right)')
    parser.add_argument('--image', dest='infile', default = None, help='image or DEM used to infer bounding coordinates for clipping. Use jpeg for the image. This argument is mutually exclusive with --coords. If a sample is provided its upper left and lower right corner will be used.')
    parser.add_argument('--prefix', dest='prefix', default = 'clipped', help='prefix used for output file names')
    parser.add_argument('--outformat', default='AAIGrid',help='output format, default is AAIGrid (ArcInfo ascii.')
    parser.add_argument('--verbose',action = 'store_true', default=False, help='more verbose output.')    
    parser.add_argument('--hshift', action = 'store_true', default=False,help='(deprecated) shift DEM by half cell for applications that incorrectly interpret the location of the origin and data centering of a DEM. This is a bug fix for SMS < 11.1')
    parser.add_argument('demlist', help='file containing prioritized (high to low) list of dems.')
    return parser
    
def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    if (not (args.coords[0] or args.infile)):
        raise ValueError("Either --coords or --sample argument required. See --help for usage.") 
    if (args.coords[0] and args.infile):
        raise ValueError("Arguments --coords and --sample cannot both be supplied. See --help for usage.") 
    
    if args.coords[0]:
        x0 = (args.coords[0], args.coords[1])
        x1 = (args.coords[2], args.coords[3])
    else:
        x0, x1 = bounding_coords(args.infile)
        
    clip_dem(x0,x1,args.demlist,args.outformat,args.hshift,args.prefix, args.verbose)    


        
if __name__ == "__main__":
    main()
