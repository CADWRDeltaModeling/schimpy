#!/usr/bin/env python
import math
import click
import schimpy.schism_yaml as schism_yaml

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
    """
    Compute bounding coordinates (upper-left and lower-right) of a GDAL-readable raster.

    Parameters
    ----------
    image : str
        Path to an image/DEM readable by GDAL.

    Returns
    -------
    (xlo, xhi) : tuple[tuple[float, float], tuple[float, float]]
        xlo is (ulx, uly); xhi is (lrx, lry) in the raster’s coordinate system.
    """    
    ds = gdal.Open(image, GA_ReadOnly)
    gt = ds.GetGeoTransform()
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    xlo = (gt[0], gt[3])
    xhi = (xlo[0] + ds.RasterXSize * gt[1], xlo[1] + ds.RasterYSize * gt[5])
    print("Upper left and lower right coords of requested region: %s %s" % (xlo, xhi))
    ds = None
    return xlo, xhi

def load_dem_list(dem_spec_yaml):
    """
    Load a prioritized DEM list from a SCHISM-YAML spec (via `schism_yaml`).

    The YAML must contain either:
      - `dem_list: [...]`
      - or `dem: { dem_list: [...] }`

    Parameters
    ----------
    dem_spec_yaml : str
        Path to YAML file.

    Returns
    -------
    list[str]
        DEM filenames in priority order (high → low).

    Raises
    ------
    FileNotFoundError
        If `dem_spec_yaml` does not exist.
    KeyError
        If no `dem_list` can be found in the supported locations.
    TypeError
        If the `dem_list` is not a list of strings.
    ValueError
        If the list is empty.
    """
    if not os.path.isfile(dem_spec_yaml):
        raise FileNotFoundError(f"DEM spec YAML not found: {dem_spec_yaml}")

    with open(dem_spec_yaml, "r") as f:
        data = schism_yaml.load(f)

    if isinstance(data, dict) and "dem_list" in data:
        dem_list = data["dem_list"]
    elif isinstance(data, dict) and "dem_list" in data["mesh"]:
        dem_list = data["mesh"]["dem_list"]
    elif isinstance(data, dict) and "dem" in data and isinstance(data["dem"], dict) and "dem_list" in data["dem"]:
        dem_list = data["dem"]["dem_list"]
    else:
        raise KeyError(
            "DEM spec YAML must contain either `dem_list: [...]` or `dem: { dem_list: [...] }`"
        )

    if not isinstance(dem_list, list) or not all(isinstance(x, str) for x in dem_list):
        raise TypeError("`dem_list` must be a list of strings")
    if len(dem_list) == 0:
        raise ValueError("`dem_list` is empty")

    return dem_list


def clip_dem(xlo, xhi, demlist="dem.txt", outformat="AAIGrid",
             hshift=False, prefix="dem_clip", verbose=False, merge=False):
    """
    Clip each DEM that intersects the requested bounding box, producing numbered outputs.

    Parameters
    ----------
    xlo, xhi : tuple[float, float]
        Upper-left (xlo) and lower-right (xhi) coordinates for the clipping box.
    demlist : str
        Path to DEM list YAML spec (see `load_dem_list`).
    outformat : str
        GDAL output format name (e.g. 'AAIGrid', 'GTiff', 'PNG', ...).
    hshift : bool
        Deprecated SMS < 11.1 half-cell shift option.
    prefix : str
        Output filename prefix; outputs are `{prefix}_{i}.{ext}`.
    verbose : bool
        If True, print details and show GDAL command.
    merge : bool
        If True, merge all clipped DEMs into a single output respecting priority order.
    """
    filelist = load_dem_list(demlist)

    iout = 0
    clipped_files = []  # Track clipped files for merging
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
            outname = os.path.abspath("%s_%s.%s" % (prefix, iout, extension))
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
            if not merge:
                print("Output file: %s" % outname)
            clipped_files.append(outname)

    if merge and len(clipped_files) > 0:
        merge_output = os.path.abspath("%s_merged.%s" % (prefix, extension))
        print("\nMerging %d clipped DEMs into: %s" % (len(clipped_files), merge_output))
        
        # Determine output grid parameters from first file
        first_ds = gdal.Open(clipped_files[0], GA_ReadOnly)
        first_gt = first_ds.GetGeoTransform()
        base_proj = first_ds.GetProjection()
        dx = first_gt[1]
        dy = first_gt[5]
        first_band = first_ds.GetRasterBand(1)
        nodata = first_band.GetNoDataValue()
        if nodata is None:
            nodata = -9999.0
        dtype = first_band.DataType
        
        # Calculate unified output grid dimensions
        out_cols = int(np.ceil((xhi[0] - xlo[0]) / dx))
        out_rows = int(np.ceil((xlo[1] - xhi[1]) / abs(dy)))
        out_gt = (xlo[0], dx, 0, xlo[1], 0, dy)
        
        # Initialize output array with nodata
        out_array = np.full((out_rows, out_cols), nodata, dtype=np.float32)
        
        # Read all clipped DEMs in reverse order (low to high priority)
        for fname in reversed(clipped_files):
            ds = gdal.Open(fname, GA_ReadOnly)
            if ds is None:
                print(f"Warning: Could not open {fname} for merging")
                continue
            
            # Get this dataset's geotransform and data
            ds_gt = ds.GetGeoTransform()
            band = ds.GetRasterBand(1)
            data = band.ReadAsArray()
            ds_nodata = band.GetNoDataValue()
            if ds_nodata is None:
                ds_nodata = nodata
            
            # Calculate position in output grid
            ds_ulx = ds_gt[0]
            ds_uly = ds_gt[3]
            
            # Column and row offset in output grid
            col_off = int(np.round((ds_ulx - xlo[0]) / dx))
            row_off = int(np.round((xlo[1] - ds_uly) / abs(dy)))
            
            # Extract the region to paste
            ds_rows, ds_cols = data.shape
            
            # Compute valid region bounds
            out_row_start = max(0, row_off)
            out_row_end = min(out_rows, row_off + ds_rows)
            out_col_start = max(0, col_off)
            out_col_end = min(out_cols, col_off + ds_cols)
            
            ds_row_start = max(0, -row_off)
            ds_row_end = ds_row_start + (out_row_end - out_row_start)
            ds_col_start = max(0, -col_off)
            ds_col_end = ds_col_start + (out_col_end - out_col_start)
            
            # Extract the valid data region
            data_region = data[ds_row_start:ds_row_end, ds_col_start:ds_col_end]
            valid_mask = data_region != ds_nodata
            
            # Overwrite output where this dataset has valid data
            out_array[out_row_start:out_row_end, out_col_start:out_col_end][valid_mask] = data_region[valid_mask]
            
            if verbose:
                print(f"  Merged {os.path.basename(fname)}: {np.sum(valid_mask)} valid cells")
            
            ds = None
        
        # Write merged output
        # AAIGrid doesn't support Create(), so create in memory/GTiff first
        driver = gdal.GetDriverByName('GTiff')
        temp_output = os.path.abspath("%s_merged_temp.tif" % prefix)
        out_ds = driver.Create(
            temp_output,
            out_cols,
            out_rows,
            1,
            dtype
        )
        if out_ds is None:
            raise Exception("Failed to create temporary output file")
        out_ds.SetGeoTransform(out_gt)
        out_ds.SetProjection(base_proj)
        out_band = out_ds.GetRasterBand(1)
        out_band.SetNoDataValue(nodata)
        out_band.WriteArray(out_array)
        out_band.FlushCache()
        out_ds = None
        first_ds = None
        
        # Convert to requested format if not GTiff
        if outformat != 'GTiff':
            translate_cmd = "gdal_translate -of %s %s %s" % (outformat, temp_output, merge_output)
            if verbose:
                print("Converting to %s format: %s" % (outformat, translate_cmd))
            p = subprocess.Popen(translate_cmd.split(), shell=True)
            err = p.wait()
            if err:
                raise Exception("Format conversion failed:\n %s" % translate_cmd)
            # Remove temp file
            try:
                os.remove(temp_output)
            except:
                pass
        else:
            # Already in GTiff, just rename
            try:
                if os.path.exists(merge_output):
                    os.remove(merge_output)
                os.rename(temp_output, merge_output)
            except Exception as e:
                print("Warning: Could not rename temp file: %s" % e)
                merge_output = temp_output
        
        print("Merged DEM written to: %s" % merge_output)
        
        # Clean up temporary clipped files when merging
        print("Cleaning up %d temporary files..." % len(clipped_files))
        for fname in clipped_files:
            try:
                os.remove(fname)
                if verbose:
                    print(f"  Removed: {fname}")
            except Exception as e:
                print(f"  Warning: Could not remove {fname}: {e}")


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
@click.option(
    "--merge",
    is_flag=True,
    default=False,
    help="Merge all clipped DEMs into a single output file, respecting priority order (higher priority DEMs overwrite lower priority where they overlap).",
)
@click.help_option("-h", "--help")
def clip_dems_cli(coords, infile, prefix, outformat, verbose, hshift, merge, demlist):
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

    clip_dem(x0, x1, demlist, outformat, hshift, prefix, verbose, merge)


if __name__ == "__main__":
    clip_dems_cli()
