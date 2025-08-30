"""Routines to fill elevation (or other scalar) values at points from a prioritized list of rasters"""

try:
    from osgeo import gdal
    gdal.UseExceptions()
    from osgeo.gdalconst import *

    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal
    from gdalconst import *

import yaml
import hashlib
from xml.etree.ElementInclude import include
from schimpy.schism_setup import ensure_outdir
from typing import Callable, Dict, List, Optional, Tuple, Union
from datetime import datetime
import numpy as np
import sys
import os
import diskcache as dc
import json



DEFAULT_NA_FILL = 2.0
_MAX_CACHE_KEYS = 5  # hard cap on total keys in this cache namespace

def _stacked_dem_fill(
    files, points, out_dir, values=None, negate=False, require_all=True, na_fill=None
):
    """
    Implementational core for stacked_dem_fill. See that function for details
    """
    if values is None:
        values = np.full(points.shape[0], np.nan)
    else:
        values = np.array(values, copy=True)
    points = points[:, :2]
    if require_all and na_fill:
        raise ValueError("Incompatible: require_all and na_fill")
    for infile in files:
        if not os.path.exists(infile):
            raise ValueError("File does not exist: %s" % infile)
        print("Using DEM: %s" % infile)
        indataset = gdal.Open(infile, GA_ReadOnly)
        gt = indataset.GetGeoTransform()
        rb = indataset.GetRasterBand(1)
        dem = rb.ReadAsArray()
        nd = rb.GetNoDataValue()
        if not (nd is None):
            dem[dem == nd] = np.nan
        filled = bilinear(points[np.isnan(values)], gt, dem)
        filled[filled < -20000] = np.nan
        values[np.isnan(values)] = filled
        indataset = None
        dem = None
        rb = None

    if negate:
        values = np.negative(values)
    if any(np.isnan(values)):
        f_out = ensure_outdir(out_dir, "dem_misses.txt")
        miss_idx = np.isnan(values)
        miss_points = points[miss_idx]
        f = open(f_out, "w")
        for p in miss_points:
            f.write("%s,%s\n" % tuple(p))
            message = "DEMs provided do not cover all the points. See file dem_misses.txt for locations"
        f.close()

        gdf = gpd.GeoDataFrame(
            {"x": miss_points[:, 0], "y": miss_points[:, 1]},
            geometry=[Point(xy) for xy in miss_points],
            crs="EPSG:4326",  # Change if your coordinates are not WGS84
        )
        gdf.to_file(os.path.join(out_dir, "dem_misses.shp"))
        if require_all:
            raise ValueError(message)
        else:
            print(message)
    if na_fill is not None:
        values[miss_idx] = na_fill
    return values



def _cache_dir_for(out_dir: str) -> str:
    # Keep cache colocated with outputs unless you decide otherwise
    cdir = os.path.join(out_dir, ".dem_cache")
    os.makedirs(cdir, exist_ok=True)
    return cdir

def _make_cache_key(
    files, points, values, require_all, na_fill
) -> str:
    """
    Build a stable key that:
      - Includes the DEM file list (absolute paths, ordered)
      - Includes the *content* of points
      - Includes values mask & non-NaN values (to honor prefilled inputs)
      - Includes options that affect the numeric result (require_all, na_fill)
    NOTE: 'negate' is intentionally excluded per requirement.
    """
    files_abs = tuple(os.path.abspath(f) for f in files)
    pts = np.asarray(points, dtype=np.float64)
    if values is None:
        vals = np.full(pts.shape[0], np.nan, dtype=np.float64)
    else:
        vals = np.asarray(values, dtype=np.float64)

    mask = np.isnan(vals)
    packed = {
        "files": files_abs,
        "require_all": bool(require_all),
        "na_fill": None if na_fill is None else float(na_fill),
        # Hash large arrays by content to avoid JSON bloat
        "points_hash": hashlib.blake2b(pts.tobytes(), digest_size=16).hexdigest(),
        "mask_hash": hashlib.blake2b(mask.tobytes(), digest_size=16).hexdigest(),
        "vals_non_nan_hash": hashlib.blake2b(vals[~mask].tobytes(), digest_size=16).hexdigest(),
        "n_points": int(pts.shape[0]),
    }
    key_src = json.dumps(packed, sort_keys=True).encode("utf-8")
    return hashlib.blake2b(key_src, digest_size=16).hexdigest()

def _touch_lru(cache: dc.Cache, key: str):
    """
    Maintain a simple LRU list inside the cache under a reserved meta-key.
    Ensures at most _MAX_CACHE_KEYS items exist in the cache namespace.
    """
    META = "__order__"
    order = cache.get(META, [])
    if key in order:
        order.remove(key)
    order.append(key)
    # Cull least-recently-used beyond cap
    while len(order) > _MAX_CACHE_KEYS:
        old = order.pop(0)
        try:
            del cache[old]
        except KeyError:
            pass
    cache[META] = order

def stacked_dem_fill(
    files, points, out_dir, values=None, negate=False, require_all=True, na_fill=None
):
    """
    Fill values at an array of points using bilinear interpolation from a prioritized stack of DEMs.
    This function manages prioritization of the DEMs, gradually filling points
    that are still marked as `nan`. Uses a cache to speed up repeated queries.

    Parameters
    ----------
    files : list of str
        List of DEM file paths. Existence is checked; raises ValueError if any file is missing.
    points : np.ndarray
        Numpy array of shape (n_points, 2) containing (x, y) coordinates for interpolation.
    out_dir : str or Path
        Directory to write output files and cache.
    values : np.ndarray or None, optional
        Optional array of initial values at each point. If provided, must be same length as points.
        Points with NaN values will be filled.
    negate : bool, optional
        If True, output values will be negated (e.g., convert elevation to depth).
    require_all : bool, optional
        If True, raises ValueError if the DEMs do not cover all points. In either case,
        a file `dem_misses.txt` is created or overwritten with locations of uncovered points.
    na_fill : float or None, optional
        Value to substitute for points that remain NaN after all DEMs are processed.
        If `require_all` is True, `na_fill` must be None.

    Returns
    -------
    np.ndarray
        Array of values at the input points, filled from the DEMs.

    Raises
    ------
    ValueError
        If any DEM file does not exist, or if `require_all` is True and not all points are covered.
    """
    cache = dc.Cache(_cache_dir_for(out_dir))
    key = _make_cache_key(files, points, values, require_all, na_fill)

    if key in cache:
        result = np.asarray(cache[key], dtype=np.float64)
        _touch_lru(cache, key)
    else:
        # IMPORTANT: negate=False here so cache stores the canonical (non-negated) result
        result = _stacked_dem_fill(
            files=files,
            points=np.asarray(points, dtype=np.float64),
            out_dir=out_dir,
            values=values,
            negate=False,
            require_all=require_all,
            na_fill=na_fill,
        )
        cache[key] = np.asarray(result, dtype=np.float64)
        _touch_lru(cache, key)

    if negate:
        result = np.negative(result)
    return result



# =============================================================================

def create_dem_sampler(
    dem_list,
    # 'cache_dir' and 'q' kept only for backward compatibility; no longer used here
    cache_dir: str = None,
    q: float = None,
    *,
    out_dir: str = "logs",
    require_all: bool = False,
    na_fill: float = DEFAULT_NA_FILL,
    negate: bool = True,
):
    """
    Return a callable dem_sampler(points_xy) that uses stacked_dem_fill.
    Caching now lives inside stacked_dem_fill; this factory is intentionally thin.

    Parameters
    ----------
    dem_list : str | list[str]
        YAML path (list of DEM files) or a list of DEM file paths.
    cache_dir, q : deprecated
        Ignored; retained to avoid breaking older call sites.
    out_dir : str
        Output/log directory passed through to stacked_dem_fill (also hosts the cache).
    require_all : bool
        See stacked_dem_fill.
    na_fill : float | None
        See stacked_dem_fill.
    negate : bool
        If True, return depth (+down). Negation happens *outside* the cache in stacked_dem_fill.

    Returns
    -------
    callable
        dem_sampler(points_xy[, values=None]) -> np.ndarray
    """
    if isinstance(dem_list, str):
        with open(dem_list, "r") as f:
            dem_files = yaml.safe_load(f)
    else:
        dem_files = list(dem_list)

    def dem_sampler(points_xy, values=None):
        return stacked_dem_fill(
            files=dem_files,
            points=np.asarray(points_xy, dtype=float),
            out_dir=out_dir,
            values=values,
            negate=negate,
            require_all=require_all,
            na_fill=na_fill,
        )

    return dem_sampler


# =============================================================================


def ParseType(type):
    gdal_dt = gdal.GetDataTypeByName(type)
    if gdal_dt is GDT_Unknown:
        gdal_dt = GDT_Byte
    return gdal_dt


# =============================================================================


def bilinear(points, gt, raster):
    """Bilinear interpolated point  using data on raster
    points: npoint x 3 array of (x,y,z points)
    gt: geotransform information from gdal
    raster: the data
    """

    px = points[:, 0]
    py = points[:, 1]
    ny, nx = raster.shape

    # Half raster element widths
    hx = gt[1] / 2.0
    hy = gt[5] / 2.0
    # Calculate raster lower bound indices from point
    fx = (px - (gt[0] + hx)) / gt[1]
    fy = (py - (gt[3] + hy)) / gt[5]
    ix1 = np.floor(fx).astype(int)
    iy1 = np.floor(fy).astype(int)

    # Special case where point is on upper bounds
    # ix1=np.where(fx == float(nx - 1), ix1 - 1, ix1)
    # iy1=np.where(fy == float(ny - 1), iy1 -1, iy1)
    # Upper bound indices on raster
    ix2 = ix1 + 1
    iy2 = iy1 + 1

    # Test array bounds to ensure point is within raster midpoints
    inbounds = (ix1 >= 0) & (iy1 >= 0) & (ix2 <= nx - 1) & (iy2 <= ny - 1)

    #    return no_data
    inboundx = points[inbounds, 0]
    inboundy = points[inbounds, 1]

    ix1 = ix1[inbounds]
    ix2 = ix2[inbounds]
    iy1 = iy1[inbounds]
    iy2 = iy2[inbounds]

    # nearx = (np.abs(inboundx - 6.25143225e+005) <
    #          1.) & (np.abs(inboundy-4.20178364e+006) < 1.)

    # Calculate differences from point to bounding raster midpoints
    dx1 = inboundx - (gt[0] + ix1 * gt[1] + hx)
    dy1 = inboundy - (gt[3] + iy1 * gt[5] + hy)
    dx2 = (gt[0] + ix2 * gt[1] + hx) - inboundx
    dy2 = (gt[3] + iy2 * gt[5] + hy) - inboundy
    # Use the differences to weigh the four raster values
    div = gt[1] * gt[5]

    res = np.zeros_like(points[:, 0])
    res[:] = np.nan
    res[inbounds] = (
        raster[iy1, ix1] * dx2 * dy2 / div
        + raster[iy1, ix2] * dx1 * dy2 / div
        + raster[iy2, ix1] * dx2 * dy1 / div
        + raster[iy2, ix2] * dx1 * dy1 / div
    )
    ix1 = None
    ix2 = None
    iy1 = None
    iy2 = None

    return res

# ----------------------------- DEM sampler & cache ----------------------------


def create_dem_sampler(
    dem_list: Union[str, List[str]], cache_dir: str = "./dem_cache", q: float = 1e-6
):
    """Create a DEM sampler ``dem_sampler(xy)`` with disk caching.

    The cache key includes:
    - The list of DEM source names (from the YAML or the provided list).
    - A daily stamp (YYYYMMDD) to force refresh at least daily.
    - A compact signature of the request batch: length, quantized centroid, and
      quantized bounding box (min/max) of XY.

    Parameters
    ----------
    dem_list : str or list of str
        Either a YAML file path consumable by :func:`stacked_dem_fill` or a list
        of DEM file paths.
    cache_dir : str, optional
        Directory for :mod:`diskcache` storage, by default ``"./dem_cache"``.
    q : float, optional
        Quantization (meters) for the request signature, by default 1e-6.

    Returns
    -------
    callable
        ``dem_sampler(points_xy)`` mapping ``(k,2)`` XY→depth(+down) array.
    """
    if isinstance(dem_list, str):
        with open(dem_list, "r") as f:
            dem_files = yaml.safe_load(f)
    else:
        dem_files = list(dem_list)

    dem_names = tuple(map(str, dem_files))  # stable order if YAML preserves it
    dem_hash = hashlib.blake2b(
        "|".join(dem_names).encode("utf-8"), digest_size=12
    ).hexdigest()

    # Ensure cache_dir exists
    os.makedirs(cache_dir, exist_ok=True)
    # Ensure ./logs exists at the present level
    os.makedirs("logs", exist_ok=True)
    cache = dc.Cache(cache_dir)

    def _batch_signature(points_xy: np.ndarray, q=q) -> str:
        pts = np.asarray(points_xy, dtype=np.float64)
        n = pts.shape[0]
        if n == 0:
            return f"{dem_hash}|empty"
        # quantize
        pts_q = np.round(pts / q) * q
        cen = np.mean(pts_q, axis=0)
        mn = np.min(pts_q, axis=0)
        mx = np.max(pts_q, axis=0)
        sig_bytes = f"{n}|{cen[0]:.6f},{cen[1]:.6f}|{mn[0]:.6f},{mn[1]:.6f}|{mx[0]:.6f},{mx[1]:.6f}".encode(
            "utf-8"
        )
        sig_hash = hashlib.blake2b(sig_bytes, digest_size=12).hexdigest()
        # daily stamp ensures at least daily refresh
        day = datetime.utcnow().strftime("%Y%m%d")
        return f"{dem_hash}|{day}|{sig_hash}"

    def dem_sampler(points_xy: np.ndarray) -> np.ndarray:
        key = _batch_signature(points_xy)
        if key in cache:
            z = cache[key]
            return np.asarray(z, dtype=float)

        z = stacked_dem_fill(
            files=dem_files,
            points=np.asarray(points_xy, dtype=float),
            out_dir="logs",
            require_all=False,
            na_fill=2.0,
            negate=True,  # depth +down as used in the workflow
        )
        cache[key] = np.asarray(z, dtype=float)
        return np.asarray(z, dtype=float)

    return dem_sampler


# =================================================================
def test_main():
    points = np.zeros((4, 3), dtype=np.float64)
    points[0, :] = (546679, 4184973, np.nan)
    points[1, :] = (542385, 4184978, np.nan)
    points[2:,] = (616877, 4226454, np.nan)
    points[3:,] = (816877, 4226454, np.nan)

    print(
        stacked_dem_fill(
            ["bay_delta_mini.tif", "bay_delta_west.tif", "bay_delta_east.tif"],
            points[:, :2],
            "./",
            points[:, 2],
        )
    )


def filelist_main(demlistfile, pointlistfile, out_dir="./", sep=""):
    """higher level driver routine that parses the names of the DEMs from demfilelist
    points ( x,y or x,y,z ) from pointlistfile
    """

    filelist = [
        f.strip() for f in open(demlistfile, "r").readlines() if f and len(f) > 1
    ]
    bad_files = [f for f in filelist if not os.path.exists(f)]
    if len(bad_files) > 0:
        for fname in bad_files:
            print("Does not exist: %s " % fname)
        raise ValueError("Some file(s) not found")
    points = np.loadtxt(pointlistfile, delimiter=sep).copy()
    print(points)
    pshape = points.shape
    if pshape[1] == 2:
        z = np.zeros((points.shape[0], 1), dtype=np.float64)
        points = np.hstack((points, z))
        # points.resize((pshape[0],3),refcheck=False)
    points[:, 2] = np.nan
    print(points)
    points = stacked_dem_fill(filelist, points[:, :2], out_dir, points[:, 2])
    print(points)
    return points


def fill_2dm(infile, outfile, files, na_fill=DEFAULT_NA_FILL):
    import string

    print("Filling elevations in %s and writing to %s" % (infile, outfile))
    f = open(infile, "r")
    all_lines = f.readlines()
    f.close()
    firstndx = next((i for i, v in enumerate(all_lines) if v.startswith("ND")), -1)
    endndx = (
        next(
            (i for i, v in enumerate(all_lines[firstndx:]) if not v.startswith("ND")),
            -1,
        )
        + firstndx
    )
    pts = np.array(
        [
            [float(x) for x in string.split(s.strip())[2:5]]
            for s in all_lines[firstndx:endndx]
        ]
    )
    pts[:, 2] = np.nan
    values = stacked_dem_fill(
        files,
        pts[:, :2],
        os.path.dirname(outfile),
        pts[:, 2],
        require_all=False,
        na_fill=na_fill,
    )

    all_lines[firstndx:endndx] = (
        "ND %i  %18.8f %18.8f  %18.8f\n" % ((i + 1), pt[0], pt[1], val)
        for i, (pt, val) in enumerate(zip(pts, values))
    )
    fout = open(outfile, "w")
    for line in all_lines:
        fout.write(line)
    return None


def fill_gr3(infile, outfile, files, elev2depth=True, na_fill=DEFAULT_NA_FILL):
    import string
    import math

    print("Filling elevations in %s and writing to %s" % (infile, outfile))
    f = open(infile, "r")
    all_lines = f.readlines()
    f.close()
    # read the first line of text and convert nelement and nnode to integers
    nelement, nnode = [int(x) for x in all_lines[1].split()[0:2]]
    firstndx = 2
    stopndx = nnode + 2
    pts = np.array(
        [
            [float(x) for x in string.split(s.strip())[1:4]]
            for s in all_lines[firstndx:stopndx]
        ]
    )
    pts[:, 2] = np.nan
    values = stacked_dem_fill(
        files,
        pts[:, :2],
        os.path.dirname(outfile),
        pts[:, 2],
        require_all=False,
        na_fill=na_fill,
        negate=elev2depth,
    )
    padding = 2
    maxnum = int(math.log10(max(nelement, nnode))) + padding
    # ifmt = "%" + ("%ii" % maxnum)
    ifmtj = "%-" + ("%ii" % maxnum)
    ffmt = "%18.8f"
    nfmt = ifmtj + ffmt * 3 + "\n"

    all_lines[firstndx:stopndx] = (
        nfmt % ((i + 1), pt[0], pt[1], val)
        for i, (pt, val) in enumerate(zip(pts, values))
    )

    fout = open(outfile, "w")
    for line in all_lines:
        fout.write(line)
    return None


def create_arg_parser():
    import argparse

    parser = argparse.ArgumentParser(
        description="Fill node elevations in a .2dm SMS mesh or gr3 file using a prioritized list of DEMs."
    )
    parser.add_argument(dest="filename", default=None, help="name of 2dm or gr3 file")
    parser.add_argument(
        "demfile",
        default=None,
        help="file containing list of DEMs. These can be in any form that gdal accepts, which includes ESRI ascii format and GeoTiffs",
    )
    parser.add_argument(
        "--elev2depth",
        action="store_true",
        default=False,
        help="Convert elevation to depth by flipping sign. This is typical when using gr3 format, less so with 2dm.",
    )
    parser.add_argument(
        "--fill",
        default=DEFAULT_NA_FILL,
        help="Fill value for areas not covered by supplied rasters.",
    )
    return parser


def main():
    import shutil

    parser = create_arg_parser()
    args = parser.parse_args()
    infile = args.filename
    bakfile = infile + ".bak"
    print("Backing up file to %s" % bakfile)
    shutil.copyfile(infile, bakfile)
    demlistfile = args.demfile

    if demlistfile.endswith(".txt"):
        files = [
            x.strip()
            for x in open(demlistfile, "r").readlines()
            if x and len(x) > 1 and not x.startswith("#")
        ]
    elif demlistfile.endswith(".yaml"):
        import yaml

        with open(demlistfile, "r") as f:
            files = yaml.load(f)
    else:
        raise ValueError("Not supported DEM list file extension")

    elev2depth = args.elev2depth
    na_fill = args.fill
    if infile.endswith(".2dm"):
        if elev2depth:
            print("Warning: --elev2depth is an unusual option for 2dm files")
        fill_2dm(bakfile, infile, files, na_fill=na_fill)
    elif infile.endswith(".gr3"):
        if not elev2depth:
            print("Warning: omitting --elev2depth is an unusual option for gr3 files")
        fill_gr3(bakfile, infile, files, elev2depth, na_fill=na_fill)
    else:
        raise ValueError("Input file format not recognized (no gr3 or 2dm extension")


if __name__ == "__main__":
    main()
