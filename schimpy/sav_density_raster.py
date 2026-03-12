"""
Raster-based sampling utilities for submerged and emergent vegetation attributes.

This module supports two related preprocessing tasks for SCHISM nodes:

1. Assign a vegetation density value from categorical raster products.
2. Assign a vegetation canopy height from categorical raster products.

The raster products are treated as integer-valued thematic layers rather than
continuous scalar fields. This distinction is important.

For density, the algorithm uses two rasters:

- a vegetation type raster
- a vegetation density-bin raster

The vegetation type at a node is taken from the pixel containing the node.
The density-bin at a node is also taken from the pixel containing the node.
A local stencil, typically 3x3 pixels, is then sampled around the node.
Within that stencil, only pixels that match both the center vegetation type
and the center density bin are considered "like-classified" and contribute
to the final density value. This avoids averaging unlike categories together.

For height, the categorical rasters are sampled with 1x1 support. Height is
derived from the node's vegetation type, the node's density bin, and the node
depth. No neighborhood averaging is applied to height because the type and
density rasters are categorical and averaging categories produces ambiguous
values.

Coverage behavior is controlled by ``strict``:

- ``strict=True`` raises an exception when requested samples fall outside raster
  coverage.
- ``strict=False`` replaces uncovered or nodata samples with zero. In practice,
  this means "no vegetation" for the categorical products used here.

The module is intended for use from SCHISM preprocessing YAML expressions, but
the lower-level functions are also suitable for direct testing and scripted use.

Notes
-----
The algorithms in this module assume:

- the vegetation rasters use integer category codes;
- type code 0 means no vegetation;
- density bin 0 means no density / no vegetation;
- mesh node ``z`` is positive downward when height or depth-based screening
  is applied.

The density mapping from ``(type, density_bin)`` to final density value is
supplied either as a CSV file or as a 2D lookup table.
"""

import csv
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import rowcol
from rasterio.windows import Window


def _normalize_stencil(stencil):
    """
    Normalize a stencil specification to odd integer dimensions.

    Parameters
    ----------
    stencil : int or iterable of length 2
        Neighborhood shape specification. An integer ``n`` means an ``n x n``
        stencil. A two-element iterable gives ``(n_rows, n_cols)``.

    Returns
    -------
    sy : int
        Number of stencil rows.
    sx : int
        Number of stencil columns.

    Raises
    ------
    ValueError
        Raised when the stencil is not an integer or length-2 iterable, when
        either dimension is smaller than 1, or when either dimension is even.

    Notes
    -----
    Odd dimensions are required so that the stencil has a single, unambiguous
    center cell aligned with the pixel containing the node.

    Examples
    --------
    ``3`` becomes ``(3, 3)``.

    ``(3, 5)`` is kept as ``(3, 5)``.
    """    
    if isinstance(stencil, int):
        sy = sx = int(stencil)
    else:
        if len(stencil) != 2:
            raise ValueError("stencil must be an int or length-2 iterable")
        sy, sx = int(stencil[0]), int(stencil[1])

    if sy < 1 or sx < 1:
        raise ValueError("stencil dimensions must be >= 1")
    if (sy % 2) == 0 or (sx % 2) == 0:
        raise ValueError("stencil dimensions must be odd")

    return sy, sx



def _lut_from_csv(path):
    """
    Read a complete density lookup table from CSV.

    Parameters
    ----------
    path : str or path-like
        CSV file with exactly the columns ``type``, ``density_bin``, and
        ``density``.

    Returns
    -------
    lut : ndarray of float, shape (n_type, n_bin)
        Lookup table such that ``lut[type_code, density_bin]`` returns the
        mapped density value.

    Raises
    ------
    FileNotFoundError
        Raised when the CSV file does not exist.
    ValueError
        Raised when the CSV columns are incorrect, when the file is empty,
        when duplicate rows are present, or when the grid of type/bin
        combinations is incomplete.

    Notes
    -----
    The CSV is interpreted as a complete table, not a sparse one. Every
    ``(type, density_bin)`` combination from zero up to the maximum values
    present in the file must appear exactly once. This avoids silent defaults
    and keeps the mapping explicit.

    A typical table might encode rules such as:

    - no vegetation -> 0
    - submerged, light -> 5
    - submerged, medium -> 10
    - submerged, heavy -> 25
    - floating, heavy -> 15
    - emergent, heavy -> 30
    """    
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"mapping csv not found: {path}")

    rows = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        expected = {"type", "density_bin", "density"}
        if set(reader.fieldnames or []) != expected:
            raise ValueError(
                f"mapping csv must have exactly columns {sorted(expected)}, "
                f"got {reader.fieldnames}"
            )
        for row in reader:
            rows.append(
                (
                    int(row["type"]),
                    int(row["density_bin"]),
                    float(row["density"]),
                )
            )

    if not rows:
        raise ValueError("mapping csv is empty")

    max_type = max(r[0] for r in rows)
    max_bin = max(r[1] for r in rows)

    lut = np.zeros((max_type + 1, max_bin + 1), dtype=float)
    seen = set()
    for t, b, d in rows:
        key = (t, b)
        if key in seen:
            raise ValueError(f"duplicate mapping row for type={t}, density_bin={b}")
        seen.add(key)
        lut[t, b] = d

    for t in range(max_type + 1):
        for b in range(max_bin + 1):
            if (t, b) not in seen:
                raise ValueError(f"missing mapping row for type={t}, density_bin={b}")

    return lut


def _lut_from_matrix(mapping):
    """
    Convert an in-memory mapping matrix to a lookup table.

    Parameters
    ----------
    mapping : array-like, shape (n_type, n_bin)
        Two-dimensional lookup table indexed by ``[type, density_bin]``.

    Returns
    -------
    lut : ndarray of float
        The mapping converted to a NumPy array.

    Raises
    ------
    ValueError
        Raised when the supplied mapping is not two-dimensional.

    Notes
    -----
    This form is convenient in tests and in direct Python use. In YAML-driven
    preprocessing, the CSV form is often easier because the expression string
    can pass a scalar file path.
    """    
    lut = np.asarray(mapping, dtype=float)
    if lut.ndim != 2:
        raise ValueError("mapping matrix must be 2D")
    return lut


def _read_bounding_window(path, rows, cols, *, pad=0, band=1, strict=True):
    """
    Read the smallest raster window that covers a set of row/column locations.

    Parameters
    ----------
    path : str or path-like
        Raster path.
    rows : array-like of int
        Raster row indices to cover.
    cols : array-like of int
        Raster column indices to cover.
    pad : int, optional
        Extra cells added on all sides of the requested bounding box. This is
        typically the stencil half-width.
    band : int, optional
        Raster band number.
    strict : bool, optional
        Controls behavior when the requested window lies outside raster
        coverage. If True, the function raises an exception. If False, it
        returns ``(None, None, None)``.

    Returns
    -------
    arr : ndarray or None
        Raster data for the clipped window. Nodata values are converted to
        ``np.nan``. If ``strict=False`` and the requested window has no
        overlap with the raster, ``arr`` is ``None``.
    r0 : int or None
        Row index of the window origin in full-raster coordinates.
    c0 : int or None
        Column index of the window origin in full-raster coordinates.

    Raises
    ------
    ValueError
        Raised when no points are supplied or when the requested window has no
        overlap with the raster and ``strict=True``.

    Notes
    -----
    This helper is the basis for chunked sampling. The calling code converts
    full-raster row/column indices to local indices within ``arr`` by
    subtracting ``r0`` and ``c0``.

    The window is clipped to raster bounds before reading. This keeps I/O local
    to the points being sampled and avoids reading the full raster into memory.
    """    
    rows = np.asarray(rows, dtype=int)
    cols = np.asarray(cols, dtype=int)

    if rows.size == 0 or cols.size == 0:
        if strict:
            raise ValueError("no points provided to window reader")
        return None, None, None

    with rasterio.open(path) as src:
        rmin = int(rows.min()) - pad
        rmax = int(rows.max()) + pad
        cmin = int(cols.min()) - pad
        cmax = int(cols.max()) + pad

        rmin_clip = max(rmin, 0)
        cmin_clip = max(cmin, 0)
        rmax_clip = min(rmax, src.height - 1)
        cmax_clip = min(cmax, src.width - 1)

        h = rmax_clip - rmin_clip + 1
        w = cmax_clip - cmin_clip + 1

        if h <= 0 or w <= 0:
            if strict:
                raise ValueError("empty bounding window; points may be outside raster extent")
            return None, None, None

        win = Window(cmin_clip, rmin_clip, w, h)
        arr = src.read(band, window=win)

        nodata = src.nodata
        if nodata is not None:
            arr = arr.astype(float, copy=False)
            arr[arr == nodata] = np.nan

    return arr, rmin_clip, cmin_clip


def sample_raster_stencil_int(
    xy,
    raster_file,
    *,
    stencil=(3, 3),
    nan_to=0,
    strict=False,
    chunk_size=200_000,
):
    """
    Sample an integer-valued raster at node centers and over a local stencil.

    Parameters
    ----------
    xy : ndarray, shape (n, 2)
        Node coordinates in the raster coordinate system, given as ``(x, y)``.
    raster_file : str or path-like
        Path to the raster to sample.
    stencil : int or (int, int), optional
        Stencil size centered on the containing pixel. An integer gives a
        square stencil. A two-element iterable gives ``(n_rows, n_cols)``.
        Dimensions must be odd.
    nan_to : int, optional
        Integer value used in place of nodata or uncovered samples when
        ``strict=False``.
    strict : bool, optional
        Controls behavior outside raster coverage.

        - If True, the function raises an exception when a requested sample
          cannot be represented by the raster.
        - If False, uncovered center pixels and uncovered stencil cells are
          replaced with ``nan_to``.

    chunk_size : int, optional
        Number of points processed per chunk.

    Returns
    -------
    center_out : ndarray of int, shape (n,)
        Raster values for the pixel containing each point.
    stencil_out : ndarray of int, shape (n, sy*sx)
        Flattened stencil values in row-major order for each point.
    rows_out : ndarray of int, shape (n,)
        Raster row index of the containing pixel for each point.
    cols_out : ndarray of int, shape (n,)
        Raster column index of the containing pixel for each point.

    Raises
    ------
    ValueError
        Raised when ``xy`` does not have shape ``(n, 2)``, when stencil
        dimensions are invalid, or when raster coverage is insufficient and
        ``strict=True``.

    Notes
    -----
    Sampling is based on the pixel containing the point, not bilinear
    interpolation. This is intentional.

    The rasters used here are thematic rasters with integer class codes.
    Interpolating between codes such as "submerged", "floating", and
    "emergent" would produce values with no categorical meaning. The correct
    operation for these products is to identify the containing pixel and then
    inspect neighboring pixels explicitly.

    The returned stencil is flattened in row-major order. For a 3x3 stencil,
    index 4 is the center cell.

    Behavior when ``strict=False``
    ------------------------------
    When coverage is incomplete, the function substitutes ``nan_to``. In the
    vegetation workflows in this module, ``nan_to=0`` is used so that missing
    coverage is treated as "no vegetation".

    Performance
    -----------
    Points are processed in chunks, and each chunk reads a single bounding
    window from the raster. This keeps disk access local and avoids reading the
    full raster for every point.
    """
    xy = np.asarray(xy, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError("xy must be shape (n, 2)")

    sy, sx = _normalize_stencil(stencil)
    pad = max(sy // 2, sx // 2)
    n = xy.shape[0]

    center_out = np.full(n, nan_to, dtype=np.int64)
    stencil_out = np.full((n, sy * sx), nan_to, dtype=np.int64)
    rows_out = np.empty(n, dtype=int)
    cols_out = np.empty(n, dtype=int)

    with rasterio.open(raster_file) as src:
        xs = xy[:, 0]
        ys = xy[:, 1]
        rows, cols = rowcol(src.transform, xs, ys, op=np.floor)
        rows = np.asarray(rows, dtype=int)
        cols = np.asarray(cols, dtype=int)

    rows_out[:] = rows
    cols_out[:] = cols

    oy = np.arange(-(sy // 2), (sy // 2) + 1, dtype=int)
    ox = np.arange(-(sx // 2), (sx // 2) + 1, dtype=int)

    for i0 in range(0, n, chunk_size):
        i1 = min(i0 + chunk_size, n)

        rr = rows[i0:i1]
        cc = cols[i0:i1]

        arr, r0, c0 = _read_bounding_window(
            raster_file, rr, cc, pad=pad, strict=strict
        )

        if arr is None:
            continue

        rr_local = rr - r0
        cc_local = cc - c0

        h, w = arr.shape
        center_inb = (
            (rr_local >= 0) & (rr_local < h) &
            (cc_local >= 0) & (cc_local < w)
        )

        if strict and not np.all(center_inb):
            raise ValueError("some center points lie outside raster coverage")

        cen = np.full(rr_local.shape, nan_to, dtype=float)
        cen[center_inb] = arr[rr_local[center_inb], cc_local[center_inb]]
        cen = np.where(np.isfinite(cen), cen, nan_to).astype(np.int64, copy=False)
        center_out[i0:i1] = cen

        rr3 = rr_local[:, None, None] + oy[None, :, None]
        cc3 = cc_local[:, None, None] + ox[None, None, :]

        rr3_full = np.broadcast_to(rr3, (rr_local.shape[0], sy, sx))
        cc3_full = np.broadcast_to(cc3, (cc_local.shape[0], sy, sx))

        inb = (
            (rr3_full >= 0) & (rr3_full < h) &
            (cc3_full >= 0) & (cc3_full < w)
        )

        if strict:
            center_pos = (sy // 2, sx // 2)
            if not np.all(inb[:, center_pos[0], center_pos[1]]):
                raise ValueError("some center stencil cells lie outside raster coverage")

        vals = np.full((rr_local.shape[0], sy, sx), np.nan, dtype=float)
        vals[inb] = arr[rr3_full[inb], cc3_full[inb]]
        vals = np.where(np.isfinite(vals), vals, nan_to).astype(np.int64, copy=False)

        stencil_out[i0:i1, :] = vals.reshape(vals.shape[0], sy * sx)

    return center_out, stencil_out, rows_out, cols_out


def compute_like_class_mean_density(center_type, center_bin, type_stencil, bin_stencil, lut):
    """
    Compute vegetation density from center classes and local stencil samples.

    Parameters
    ----------
    center_type : ndarray of int, shape (n,)
        Vegetation type at the node center.
    center_bin : ndarray of int, shape (n,)
        Density bin at the node center.
    type_stencil : ndarray of int, shape (n, k)
        Vegetation type sampled over the local stencil for each node.
    bin_stencil : ndarray of int, shape (n, k)
        Density bin sampled over the local stencil for each node.
    lut : ndarray of float, shape (n_type, n_bin)
        Lookup table mapping ``(type, density_bin)`` to final density values.

    Returns
    -------
    dens : ndarray of float, shape (n,)
        Density assigned to each node.
    like : ndarray of bool, shape (n, k)
        Mask indicating which stencil cells matched both the center vegetation
        type and the center density bin.

    Raises
    ------
    ValueError
        Raised when the input array dimensions do not align, when negative class
        codes are supplied, or when class codes exceed the bounds of ``lut``.

    Notes
    -----
    This function implements the central density algorithm in the module.

    The node is first classified by its center vegetation type and center
    density bin. A stencil around the node is then examined. Only stencil cells
    that match both the center type and the center density bin are included in
    the neighborhood average.

    This "like-classified" average avoids mixing unlike vegetation classes.
    For example, a node classified as submerged vegetation with medium density
    is not influenced by adjacent floating or emergent pixels, and it is not
    influenced by submerged pixels in a different density bin.

    Nodes with either ``center_type == 0`` or ``center_bin == 0`` are treated
    as inactive and assigned density 0.

    The lookup table is applied to stencil classes before averaging. This makes
    the local average an average of final density values rather than an average
    of raw density-bin integers.
    """
    center_type = np.asarray(center_type, dtype=int)
    center_bin = np.asarray(center_bin, dtype=int)
    type_stencil = np.asarray(type_stencil, dtype=int)
    bin_stencil = np.asarray(bin_stencil, dtype=int)
    lut = np.asarray(lut, dtype=float)

    if center_type.ndim != 1 or center_bin.ndim != 1:
        raise ValueError("center_type and center_bin must be 1D")
    if type_stencil.ndim != 2 or bin_stencil.ndim != 2:
        raise ValueError("type_stencil and bin_stencil must be 2D")
    if type_stencil.shape != bin_stencil.shape:
        raise ValueError("type_stencil and bin_stencil must have same shape")
    if type_stencil.shape[0] != center_type.shape[0]:
        raise ValueError("stencil arrays must match number of centers")

    if np.any(center_type < 0) or np.any(center_bin < 0):
        raise ValueError("negative class/bin values are not allowed")
    if np.any(type_stencil < 0) or np.any(bin_stencil < 0):
        raise ValueError("negative stencil class/bin values are not allowed")

    if center_type.max(initial=0) >= lut.shape[0]:
        raise ValueError("center_type exceeds mapping LUT bounds")
    if center_bin.max(initial=0) >= lut.shape[1]:
        raise ValueError("center_bin exceeds mapping LUT bounds")
    if type_stencil.max(initial=0) >= lut.shape[0]:
        raise ValueError("type_stencil exceeds mapping LUT bounds")
    if bin_stencil.max(initial=0) >= lut.shape[1]:
        raise ValueError("bin_stencil exceeds mapping LUT bounds")

    active = (center_type > 0) & (center_bin > 0)

    like = (
        active[:, None]
        & (type_stencil == center_type[:, None])
        & (bin_stencil == center_bin[:, None])
    )

    dens_stencil = lut[type_stencil, bin_stencil]
    dens_like = np.where(like, dens_stencil, np.nan)

    mask_like = np.isfinite(dens_like)
    count = mask_like.sum(axis=1)
    sum_ = np.nansum(dens_like, axis=1)

    mean_like = np.full(sum_.shape, np.nan, dtype=float)
    have = count > 0
    mean_like[have] = sum_[have] / count[have]

    center_density = lut[center_type, center_bin]
    dens = np.where(np.isfinite(mean_like), mean_like, center_density)
    dens = np.where(active, dens, 0.0)

    return dens, like

def sav_node_classes(
    mesh,
    nodes_sel,
    veg_type_tif,
    veg_density_tif,
    *,
    strict=False,
):
    """
    Return the vegetation type and density-bin seen at selected mesh nodes.

    Parameters
    ----------
    mesh : object
        Mesh-like object providing ``mesh.nodes[:, :3]`` as ``x, y, z``.
    nodes_sel : array-like of int
        Indices of nodes to sample.
    veg_type_tif : str or path-like
        Vegetation type raster.
    veg_density_tif : str or path-like
        Vegetation density-bin raster.
    strict : bool, optional
        Controls how uncovered raster samples are handled.

    Returns
    -------
    type_center : ndarray of int, shape (n,)
        Vegetation type at each selected node.
    bin_center : ndarray of int, shape (n,)
        Density bin at each selected node.

    Notes
    -----
    This helper exposes the node-centered categorical values directly. It is
    useful when checking which class a node receives, and it is the natural
    basis for canopy-height calculations where classes should remain discrete.
    """
    nodes_sel = np.asarray(nodes_sel, dtype=int)
    xyz = np.asarray(mesh.nodes[nodes_sel, :3], dtype=float)

    xs = xyz[:, 0]
    ys = xyz[:, 1]

    return sample_center_pixels(
        xs, ys, veg_type_tif, veg_density_tif, strict=strict
    )

def compute_height(
    z,
    center_type,
    center_bin,
    *,
    sav_fractions=(0.0, 0.3, 0.65, 0.9),
    sav_min_height=0.0,
    sav_max_height=None,
    emergent_height=100.0,
    floating_height=100.0,
    return_debug=False,
):
    """
    Compute canopy height from node depth and node-centered vegetation classes.

    Parameters
    ----------
    z : array-like of float, shape (n,)
        Node depth values. Positive values indicate deeper nodes.
    center_type : array-like of int, shape (n,)
        Vegetation type at the node center.
    center_bin : array-like of int, shape (n,)
        Density bin at the node center.
    sav_fractions : sequence of float, optional
        Multipliers used for submerged aquatic vegetation, indexed by density
        bin. For example, ``(0.0, 0.3, 0.65, 0.9)`` means light, medium, and
        heavy SAV heights are 30, 65, and 90 percent of node depth.
    sav_min_height : float, optional
        Lower bound applied to submerged vegetation height after the fraction
        rule is applied.
    sav_max_height : float or None, optional
        Upper bound applied to submerged vegetation height after the fraction
        rule is applied. If None, no upper bound is imposed.
    emergent_height : float, optional
        Height assigned to emergent vegetation.
    floating_height : float, optional
        Height assigned to floating vegetation.
    return_debug : bool, optional
        If True, return both the height array and a dictionary of intermediate
        masks and inputs.

    Returns
    -------
    height : ndarray of float, shape (n,)
        Canopy height at each node.
    debug : dict, optional
        Returned only when ``return_debug=True``.

    Raises
    ------
    ValueError
        Raised when input shapes do not match, when negative class codes are
        supplied, or when ``sav_fractions`` does not cover the supplied
        density bins.

    Notes
    -----
    Height is derived from node-centered classes using discrete rules.

    Type handling
    -------------
    - type 0: no vegetation -> height 0
    - type 1: submerged vegetation -> fraction of depth, optionally bounded by
      ``sav_min_height`` and ``sav_max_height``
    - type 2: floating vegetation -> ``floating_height``
    - type 3: emergent vegetation -> ``emergent_height``

    This function uses only the node-centered class values. No neighborhood
    averaging is used for height because height depends on the local category
    assignment, not on a mixture of nearby categories.
    """
    z = np.asarray(z, dtype=float)
    center_type = np.asarray(center_type, dtype=int)
    center_bin = np.asarray(center_bin, dtype=int)

    if z.ndim != 1 or center_type.ndim != 1 or center_bin.ndim != 1:
        raise ValueError("z, center_type, center_bin must be 1D")
    if not (z.shape == center_type.shape == center_bin.shape):
        raise ValueError("z, center_type, center_bin must have same shape")
    if np.any(center_type < 0) or np.any(center_bin < 0):
        raise ValueError("negative class/bin values are not allowed")

    frac = np.asarray(sav_fractions, dtype=float)
    if frac.ndim != 1:
        raise ValueError("sav_fractions must be 1D")
    if np.any(frac < 0.0):
        raise ValueError("sav_fractions must be nonnegative")
    if frac.shape[0] <= center_bin.max(initial=0):
        raise ValueError("sav_fractions does not cover all density bins")

    if sav_min_height is not None:
        sav_min_height = float(sav_min_height)
        if sav_min_height < 0.0:
            raise ValueError("sav_min_height must be nonnegative")

    if sav_max_height is not None:
        sav_max_height = float(sav_max_height)
        if sav_max_height < 0.0:
            raise ValueError("sav_max_height must be nonnegative")

    if (
        sav_min_height is not None
        and sav_max_height is not None
        and sav_min_height > sav_max_height
    ):
        raise ValueError("sav_min_height cannot exceed sav_max_height")

    height = np.zeros_like(z, dtype=float)

    is_none = center_type == 0
    is_sav = center_type == 1
    is_float = center_type == 2
    is_emergent = center_type == 3

    unknown_type = ~(is_none | is_sav | is_float | is_emergent)
    if np.any(unknown_type):
        bad = np.unique(center_type[unknown_type])
        raise ValueError(f"unsupported vegetation type code(s): {bad.tolist()}")

    # Floating and emergent are assigned fixed large canopy heights.
    height[is_float] = float(floating_height)
    height[is_emergent] = float(emergent_height)

    # SAV height is a fraction of depth, optionally bounded.
    if np.any(is_sav):
        sav_h = frac[center_bin[is_sav]] * z[is_sav]

        # Negative depths do not produce negative canopy heights.
        sav_h = np.maximum(sav_h, 0.0)

        if sav_min_height is not None:
            active_sav = center_bin[is_sav] > 0
            sav_h = np.where(active_sav, np.maximum(sav_h, sav_min_height), 0.0)

        if sav_max_height is not None:
            sav_h = np.minimum(sav_h, sav_max_height)

        height[is_sav] = sav_h

    # Keep all heights nonnegative.
    height = np.maximum(height, 0.0)

    if return_debug:
        return height, {
            "type_center": center_type,
            "bin_center": center_bin,
            "is_none": is_none,
            "is_sav": is_sav,
            "is_float": is_float,
            "is_emergent": is_emergent,
            "sav_fractions": frac,
            "sav_min_height": sav_min_height,
            "sav_max_height": sav_max_height,
        }

    return height

def sav_height(
    mesh,
    nodes_sel,
    veg_type_tif,
    veg_density_tif,
    *,
    sav_fractions=(0.0, 0.3, 0.65, 0.9),
    sav_min_height=0.0,
    sav_max_height=None,
    emergent_height=100.0,
    floating_height=100.0,
    strict=False,
    return_debug=False,
):
    """
    Compute canopy height for selected mesh nodes from categorical vegetation rasters.

    Parameters
    ----------
    mesh : object
        Mesh-like object providing ``mesh.nodes[:, :3]`` as ``x, y, z``.
    nodes_sel : array-like of int
        Indices of nodes to sample.
    veg_type_tif : str or path-like
        Vegetation type raster.
    veg_density_tif : str or path-like
        Vegetation density-bin raster.
    sav_fractions : sequence of float, optional
        Multipliers applied to depth for submerged vegetation, indexed by
        density bin.
    sav_min_height : float, optional
        Lower bound for submerged vegetation height.
    sav_max_height : float or None, optional
        Upper bound for submerged vegetation height.
    emergent_height : float, optional
        Height assigned to emergent vegetation.
    floating_height : float, optional
        Height assigned to floating vegetation.
    strict : bool, optional
        Controls how uncovered raster samples are handled.
    return_debug : bool, optional
        If True, return both the height array and a dictionary of intermediate
        values.

    Returns
    -------
    height : ndarray of float, shape (n,)
        Height assigned to the selected nodes.
    debug : dict, optional
        Returned only when ``return_debug=True``.

    Notes
    -----
    This function samples the type and density-bin rasters with 1x1 support and
    then computes height from those node-centered classes. This is intentional.

    The vegetation rasters are categorical. For canopy height, the correct
    question is "what class does this node see?" rather than "what is the
    average category nearby?" Averaging categories over a larger support would
    blur boundaries and is not appropriate for this use.

    This function is the height analogue of :func:`sav_density`, but unlike
    density it does not use a neighborhood average.
    """
    nodes_sel = np.asarray(nodes_sel, dtype=int)
    xyz = np.asarray(mesh.nodes[nodes_sel, :3], dtype=float)

    xs = xyz[:, 0]
    ys = xyz[:, 1]
    z = xyz[:, 2]

    center_type, center_bin = sample_center_pixels(
        xs, ys, veg_type_tif, veg_density_tif, strict=strict
    )

    out = compute_height(
        z,
        center_type,
        center_bin,
        sav_fractions=sav_fractions,
        emergent_height=emergent_height,
        floating_height=floating_height,
        sav_max_height=sav_max_height,
        return_debug=return_debug,
    )

    return out


def compute_density(
    xs,
    ys,
    z,
    veg_type_tif,
    veg_density_tif,
    *,
    mapping,
    stencil=3,
    depth_limit=3.0,
    strict=False,
    chunk_size=200_000,
    return_debug=False,
):
    """
    Compute vegetation density at points from categorical vegetation rasters.

    Parameters
    ----------
    xs, ys : array-like of float
        Point coordinates in the raster coordinate system.
    z : array-like of float
        Node depth values. Positive values indicate deeper nodes.
    veg_type_tif : str or path-like
        Vegetation type raster.
    veg_density_tif : str or path-like
        Vegetation density-bin raster.
    mapping : str, path-like, or 2D array-like
        Mapping from ``(type, density_bin)`` to final density. This may be a
        CSV file or an in-memory lookup table.
    stencil : int or (int, int), optional
        Neighborhood size used for the local density average.
    depth_limit : float, optional
        Nodes with ``z > depth_limit`` are assigned density 0 after raster
        sampling and lookup.
    strict : bool, optional
        Controls how uncovered raster samples are handled.
    chunk_size : int, optional
        Number of points processed per chunk during raster sampling.
    return_debug : bool, optional
        If True, return both the density array and a dictionary containing the
        sampled center classes, stencil classes, like-mask, depth, and lookup
        table.

    Returns
    -------
    dens : ndarray of float, shape (n,)
        Density at each point.
    debug : dict, optional
        Returned only when ``return_debug=True``.

    Raises
    ------
    ValueError
        Raised when ``xs``, ``ys``, and ``z`` do not have the same shape, when
        stencil dimensions are invalid, or when mapping inputs are invalid.

    Notes
    -----
    The raster sampling and the density calculation are intentionally separated:

    - vegetation type and density bin are sampled as categorical integers;
    - the final density value is computed from the lookup table and the
      like-classified stencil average.

    This function is appropriate when coordinates are already available as
    arrays. For mesh-based workflows, use :func:`sav_density`.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    z = np.asarray(z, dtype=float)

    if xs.shape != ys.shape or xs.shape != z.shape:
        raise ValueError("xs, ys, z must have same shape")

    xy = np.column_stack([xs, ys])

    if isinstance(mapping, (str, Path)):
        lut = _lut_from_csv(mapping)
    else:
        lut = _lut_from_matrix(mapping)

    sy, sx = _normalize_stencil(stencil)

    type_center, type_stencil, _, _ = sample_raster_stencil_int(
        xy,
        veg_type_tif,
        stencil=(sy, sx),
        nan_to=0,
        strict=strict,
        chunk_size=chunk_size,
    )
    bin_center, bin_stencil, _, _ = sample_raster_stencil_int(
        xy,
        veg_density_tif,
        stencil=(sy, sx),
        nan_to=0,
        strict=strict,
        chunk_size=chunk_size,
    )

    dens, like = compute_like_class_mean_density(
        type_center,
        bin_center,
        type_stencil,
        bin_stencil,
        lut,
    )

    dens = np.where(z > float(depth_limit), 0.0, dens)

    if return_debug:
        return dens, {
            "type_center": type_center,
            "bin_center": bin_center,
            "type_stencil": type_stencil,
            "bin_stencil": bin_stencil,
            "like_mask": like,
            "depth": z,
            "lut": lut,
        }

    return dens


def sav_density(
    mesh,
    nodes_sel,
    veg_type_tif,
    veg_density_tif,
    *,
    stencil=3,
    depth_limit=3.0,
    mapping,
    strict=False,
    chunk_size=200_000,
    return_debug=False,
):
    """
    Compute vegetation density for a selected set of mesh nodes.

    Parameters
    ----------
    mesh : object
        Mesh-like object providing ``mesh.nodes[:, :3]`` as ``x, y, z``.
    nodes_sel : array-like of int
        Indices of nodes to sample.
    veg_type_tif : str or path-like
        Vegetation type raster.
    veg_density_tif : str or path-like
        Vegetation density-bin raster.
    stencil : int or (int, int), optional
        Neighborhood size used for the local density average.
    depth_limit : float, optional
        Nodes deeper than this threshold are assigned density 0.
    mapping : str, path-like, or 2D array-like
        Mapping from ``(type, density_bin)`` to final density.
    strict : bool, optional
        Controls how uncovered raster samples are handled.
    chunk_size : int, optional
        Number of points processed per chunk during raster sampling.
    return_debug : bool, optional
        If True, return both the density array and a dictionary of intermediate
        values.

    Returns
    -------
    dens : ndarray of float, shape (n,)
        Density assigned to the selected nodes.
    debug : dict, optional
        Returned only when ``return_debug=True``.

    Notes
    -----
    This is the YAML-facing wrapper used in SCHISM preprocessing. It extracts
    the selected node coordinates from the mesh and delegates the calculation
    to :func:`compute_density`.

    Typical use from preprocessing YAML is to call this function directly in an
    ``attribute`` expression and pass the mapping as a CSV file path.
    """
    nodes_sel = np.asarray(nodes_sel, dtype=int)
    xyz = np.asarray(mesh.nodes[nodes_sel, :3], dtype=float)

    xs = xyz[:, 0]
    ys = xyz[:, 1]
    z = xyz[:, 2]

    return compute_density(
        xs,
        ys,
        z,
        veg_type_tif,
        veg_density_tif,
        mapping=mapping,
        stencil=stencil,
        depth_limit=depth_limit,
        strict=strict,
        chunk_size=chunk_size,
        return_debug=return_debug,
    )

def sample_center_pixels(xs, ys, veg_type_tif, veg_density_tif, *, strict=False):
    """
    Sample vegetation type and density-bin at point centers using 1x1 support.

    Parameters
    ----------
    xs, ys : array-like of float
        Point coordinates in the raster coordinate system.
    veg_type_tif : str or path-like
        Vegetation type raster.
    veg_density_tif : str or path-like
        Vegetation density-bin raster.
    strict : bool, optional
        Controls how uncovered raster samples are handled.

    Returns
    -------
    type_center : ndarray of int, shape (n,)
        Vegetation type at each point.
    dens_center : ndarray of int, shape (n,)
        Density bin at each point.

    Notes
    -----
    This is a convenience function for workflows where only the center class is
    needed.

    It is the appropriate sampling method for categorical inputs used directly
    as class labels, for example canopy-height logic based on node type and node
    density bin. In that setting, averaging categories over a neighborhood would
    blur category boundaries and produce values with no categorical meaning.
    """
    xy = np.column_stack([np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)])

    type_center, _, _, _ = sample_raster_stencil_int(
        xy, veg_type_tif, stencil=(1, 1), strict=strict, nan_to=0
    )
    dens_center, _, _, _ = sample_raster_stencil_int(
        xy, veg_density_tif, stencil=(1, 1), strict=strict, nan_to=0
    )
    return type_center, dens_center


def sample_stencil(xs, ys, veg_type_tif, veg_density_tif, *, stencil=3, strict=False):
    """
    Sample vegetation type and density-bin stencils around point centers.

    Parameters
    ----------
    xs, ys : array-like of float
        Point coordinates in the raster coordinate system.
    veg_type_tif : str or path-like
        Vegetation type raster.
    veg_density_tif : str or path-like
        Vegetation density-bin raster.
    stencil : int or (int, int), optional
        Neighborhood size centered on the containing pixel.
    strict : bool, optional
        Controls how uncovered raster samples are handled.

    Returns
    -------
    type_stencil : ndarray of int, shape (n, sy*sx)
        Flattened type stencil for each point.
    dens_stencil : ndarray of int, shape (n, sy*sx)
        Flattened density-bin stencil for each point.

    Notes
    -----
    This is mainly a convenience function for testing and inspection. It is
    useful when verifying that a node sees the expected local categorical
    neighborhood.

    The flattening order is row-major. For a 3x3 stencil, the center cell is
    index 4.
    """
    xy = np.column_stack([np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)])

    _, type_stencil, _, _ = sample_raster_stencil_int(
        xy, veg_type_tif, stencil=stencil, strict=strict, nan_to=0
    )
    _, dens_stencil, _, _ = sample_raster_stencil_int(
        xy, veg_density_tif, stencil=stencil, strict=strict, nan_to=0
    )
    return type_stencil, dens_stencil

