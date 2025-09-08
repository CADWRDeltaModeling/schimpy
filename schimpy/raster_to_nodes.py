import numpy as np
import itertools
import rasterio
from rasterio.features import rasterize
from rasterio.windows import from_bounds
from shapely.geometry import Polygon, box

# ---------- helpers (minimal) ----------

def _union_missing_mask(ma, extra_mask_values):
    """Union raster mask with explicit sentinels (e.g., 0.0, -9999)."""
    if extra_mask_values is None:
        return ma
    if not isinstance(extra_mask_values, (list, tuple, np.ndarray)):
        extra_mask_values = [extra_mask_values]
    data = ma.data
    add_mask = np.zeros(data.shape, dtype=bool)
    for mv in extra_mask_values:
        add_mask |= (data == mv)
    return np.ma.array(data, mask=(np.asarray(ma.mask, bool) | add_mask), copy=False)

def _elements_bounds_in_raster(elements_polys, raster_bounds):
    """BBox of union(elements ∩ raster) or None if empty."""
    rb = box(*raster_bounds)
    any_hit = False
    xmin, ymin, xmax, ymax = +np.inf, +np.inf, -np.inf, -np.inf
    for poly in elements_polys:
        if poly.is_empty: 
            continue
        inter = poly.intersection(rb)
        if inter.is_empty: 
            continue
        bxmin, bymin, bxmax, bymax = inter.bounds
        xmin, ymin = min(xmin, bxmin), min(ymin, bymin)
        xmax, ymax = max(xmax, bxmax), max(ymax, bymax)
        any_hit = True
    return (xmin, ymin, xmax, ymax) if any_hit else None

# ---------- main (rasterize-only) ----------

def raster_to_nodes(
    mesh,
    nodes_sel,
    path_raster,
    bins=None,
    mapped_values=None,
    mask_value=None,
    fill_value=None,
    band=1,
    missing_policy="exclude",
    rasterize_all_touched=True,
):
    """
    Map raster values onto mesh nodes by averaging values over surrounding elements.

    This function samples a raster onto an unstructured mesh. Each element polygon
    is intersected with the raster grid, values are averaged per element, and then
    element values are aggregated to nodes using an area-weighted mean.

    Parameters
    ----------
    mesh : object
        Mesh object with the following interface:
          - ``mesh.nodes`` : array of node coordinates (N, >=2), x and y in the first two columns.
          - ``mesh.elem(i)`` : return node indices for element ``i``.
          - ``mesh.get_elems_i_from_node(node_i)`` : return element indices connected to node ``node_i``.
          - ``mesh.areas()`` : array of element areas indexed by element id.
    nodes_sel : sequence of int
        Node indices at which values will be computed.
    path_raster : str
        Path to a raster file (e.g., GeoTIFF) readable by rasterio.
    bins : sequence of float, optional
        Bin edges used to classify raster values. If provided, the raster values
        are digitized into bins before aggregation. Length must be ``M``.
    mapped_values : sequence of float, optional
        Values assigned to each bin index. Must have length ``M+1`` if ``bins`` is provided.
        Typically used to map classes into the interval [0, 1].
    mask_value : scalar or sequence of scalars, optional
        Additional raster values to treat as missing (in addition to raster
        NoData). Example: ``[0.0, -9999]``.
    fill_value : float, optional
        Default value for elements with no valid raster pixels in the continuous
        (non-classified) case. If omitted, defaults to 0.0. Ignored when
        classification is active.
    band : int, default=1
        1-based band index in the raster to sample.
    missing_policy : {"exclude", "as_zero"}, default="exclude"
        Policy for handling missing pixels in the classified case:
        
        - ``"exclude"`` : compute class averages using only valid pixels. Elements
          with no valid pixels receive the default class value (0.0).
        - ``"as_zero"`` : missing pixels count as class 0 by including them in the
          denominator. This dilutes averages toward 0 when coverage is sparse.
          
        Has no effect in the continuous (non-classified) case.
    rasterize_all_touched : bool, default=True
        If True, count all raster cells touched by element polygons. If False,
        count only cells whose centers fall within the polygons.

    Returns
    -------
    ndarray of float
        Array of values at each node in ``nodes_sel``. Length matches ``len(nodes_sel)``.

    Notes
    -----
    - Continuous mode (no ``bins``): averages raster values directly.
    - Classified mode (with ``bins`` and ``mapped_values``): raster values are
      digitized into bins, mapped to user-provided class values, then averaged.
    - All aggregation is area-weighted: element values are computed from raster
      pixels, and node values are weighted by connected element areas.
    """

    classify = bins is not None
    if classify:
        if mapped_values is None:
            raise ValueError("mapped_values must be provided if bins are used.")
        if len(mapped_values) != len(bins) + 1:
            raise ValueError("mapped_values must be one longer than bins.")
        if missing_policy not in ("exclude", "as_zero"):
            raise ValueError("missing_policy must be 'exclude' or 'as_zero'.")

    # Node -> elements and unique element ids
    elements_in_balls = {n: list(mesh.get_elems_i_from_node(n)) for n in nodes_sel}
    elements_all = sorted(set(itertools.chain.from_iterable(elements_in_balls.values())))

    # Build element polygons (XY only)
    element_polygons = {}
    for elem_i in elements_all:
        nodes_idx = mesh.elem(elem_i)
        coords = mesh.nodes[nodes_idx, :2]
        element_polygons[elem_i] = Polygon(coords)

    # Defaults per path
    elem_default_cont  = (fill_value if fill_value is not None else 0.0)
    elem_default_class = 0.0

    with rasterio.open(path_raster) as src:
        # Minimal window that covers mesh ∩ raster
        bbox = _elements_bounds_in_raster(list(element_polygons.values()), src.bounds)
        # If nothing overlaps: return defaults
        if bbox is None:
            default_val = elem_default_class if classify else elem_default_cont
            return np.full(len(nodes_sel), default_val, dtype=float)

        win = from_bounds(*bbox, transform=src.transform)
        full_win = from_bounds(*src.bounds, transform=src.transform)
        win = win.intersection(full_win)

        # Read data as masked array over that window
        data = src.read(band, window=win, masked=True)  # (H,W) masked
        if data.ndim == 3:  # (1,H,W) -> (H,W) if env returns 3D
            data = data[0]
        data = _union_missing_mask(data, mask_value)

        # Rasterize element IDs on same window grid
        window_transform = rasterio.windows.transform(win, src.transform)
        shapes = [(poly, int(eid)) for eid, poly in element_polygons.items()]
        labels = rasterize(
            shapes,
            out_shape=data.shape,
            transform=window_transform,
            fill=0,           # background label
            dtype="int32",
            all_touched=rasterize_all_touched,
        )

        # Prepare result per-element
        elem_values = {eid: (elem_default_class if classify else elem_default_cont)
                       for eid in element_polygons.keys()}

        # VALID pixels only
        valid = ~data.mask
        labs_valid = labels[valid]
        vals_valid = data.data[valid]
        sel_valid = labs_valid != 0

        if classify:
            # If there are valid class pixels, aggregate them
            if sel_valid.any():
                bins_arr = np.asarray(bins)
                mapped = np.asarray(mapped_values)
                idx = np.digitize(vals_valid[sel_valid], bins_arr, right=False)  # 0..len(bins)
                class_vals = mapped[idx]

                max_id = int(labels.max()) if labels.size else 0
                sums_valid = np.bincount(labs_valid[sel_valid], weights=class_vals, minlength=max_id + 1)
                cnts_valid = np.bincount(labs_valid[sel_valid], minlength=max_id + 1)
            else:
                max_id = int(labels.max()) if labels.size else 0
                sums_valid = np.zeros(max_id + 1, dtype=float)
                cnts_valid = np.zeros(max_id + 1, dtype=float)

            if missing_policy == "as_zero":
                # Denominator: all pixels of element (valid+masked), excluding background
                labs_all = labels.ravel()
                sel_all = labs_all != 0
                cnts_total = np.bincount(labs_all[sel_all], minlength=max_id + 1)
                present = np.nonzero(cnts_total)[0]
                for eid in present:
                    denom = cnts_total[eid]
                    num = sums_valid[eid]  # masked contribute 0
                    elem_values[eid] = float(num / denom) if denom > 0 else elem_default_class
            else:
                # "exclude": denominator is count of valid only
                present = np.nonzero(cnts_valid)[0]
                for eid in present:
                    denom = cnts_valid[eid]
                    elem_values[eid] = float(sums_valid[eid] / denom) if denom > 0 else elem_default_class

        else:
            # Continuous: mean over valid pixels per element
            if sel_valid.any():
                max_id = int(labels.max()) if labels.size else 0
                sums = np.bincount(labs_valid[sel_valid], weights=vals_valid[sel_valid], minlength=max_id + 1)
                cnts = np.bincount(labs_valid[sel_valid], minlength=max_id + 1)
                present = np.nonzero(cnts)[0]
                for eid in present:
                    elem_values[eid] = float(sums[eid] / cnts[eid]) if cnts[eid] > 0 else elem_default_cont
            # else: keep defaults (no valid pixels anywhere)

    # Area-weighted average from elements to nodes
    elem_areas = mesh.areas()
    out_vals = np.empty((len(nodes_sel),), dtype=float)
    for i, node_i in enumerate(nodes_sel):
        ball = elements_in_balls[node_i]
        if not ball:
            out_vals[i] = 0.0
            continue
        values = np.array([elem_values[e] for e in ball], dtype=float)
        weights = np.array(elem_areas[ball], dtype=float)
        bad = ~np.isfinite(values)
        if bad.any():
            weights = weights.copy(); weights[bad] = 0.0
            values = np.where(bad, 0.0, values)
        out_vals[i] = np.average(values, weights=weights) if weights.sum() > 0 else 0.0

    return out_vals
