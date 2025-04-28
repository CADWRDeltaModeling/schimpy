# -*- coding: utf-8 -*-
"""Optimized raster_to_nodes function without rasterstats, using rasterio directly."""

import numpy as np
import rasterio
from rasterio.mask import mask
from shapely.geometry import Polygon, mapping
import itertools

def raster_to_nodes(mesh, nodes_sel, path_raster,
                    bins=None,
                    mapped_values=None,
                    mask_value=None,
                    fill_value=None,
                    band=1):
    """
    Applies raster values to mesh by calculating means over surrounding elements.

    Parameters
    ----------
    (same as original)

    Returns
    -------
    numpy.array
        means of raster values of the element balls around the nodes
    """

    # Step 1: Prepare binning if needed
    classify_raster = False
    if bins is not None:
        if mapped_values is None:
            raise ValueError("mapped_values must be provided if bins are used.")
        if len(mapped_values) != len(bins) + 1:
            raise ValueError("mapped_values must be one longer than bins.")
        classify_raster = True

    # Step 2: Precompute node balls (cache)
    elements_in_balls = {node_i: list(mesh.get_elems_i_from_node(node_i)) for node_i in nodes_sel}
    elements_in_polygon = sorted(set(itertools.chain.from_iterable(elements_in_balls.values())))

    # Step 3: Build element polygons
    element_polygons = {}
    for elem_i in elements_in_polygon:
        nodes_idx = mesh.elem(elem_i)
        coords = mesh.nodes[nodes_idx, :2]  # x,y only
        element_polygons[elem_i] = Polygon(coords)

    # Step 4: Open raster
    with rasterio.open(path_raster) as src:
        nodata = src.nodatavals[band-1] if src.nodatavals else None
        
        # Read full raster band if binning is needed
        if classify_raster:
            raster_array = src.read(band)
            if mask_value is not None:
                raster_array = np.where(raster_array == mask_value, fill_value, raster_array)
            digitized = np.digitize(raster_array, bins, right=False)
            classified_array = np.array(mapped_values)[digitized]
        else:
            classified_array = None  # Will read on demand

        # Step 5: Calculate mean per element
        sav_in_elements = {}
        for elem_i, polygon in element_polygons.items():
            geom = [mapping(polygon)]
            if classify_raster is False:
                out_image, out_transform = mask(src, geom, crop=True, indexes=band, nodata=nodata)
                data = out_image[0]
            else:
                # Clip manually from classified_array
                bounds = polygon.bounds  # (minx, miny, maxx, maxy)
                row_min, col_min = src.index(bounds[0], bounds[3])  # (xmin, ymax)
                row_max, col_max = src.index(bounds[2], bounds[1])  # (xmax, ymin)
                rows = slice(min(row_min, row_max), max(row_min, row_max)+1)
                cols = slice(min(col_min, col_max), max(col_min, col_max)+1)
                if rows.start == rows.stop or cols.start == cols.stop:
                    sav_in_elements[elem_i] = 0.0
                    continue
                window = classified_array[rows, cols]
                # Verify that window has nonzero size
                if window.shape[0] == 0 or window.shape[1] == 0:
                    sav_in_elements[elem_i] = 0.0
                    continue
                # Build temporary profile
                transform = src.window_transform(((rows.start, rows.stop), (cols.start, cols.stop)))
                with rasterio.io.MemoryFile() as memfile:
                    with memfile.open(
                        driver='GTiff',
                        height=window.shape[0],
                        width=window.shape[1],
                        count=1,
                        dtype=window.dtype,
                        transform=transform,
                        crs=src.crs,
                        nodata=-999
                    ) as dataset:
                        dataset.write(window, 1)
                        out_image, out_transform = mask(dataset, geom, crop=True, indexes=1, nodata=-999)
                        data = out_image[0]

            masked = np.ma.masked_array(data, mask=(data == nodata))
            mean_val = masked.mean() if masked.count() > 0 else 0.0
            sav_in_elements[elem_i] = mean_val

    # Step 6: Assemble node values
    elem_areas = mesh.areas()
    sav_at_nodes = np.empty((len(nodes_sel),), dtype=float)
    for i, node_i in enumerate(nodes_sel):
        ball = elements_in_balls[node_i]
        values = [sav_in_elements[e] for e in ball]
        weights = elem_areas[ball]
        sav_at_nodes[i] = np.average(values, weights=weights)

    return sav_at_nodes
