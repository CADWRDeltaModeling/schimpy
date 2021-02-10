# -*- coding: utf-8 -*-
""" Functions to process raster data for schimpy pre-processing.

    This add rasterstats, boto3 dependency to calculate raster stastistics such as
    means.
"""
import os
import tempfile
import itertools
from copy import deepcopy
import numpy as np
import gdal
from shapely.geometry import Polygon
import rasterstats

__all__ = ['raster_to_nodes', ]


def raster_to_nodes(mesh, nodes_i, path_raster,
                    bins=None,
                    mapped_values=None,
                    mask_value=None,
                    fill_value=None,
                    band=1):
    """ Calculate the means of raster values of the element balls
        around the nodes.

        When an element is completely outside of the raster data, zero will be
        assigned for the it.
        When an element is partly covered by the raster data, an masked average
        is calculated.

        If the bins are provided, the raster data will be binned or mapped to
        the mapped_values first and processed. The left end of the bins is
        included but the right end is not for binning.
        The length of the mapped values should be larger by one than
        that of the bins. The binning uses numpy.digitize to classification.

        If the mask_value is given, the raster data with the given mask value
        will be swapped with the fill_value before binning the data. This
        masking will work only when the bins are provided.

        An example of this function in the pre-processor yaml is below:

        gr3:
          sav_N.gr3:
            default: 0.
            polygons:
              - attribute: raster_to_nodes.raster_to_nodes(mesh, nodes_i, 'NDVI.tif', bins=[-998., 0.0001, 0.3, 0.6, 1.0], mapped_values=[-999., 0., 40., 70., 100., -999.], maske_value=-10., fill_value=0.1)
                imports: schimpy.raster_to_nodes
                type: none
                vertices:
                ...

        The variables "mesh" and "node_i" in the example above are
        'magic words' that the pre-processor understands.
        The exmaple bins (or classifies) the raster values from 0.0001 to 0.3
        to 0. and so forth. The outside of the bins are binned to -999.
        Also raster points with the value of -10. will be switched to 0.2 first
        before binning, so they will be mapped to 40.

        Parameters
        ----------
        mesh: schimpy.SchismMesh
            The base mesh to work on
        nodes_i: array-like
            The array or list of the nodes to process
        path_raster: string-like
            File fath to the raster data
        bins: array-like, optional
            Array of bins. It has to be 1-dimensional and monotonic
        mapped_values: array-like, optional
           The values of the classes. It should be bigger by one than the bins.
        band: int, optional
            The band index to use from the raster data.
            The default value is 1.

        Returns
        -------
        numpy.array
            means of raster values of the element balls around the nodes
    """
    if bins is not None:
        if mapped_values is None:
            raise ValueError(
                "The mapped values need")
        if len(mapped_values) != len(bins) + 1:
            raise ValueError(
                "The number of mapped_values should be smaller by one than the number of the bins.")

        # Read the raster data
        raster = gdal.Open(path_raster)
        band = raster.GetRasterBand(band)
        raster_array = band.ReadAsArray()

        # TODO:This converts the whole raster and not very efficient.
        if mask_value is not None:
            if fill_value is None:
                raise ValueError(
                    "The arugment fill_value must be given when mask_value is set.")
            raster_array_masked = deepcopy(raster_array)
            raster_array_masked[raster_array == mask_value] = fill_value
        else:
            raster_array_masked = raster_array
        digitized = np.digitize(raster_array_masked,
                                bins, right=False).reshape((-1,))
        classified_array = np.array(mapped_values)[
            digitized].reshape(raster_array.shape)

        # Save the classifed array
        driver = gdal.GetDriverByName("GTiff")
        path_temp = os.path.join(tempfile.mkdtemp(), 'temp.tif')
        outdata = driver.CreateCopy(path_temp, raster)
        outdata.GetRasterBand(1).WriteArray(classified_array)
        outdata.GetRasterBand(1).SetNoDataValue(-999)
        outdata.FlushCache()  # saves to disk!

    # Get the zonal data
    elements_in_balls = [list(mesh.get_elems_i_from_node(node_i))
                         for node_i in nodes_i]
    elements_in_polygon = list(
        set(list(itertools.chain.from_iterable(elements_in_balls))))
    element_polygons = [Polygon(mesh.nodes[mesh.elem(elem_i), :2])
                        for elem_i in elements_in_polygon]
    path_data = path_raster if bins is None else path_temp
    zs = rasterstats.zonal_stats(element_polygons, path_data, stats='mean')
    sav_in_elements = dict(zip(elements_in_polygon, [s['mean']
                                                     if s['mean'] is not None else 0.
                                                     for s in zs]))
    elem_areas = mesh.areas()
    sav_at_nodes = np.empty((len(nodes_i),),dtype=float)
    for i, node_i in enumerate(nodes_i):
        ball = list(mesh.get_elems_i_from_node(node_i))
        sav_at_nodes[i] = np.average(
            [sav_in_elements[e_i] for e_i in ball], weights=elem_areas[ball])

    # Remove the temporary file if exists
    if bins is not None:
        os.remove(path_temp)
    return sav_at_nodes
