
Assigning Bathymetry
====================

Converting a 2dm file to gr3 or between gr3 or 2dm and shapefile


.. argparse::
    :module: schimpy.convert_mesh
    :func: create_arg_parser
    :prog: convert_mesh

Assigning depths from prioritized DEMs
--------------------------------------

These tools are used by the preprocessor but also work with functions that can
be used generically and consistently with SCHISM.

stacked_dem_fill
^^^^^^^^^^^^^^^^

.. argparse::
    :module: schimpy.stacked_dem_fill
    :func: create_arg_parser
    :prog: stacked_dem_fill


Smoother for complex topography
-------------------------------

Often it is necessary to incorporate inudated marshy areas where elevations
are poorly sampled and contours are tortuous. The script contour_smooth.py
uses minmax curvature flow (Malladi and Sethian) to impose a minimum length scale
of change for contours, essentially unraveling the features that are most contorted.

.. argparse::
    :module: schimpy.contour_smooth
    :func: create_arg_parser
    :prog: contour_smooth

Optimizing depths for volume
----------------------------

A mesh draped over noisy bathymetry data does not necessarily represent important moments 
such as volumes and vertical face areas smoothly and realistically. 
To better represent these facets of the geometry, 
we compare the volumetric quantities that come from SCHISM shape functions 
(which are much like a TIN) to a higher order quadrature using the DEM 
with bilinear interpolation. The quadrature is more accurate, and also 
incorporates more sample points.

.. argparse::
    :module: schimpy.grid_opt
    :func: create_arg_parser
    :prog: grid_opt
