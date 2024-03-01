
Utilities for Working with SMS
==============================

The preprocessing toolkit has some utilities that can help when working with 
Aquaveo SMS, the main form of mesh generation we use until we are
able to distrubute a sufficient open source mesh generator.

Converting an SMS mesh (2dm) file to gr3
----------------------------------------

The script that converts SMS mesh to SCHISM gr3 format is available as a standalone script.

convert_mesh
^^^^^^^^^^^^

schimpy provides a single conversion utility between sms, shapefiles and SMS 2dm files

.. argparse::
    :module: schimpy.convert_mesh
    :func: create_arg_parser
    :prog: convert_mesh


Getting SMS Polygons for Preprocessor
-------------------------------------

SMS is able to create mesh coverages that are not mesh generation hints, including
spatial properties maps whose polygons have \"material properties\". If you save the
coverage as a polygon shapefile, the script can help convert the
resulting shapefile to the preprocessor yaml format. At the time of writing we do
not have a tutorial example of this process.

convert_polygons
^^^^^^^^^^^^^^^^

.. argparse::
    :module: schimpy.convert_polygons
    :func: create_arg_parser
    :prog: convert_polygons



Visualizing detailed bathymetry inside SMS
------------------------------------------

Earlier versions of SMS had bugs causing incorrect georefercing of DEMs
and very poor memory management. That is fixed but for efficiency reasons it continues to be impossible to efficiently assign elevations within the software for a large mesh. 
There are several ways of tackling this, which are covered in the Bay-Delta SCHISM
documentation. Here we have two scripts
that were designed to help with this problem: `clip_dems.py` and `stacked_dem_fill.py`. 


clip_dems
^^^^^^^^^

.. argparse::
    :module: schimpy.clip_dems
    :func: create_arg_parser
    :prog: clip_dems

The script clip_dems.py is used to create mini-DEMs on the fly that are small 
enough to be brought into the SMS environment as rasters without slowing it down too much. 
The same The :ref:`dem_list_file` is used here as elsewhere in the pre-processor.

The normal pattern of work is to zoom into the area of interst in SMS, save the
work area as a jpeg image and use it as the `image` argument to `clip_dems.py`.
One cut DEM will be produced for each DEM on the `dem_list` file that intersects the 
work area, using file names in order of priority (0 highest). The cut DEMs would then 
be imported and visualized in SMS, ideally using contours as the 
viewing choice and a carefully crafted color map designed to focus on the 
elevations that are locally most important.

stacked_dem_fill
^^^^^^^^^^^^^^^^

The script stacked_dem_fill.py is the one we use to populate a mesh with elevations.


.. argparse::
    :module: schimpy.stacked_dem_fill
    :func: create_arg_parser
    :prog: stacked_dem_fill

