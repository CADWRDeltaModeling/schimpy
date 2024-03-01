.. schimpy documentation master file, created by sphinx-quickstart.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _tools-home:

SCHIMPY Toolset
=================

The *schimpy* toolset is designed to help to pre- and post-process SCHISM
simulations more easily. In particular, the tools allow you
to specify things like boundaries and hydraulic structures using
geophysical coordinates rather than node/element numbers ... this
can greatly simplify the work involved when the mesh is revised.

Please direct questions, suggestions or bug reports to
Eli.Ateljevich@water.ca.gov or Kijin.Nam@water.ca.gov.

Usage
=====

To use schimpy tools in a project::

    import schimpy
	
Many schimpy utilities can be invoked as command line utilities and the most common of these is the `prepare_schism` command for the `preprocessor`_ ::

	$ prepare_schism main_bay_delta.yaml 

where `main_bay_delta.yaml` is the main yaml input file for your project.

Contents
---------

.. toctree::
   :maxdepth: 2
   
   Installation <installation>
   Spatial data in SCHISM <spatial>
   Pre-processing <preprocessing>
   Working with SMS <sms>
   Model time <model_time>
   Populating Elevation/Depth Data <depth>
   Creating hotstart <hotstart>
   Creating nudging <nudging>
   Post-processing <postprocessing>
   Evaluating Model Performance <performance>
   Using Python Pre-processing Script Programatically <direct_use>
   Modules <modules>
   History <history>
   License <license_and_disclaimer>
   Contributing <contributing>
   Authors <authors>
   Project README <readme>
   

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
