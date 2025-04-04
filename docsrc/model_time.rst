Model time conversion
=====================

The SCHISM ``*.th`` file has two flavors, binary and text. Both are multicolumn formats,
with elapsed time in seconds in the left column since the start of the simulation
and the other columns representing values. 
There is no datum within the file to link elapsed time to calendar or geophysical time, 
and it is hard to search for events, re-use the ``*.th`` file with a different 
start date or make sense of orphaned files.

The model_time script utility performs conversions between geophysical
and elapsed times on flat files. This is a pretty easy task on simple files using
just pandas, but this utility preserves formatting and floating precision of the source, 
copies over comments and can clip the file to a new start date. 

model_time.py
-------------

.. argparse::
    :module: schimpy.model_time
    :func: create_arg_parser
    :prog: model_time
