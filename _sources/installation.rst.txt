.. highlight:: shell

============
Installation
============


Stable release
--------------

We recommend conda for installing schimpy. The simplest way to install 
schimpy is to run this command in your terminal:

.. code-block:: console

    $ conda install -c cadwr-dms schimpy

For greater flexibility and a set of tools that are commonly used with SCHISM, however, 
you may want to install a more complete schism environment using the 
schism_environment.yml file that is in the parent directory of 
`BayDeltaSCHISM <https://github.com/CADWRDeltaModeling/BayDeltaSCHISM>`_

From sources
------------

The sources for schimpy can be downloaded from the `Github repo <https://github.com/CADWRDeltaModeling/schimpy>`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone https://github.com/CADWRDeltaModeling/schimpy

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/CADWRDeltaModeling/schimpy/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    # make sure all dependencies from setup_devenv.bat are installed by conda
    $ setup_devenv.bat # or copy and paste the commands into shell
    $ conda activate dev_schimpy
    $ pip install --no-deps -e .


.. _Github repo: https://github.com/CADWRDeltaModeling/schimpy
.. _tarball: https://github.com/CADWRDeltaModeling/schimpy/tarball/master
