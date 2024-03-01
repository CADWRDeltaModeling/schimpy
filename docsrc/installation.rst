.. highlight:: shell

============
Installation
============


Stable release
--------------

We recommend conda for installing schimpy. There is a general description of conda on our 
`python configuration <https://cadwrdeltamodeling.github.io/html/python.html#python-configuration-and-usage>`_ page.

The most functional way to install schimpy is to use our complete environment for schism, schism_env.yml in the parent directory of 
`BayDeltaSCHISM <https://github.com/CADWRDeltaModeling/BayDeltaSCHISM>`_

In principle, you can do a one-off install of the just the schimpy library from the command line by running this command in your conda-enabled terminal:

.. code-block:: console

    $ conda install -c cadwr-dms schimpy



From sources
------------

The source for schimpy can be downloaded from the `Github repo <https://github.com/CADWRDeltaModeling/schimpy>`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone https://github.com/CADWRDeltaModeling/schimpy

Or, more rarely, download the `tarball`_:

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

This is sometimes called a developer install, since you can edit while your work and contribute your work backto a repository. As an aside, you can do it the conda way to get the prerequisites and follow up later with a developer install.