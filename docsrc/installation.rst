.. highlight:: shell

============
Installation
============


Stable release
--------------

To install schimpy, run this command in your terminal:

.. code-block:: console

    $ conda install -c cadwr-dms schimpy

This is the preferred method to install schimpy, as it will always install the most recent stable release.

If you don't have `conda`_ installed, this `Conda installation guide`_ can guide
you through the process.

.. _conda: https://docs.conda.io/en/latest/
.. _Conda installation guide: https://docs.conda.io/projects/conda/en/latest/user-guide/install/


From sources
------------

The sources for schimpy can be downloaded from the `Github repo`_.

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
    $ pip install -e .


.. _Github repo: https://github.com/CADWRDeltaModeling/schimpy
.. _tarball: https://github.com/CADWRDeltaModeling/schimpy/tarball/master
