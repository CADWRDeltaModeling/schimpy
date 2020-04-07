conda create -y -c defaults -n dev_schimpy python=3.7*
conda install -y -c defaults -n dev_schimpy rtree numpy gdal xarray matplotlib shapely pyyaml scipy pyproj geopandas
# Below this for dev environment, etc
conda install -y -c defaults -n dev_schimpy pytest pytest-runner
# Now for document generation
conda install -y -c defaults -c conda-forge -n dev_schimpy sphinx nbsphinx sphinx-argparse numpydoc
