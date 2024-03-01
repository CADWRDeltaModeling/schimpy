
SCHISM Spatial Data
===================

Formats
-------

SCHISM requires numerous distributed spatial inputs, including node and element-centered properties, hydraulic structures specifications, output locations and so on. 

Under the native SCHISM input system, many spatial inputs are specified based on node or cell indices rather than geophysical coordinates. The user must reinterpret geophysical information. Whenever changes are made in the model mesh, the corresponding node and cell indices changes. Doing this input conversion of spatial input information manually from geophysical coordinates to node and cell indices in a large mesh is not only very cumbersome but error-prone. Our preprocessing system, as well as VIMS pyschism, address these issues. 

This is the main job addressed by the schimpy Python library, but there are some input and output processing tools as well. For instance, the model_time module converts inputs from geophysical timestamps to elapsed time. The nudging and hotstart modules facilitate the creation of files for data assimilation (nudging or Newtonian relaxation) and initial conditions respectively. 




