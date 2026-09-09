SCHISM hotstart
===============

``schimpy.schism_hotstart`` creates a SCHISM ``hotstart.nc`` from a YAML
configuration. It supports both an initial condition assembled from observations
and formulas, and transfer of a prior hotstart onto the same or a changed mesh.

Run the YAML directly; a separate Python driver is not needed::

      create_hotstart hotstart.yaml

The same command is available through the grouped CLI::

      sch create_hotstart hotstart.yaml

Core configuration
------------------

The ``hotstart`` block defines the target grid and model clock, followed by one
initializer for every requested variable. For example::

      hotstart:
         date: 2021-10-05
         run_start: 2020-09-30
         time_step: 90
         hgrid_input_file: hgrid.gr3
         vgrid_input_file: vgrid.in.3d
         vgrid_version: "5.10"
         modules: [TEM, SAL]

         elevation:
            initializer:
               simple_trend:
                  value: max(0.97, -z-0.01)
         temperature:
            initializer:
               simple_trend:
                  value: 14.0
         salinity:
            initializer:
               simple_trend:
                  value: 0.0
         velocity_u:
            initializer:
               simple_trend:
                  value: 0.0
         velocity_v:
            initializer:
               simple_trend:
                  value: 0.0
         velocity_w:
            initializer:
               simple_trend:
                  value: 0.0

For a new initial condition, ``run_start: default`` sets the origin to ``date``.
For an ``ihot=2`` continuation, ``run_start`` is the original simulation origin
and ``date`` is the restart moment. Schimpy derives ``time``, ``iths`` and
``nsteps_from_cold`` from those values and ``time_step``.

Initializers
------------

The supported initializers are:

``simple_trend``
      A constant or an expression in target-node ``x``, ``y`` and ``z``. In this
      context ``z`` is depth, positive down.

``obs_points``
      Interpolation from station observations.

``extrude_casts``
      Interpolation and vertical extrusion from profile or cruise data.

``text_init``
      Values from a GR3-style ``.ic`` or ``.gr3`` file, or another supported
      text initializer input.

``hotstart_nc``
      Values from a prior hotstart, optionally transferred from a source grid and
      vertical grid.

``patch_init``
      Dispatch to different initializers by region. ``regions_filename`` may be a
      region shapefile, an ``.ic`` file, or a schimpy polygon YAML file.

``schout_nc`` is reserved but is not implemented.

When elevation uses ``hotstart_nc``, its initializer must set a non-negative
``max_blw_bed``. This is a lower bound on free-surface elevation for novel target
nodes. A value of ``0.01`` permits the initialized surface to sit at most one
centimetre below the target bed::

      elevation:
         initializer:
            hotstart_nc:
               data_source: source_hotstart.nc
               source_hgrid: source_hgrid.gr3
               source_vgrid: source_vgrid.in.3d
               source_vgrid_version: "5.10"
               max_blw_bed: 0.01

Wet and dry flags
-----------------

Wet/dry handling depends on the elevation source. At target nodes that coincide
with nodes in an elevation ``hotstart_nc``, schimpy retains the source ``idry``
flag. Nodes without a matched source flag, including nodes initialized by
``simple_trend`` or ``text_init``, are evaluated on the target grid. A node is
dry when

.. math::

    H = dp + eta \le h0.

Side and element flags are then derived from the completed target-node flags; a
side or element is dry if any of its nodes is dry.

Changed grids and inundation
----------------------------

For a changed-grid continuation, use ``patch_init`` for elevation. The unchanged
domain normally uses ``hotstart_nc`` while each new area receives an explicit
initializer. Tracers, velocities and turbulence variables usually transfer from
the prior hotstart.

``schimpy.inundate_island`` generates four mutually consistent inputs for a
gradually inundated restoration area:

* ``depth_enforce_inundate.yaml`` for the breach dredge;
* ``elev_inundate.yaml`` for ``elev.ic``;
* ``hydraulic_structures_inundate.yaml`` for the temporary structures; and
* ``inundate_regions.yaml`` for hotstart ``patch_init``.

Generate them with::

      sch inundate_island --config breaches.yaml --hgrid hgrid.gr3 --out-dir .

The generated regions contain ``domain`` first and one entry per restoration
area. In the hotstart YAML, set ``allow_overlap: true`` and preserve that order:
list ``domain`` first and restoration regions afterwards, because the last
matching configured region wins. A common pattern is ``hotstart_nc`` for
``domain`` and ``text_init`` from the generated ``elev.ic`` for each restoration
region.

Examples
--------

Applied Bay-Delta configurations are stored in the
`BayDeltaSCHISM hotstart examples
<https://github.com/CADWRDeltaModeling/BayDeltaSCHISM/tree/master/examples/hotstart>`_.
They demonstrate regional constants, cruise casts, observations, changed-grid
transfer, sediment, age and biology modules.

These examples are reference configurations, not self-contained test cases. The
shared target grids, vertical grids and source hotstarts are not distributed, and
some case directories retain legacy ``create_hotstart.py`` drivers or YAML that
predates current required keys. Use the CLI and current schimpy source for the
interface contract. The maintained generic inundation behavior is covered by
``tests/test_inundate_island.py`` and the inundation notebook in
``docsrc/notebooks/inundate_island.ipynb``.



