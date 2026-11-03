.. _bctides:

Boundary conditions (bctides)
=============================

The ``bctide`` tool writes the SCHISM ``bctides.in`` file from a YAML
description of the open boundaries. The YAML names each boundary and states,
variable by variable, where its forcing comes from. The tool maps each choice
onto the integer type codes SCHISM expects, writes the accompanying parameters
(constants, tidal amplitudes and phases, nudging and relaxation factors) and
logs the resulting assignment for every boundary.

The keywords and checks below were verified against ``src/Hydro/schism_init.F90``
and ``src/Hydro/schism_step.F90`` in SCHISM.

.. contents:: On this page
   :local:
   :depth: 2

Running the tool
----------------

.. click:: schimpy.bctide:bctide_cli
   :prog: bctide
   :nested: full

The preprocessor (``prepare_schism``) calls the same machinery when its main
YAML contains a ``bctides`` section.

Main YAML
---------

The main YAML names the mesh and lists one entry per ``bctides`` file. Entries
are usually kept in separate files and pulled in with ``include``:

.. code-block:: yaml

   mesh:
     mesh_inputfile: ./hgrid.gr3
     open_boundaries: !include open_boundary.yaml   # optional

   bctides:
     include:
       - bctides.in.2d.yaml
       - bctides.in.3d.yaml

When ``mesh.open_boundaries`` is given, its ``linestrings`` replace the open
boundaries stored in ``hgrid.gr3``. The linestrings must have the same number
and names as the ``open_boundaries`` of each bctides entry.

Boundary names in the YAML must match the names in the ``hgrid.gr3`` boundary
comments (``! Open boundary "name"``). Entries are reordered to hgrid order, so
the YAML order does not matter.

bctides entry
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``name``
     - Output file name, e.g. ``bctides.in.3d``.
   * - ``date``
     - Start date written on the first line. Optional; defaults to 2020-01-01.
   * - ``mode``
     - Required. ``barotropic`` or ``baroclinic``.
   * - ``temperature_salinity_computation``
     - Required when ``mode`` is ``barotropic``. If ``true``, inflow and
       bidirectional boundaries must specify temperature and salinity; if
       ``false``, they must not.
   * - ``modules``
     - Tracer modules, each ``{name, num_tracers}``. Names are ``GEN``,
       ``AGE``, ``SED``, ``ECO``, ``ICM``, ``COSINE``, ``FECO``, ``TIMOR``,
       ``FABM`` and ``DVD``. They are written in that (SCHISM) order regardless
       of YAML order.
   * - ``earth_tides``
     - Tidal potential: ``cutoff_depth`` and ``constituents``, each with
       ``amplitude``, ``node_factor``, ``angular_frequency`` and
       ``equilibrium_argument``. Omit for none.
   * - ``boundary_forcing_tides``
     - Tidal forcing frequencies: ``constituents``, each with
       ``angular_frequency``, ``node_factor`` and ``equilibrium_argument``.
       Required when any boundary uses a tidal source.
   * - ``open_boundaries``
     - List of open boundaries, described next.

Open boundaries
---------------

Each open boundary has a ``name``, a required ``flow_direction`` and a list of
``variables``. Each item in ``variables`` is a single variable:

.. code-block:: yaml

   open_boundaries:
     - name: ocean
       flow_direction: bidirection
       variables:
         - elevation:
             source: th_nc
         - velocity:
             source: th_nc
         - temperature:
             source: initial_profile
             nudge: 0.01
         - salinity:
             source: initial_profile
             nudge: 0.01
     - name: sjr
       flow_direction: inflow
       variables:
         - velocity:
             source: th_text
         - temperature:
             source: th_text
             nudge: 1.0
         - salinity:
             source: 0.074318
             nudge: 1.0

``flow_direction`` is ``inflow``, ``outflow`` or ``bidirection``:

* ``outflow`` boundaries may not specify temperature, salinity or tracers.
* ``inflow`` and ``bidirection`` boundaries in ``baroclinic`` mode must
  specify temperature and salinity, and tracers when ``modules`` is present.
  In ``barotropic`` mode the requirement follows
  ``temperature_salinity_computation``; tracers are not allowed.

A variable that is left out gets SCHISM type code 0 (not specified).

Source keywords
---------------

The ``source`` of a variable is either a number, meaning a constant (code 2),
or one of the keywords below. Keywords contain no spaces and there are no
alternative spellings. ``th_text`` means the ASCII time history file and
``th_nc`` the netCDF time history file; the file name is fixed by SCHISM for
the variable, as listed.

Elevation (``iettype``)
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 8 72

   * - ``source``
     - Code
     - Meaning
   * - ``th_text``
     - 1
     - Time history from ``elev.th``, uniform along the boundary.
   * - number
     - 2
     - Constant elevation.
   * - ``tidal``
     - 3
     - Tidal harmonics given in ``constituents``.
   * - ``th_nc``
     - 4
     - Time and space varying history from ``elev2D.th.nc``.
   * - ``tidal_th_nc``
     - 5
     - ``elev2D.th.nc`` plus the tidal harmonics in ``constituents``.

Velocity (``ifltype``)
^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 24 8 68

   * - ``source``
     - Code
     - Meaning
   * - ``th_text``
     - 1
     - Discharge history from ``flux.th`` (negative is inflow).
   * - number
     - 2
     - Constant discharge (negative is inflow).
   * - ``tidal``
     - 3
     - Tidal u, v harmonics given in ``constituents``.
   * - ``th_nc``
     - 4
     - Velocity history from ``uv3D.th.nc``.
   * - ``tidal_th_nc``
     - 5
     - ``uv3D.th.nc`` plus the tidal harmonics in ``constituents``.
   * - ``th_nc_relaxed``
     - -4
     - Model velocity relaxed toward ``uv3D.th.nc`` with separate factors for
       incoming and outgoing flow. Requires ``inflow_relax``,
       ``outflow_relax`` and an elevation boundary.
   * - ``tidal_th_nc_relaxed``
     - -5
     - As ``th_nc_relaxed``, with the target including the tidal harmonics in
       ``constituents``.
   * - ``flather_1``
     - -1
     - Flather radiation condition. Requires ``eta_m0`` and ``vnormal_nodes``
       and no elevation boundary.
   * - ``discharge_relation``
     - -2
     - Outgoing discharge from the polynomial
       :math:`Q = \sum_{j=1}^{4} c_j \eta^j`. Requires ``coefficients`` and no
       elevation boundary.

SCHISM's own manual does not list codes -5 and -2; they are read and applied by
the SCHISM source.

Temperature, salinity and tracers (``itetype``, ``isatype``, ``itrtype``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 8 72

   * - ``source``
     - Code
     - Meaning
   * - ``th_text``
     - 1
     - Time history from ``TEM_1.th``, ``SAL_1.th`` or ``[MOD]_[ID].th``.
   * - number or list
     - 2
     - Constant value. Tracers take a list with one value per tracer in the
       module (a single number is accepted when ``num_tracers`` is 1).
   * - ``initial_profile``
     - 3
     - Relax to the initial condition.
   * - ``th_nc``
     - 4
     - Time and space varying history from ``TEM_3D.th.nc``, ``SAL_3D.th.nc``
       or ``[MOD]_3D.th.nc``.

Every temperature, salinity and tracer entry requires ``nudge``, a factor in
[0, 1] with 1 the strongest nudging. SCHISM applies these values on inflow.

Tracers are listed per module:

.. code-block:: yaml

   - tracers:
       - module: SED
         source: [0.007, 0.008, 0.0]
         nudge: 0.01
       - module: GEN
         source: th_text
         nudge: 1.0

A module that is left out of a boundary's list gets code 0 on that boundary.

Additional keys by velocity type
--------------------------------

``tidal``, ``tidal_th_nc``, ``tidal_th_nc_relaxed`` (and elevation ``tidal``, ``tidal_th_nc``)
   ``constituents``: one entry per ``boundary_forcing_tides`` constituent, in
   the same order. Elevation entries give ``amplitude`` and ``phase``;
   velocity entries give ``u_amplitude``, ``u_phase``, ``v_amplitude`` and
   ``v_phase``. Each value is a scalar (uniform) or a list with one value per
   boundary node.

   .. code-block:: yaml

      - elevation:
          source: tidal
          constituents:
            - O1:
                - amplitude: 0.26
                - phase: 123.0

``th_nc_relaxed``, ``tidal_th_nc_relaxed``
   ``inflow_relax`` and ``outflow_relax``: relaxation factors in [0, 1] for
   incoming and outgoing flow (1 imposes ``uv3D.th.nc`` fully).

   .. code-block:: yaml

      - elevation:
          source: th_nc
      - velocity:
          source: th_nc_relaxed
          inflow_relax: 0.8
          outflow_relax: 0.1

``flather_1``
   ``eta_m0``: mean elevation, one value per boundary node. ``vnormal_nodes``:
   mean normal velocity, one list per boundary node with a value per vertical
   level.

``discharge_relation``
   ``coefficients``: the four polynomial coefficients :math:`c_1 \ldots c_4`.

Checks
------

In addition to the keyword and ``flow_direction`` rules above, ``bctide``
raises ``ValueError`` before writing anything SCHISM would reject at start-up:

* ``flather_1`` and ``discharge_relation`` with an elevation boundary;
* ``th_nc_relaxed`` and ``tidal_th_nc_relaxed`` without an elevation boundary;
* ``nudge``, ``inflow_relax`` or ``outflow_relax`` missing or outside [0, 1];
* ``discharge_relation`` without exactly four ``coefficients``;
* ``flather_1`` per-node lists whose length differs from the boundary;
* tidal sources without ``constituents``;
* a constant tracer list whose length differs from ``num_tracers``;
* a tracer module not declared in ``modules``, or listed twice;
* a ``variables`` item that is not a single variable (often an indentation
  error that leaves a variable empty).

Logging
-------

Each boundary's assignment is logged at INFO level, for example::

   Boundary 'ocean' (bidirection): elevation=th_nc(4), velocity=th_nc(4),
   temperature=initial_profile(3) nudge=0.01, salinity=initial_profile(3) nudge=0.01

Use ``--logdir`` to keep the log alongside the run.
