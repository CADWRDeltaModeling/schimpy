.. _hydraulic_structures:

Hydraulic Structures
====================

SCHISM represents a hydraulic structure as two paired internal boundaries. The
model computes one total discharge from water levels at reference nodes, applies
that discharge as an outflow on one boundary and an equal inflow on the other,
and transfers tracers conservatively between the two sides. This mechanism
supports prescribed transfers, weirs, culverts, rectangular orifices, radial
gates, and compound weir-culvert structures.

This chapter describes both inputs used in a SCHISM workflow:

* a geographic YAML description consumed by the schimpy preprocessor; and
* the native ``hydraulics.in`` file consumed by SCHISM.

The equations and file formats below have been checked against
``src/Core/hydraulic_structures.F90`` in SCHISM. They describe the implemented
behavior rather than a general-purpose hydraulic design model.

.. contents:: On this page
   :local:
   :depth: 2

Mesh representation
-------------------

Each structure consists of corresponding node pairs on opposite sides of a
barrier. The first member of every pair belongs to the nominal upstream side
and the second to the nominal downstream side. ``upnode`` and ``downnode`` are
reference nodes at which SCHISM samples water-surface elevation for the flow
calculation. Actual flow may reverse; the nominal orientation controls the sign
and selects the directional operating coefficient.

.. figure:: documents/image/struct.png
   :alt: Paired node strings and reference nodes defining a hydraulic structure
   :align: center

   A hydraulic structure imposed on the horizontal grid.

The paired node strings must span the channel and terminate on land boundaries.
A gap around either end bypasses the structure. At least two node pairs are
required. Reference nodes must remain wet, as is required for other SCHISM
boundary nodes.

After SCHISM calculates the total structure flow, it distributes the flow as a
uniform flux over the paired boundary. The structure is not a geometric solid
inserted into the mesh. When a time series sets ``install`` to zero, SCHISM
removes the paired-boundary treatment and restores the ordinary equations of
motion across that part of the mesh.

Flow convention and relaxation
------------------------------

Let :math:`z_u` and :math:`z_d` be water-surface elevations at the nominal
upstream and downstream reference nodes, respectively. Positive discharge is
from the nominal upstream side to the nominal downstream side. SCHISM chooses
``op_down`` for positive flow and ``op_up`` for reverse flow. An operating
coefficient of zero closes the structure in that direction.

The hydraulic equation produces a target discharge :math:`Q_s`. SCHISM relaxes
the imposed boundary flow toward that target:

.. math::

   Q\left(t+\frac{\Delta t}{2}\right)
   = (1-\chi)Q(t) + \chi Q_s(t),

where :math:`\chi` is the global ``nudging`` value in ``hydraulics.in``. SCHISM
requires :math:`0 \leq \chi < 1`. Smaller values damp abrupt changes but also
slow the response to changing water levels and controls.

For the equations below, define

.. math::

   z_h &= \max(z_u,z_d), \\
   z_l &= \min(z_u,z_d), \\
   \Delta z &= z_h-z_l, \\
   C_{op} &=
   \begin{cases}
     C_{op,down}, & z_u > z_d,\\
     C_{op,up}, & z_u \leq z_d.
   \end{cases}

The configured ``coefficient``, :math:`C_f`, is a nonnegative physical flow
coefficient. The sign of :math:`z_u-z_d` supplies the flow direction. All
geometric and coefficient inputs describe one unit; SCHISM multiplies the
result by ``n_duplicates``.

Supported structure types
-------------------------

``transfer``
^^^^^^^^^^^^

A transfer prescribes a signed flow :math:`Q_s` directly. Positive flow follows
the nominal upstream-to-downstream orientation. The hydraulic calculation does
not use water levels, operating coefficients, or ``n_duplicates``. Tracer mass
is transferred from the outflow side to the inflow side.

``weir``
^^^^^^^^

A weir can be dry, free flowing, or submerged. Its ``elevation`` is the crest
(or invert) elevation relative to the model datum and ``width`` is the width of
one unit.

.. figure:: documents/image/weir.png
   :alt: Free-flowing and submerged weir cases
   :align: center

   Free-flowing (a) and submerged (b) weir flow.

With head above the crest

.. math::

   H = z_h-z_{inv},

flow is zero when :math:`H \leq 0`. Otherwise the free-flow discharge is

.. math::

   Q_f = \operatorname{sgn}(z_u-z_d)\,C_{op}C_f A\sqrt{2gH},
   \qquad A=WH.

When the lower water surface is above the crest, SCHISM applies the Villemonte
submergence correction

.. math::

   S &= \frac{z_l-z_{inv}}{H}, \\
   Q_s &= Q_f\left(1-S^{1.5}\right)^{0.385}.

The implementation uses water-surface head only; it does not add a velocity
head at the reference node.

``orifice``
^^^^^^^^^^^

An orifice is a rectangular opening with invert ``elevation``, ``width``, and
``height``. The wetted opening depth and area are

.. math::

   d = \min\left(H_{gate},z_h-z_{inv}\right), \qquad A=Wd.

Flow is zero for :math:`d \leq 0`; otherwise SCHISM uses

.. math::

   Q_s = \operatorname{sgn}(z_u-z_d)\,C_{op}C_f A
         \sqrt{2g\min(z_h-z_{inv},\Delta z)}.

Thus an opening with tailwater below its invert uses upstream head above the
invert, while a submerged opening uses the water-surface difference. No
additional partially submerged correction is applied.

``culvert``
^^^^^^^^^^^

A culvert is modeled as a circular orifice. ``elevation`` is the pipe invert and
``radius`` is the pipe radius (stored as ``width`` in ``hydraulics.in``). For
wetted depth :math:`d` and radius :math:`r`, SCHISM uses

.. math::

   A(d) =
   \begin{cases}
     r^2\theta-r(r-d)\sin\theta,
       & 0 < d < 2r,\\
     \pi r^2, & d \geq 2r,
   \end{cases}
   \qquad
   \theta=\cos^{-1}\left(1-\frac{d}{r}\right).

The discharge equation and head term are the same as for an orifice. This is a
simple circular-orifice representation; it does not choose among inlet-control
and outlet-control culvert regimes or account explicitly for barrel friction
and entrance, exit, or bend losses. Fold those effects into a calibrated flow
coefficient when appropriate.

``radial``
^^^^^^^^^^

A radial gate has invert ``elevation``, gate ``width``, and opening ``height``.
The implementation uses the same rectangular wetted area as an orifice,
:math:`A=W\min(H_{gate},z_h-z_{inv})`.

.. figure:: documents/image/radial_gate.png
   :alt: Radial gate dimensions and upstream and downstream water levels
   :align: center

   Radial gate geometry.

Define the submergence ratio

.. math::

   S = \frac{z_l-z_{inv}}{z_h-z_{inv}}.

SCHISM uses three regimes:

* For :math:`S < S_p=2/3`, the free-flow equation is

  .. math::

     Q_s = \operatorname{sgn}(z_u-z_d)\,C_{op}C_f A
           \sqrt{2g(z_h-z_{inv})}.

* For :math:`S_p \leq S < S_f=0.8`, it transitions linearly toward the
  submerged equation:

  .. math::

     \widehat S &= \frac{S-S_p}{S_f-S_p}, \\
     m &= \sqrt{\frac{1}{1-S_p}}, \\
     Q_s &= \operatorname{sgn}(z_u-z_d)\,C_{op}C_f A
            \sqrt{2g\Delta z}\left[(1-\widehat S)m+\widehat S\right].

* For :math:`S \geq S_f`, it uses the fully submerged orifice equation

  .. math::

     Q_s = \operatorname{sgn}(z_u-z_d)\,C_{op}C_f A\sqrt{2g\Delta z}.

``radial_relheight``
^^^^^^^^^^^^^^^^^^^^

This alternative radial-gate relation varies the flow coefficient linearly
with relative gate opening. In addition to ``coefficient`` :math:`d`, it takes
``coefficient_height`` :math:`s`:

.. math::

   R &= \min\left(\frac{H_{gate}}{z_h-z_{inv}},1\right), \\
   C_f &= d+sR, \\
   Q_s &= \operatorname{sgn}(z_u-z_d)\,C_{op}C_f A\sqrt{2g\Delta z}.

It uses the same rectangular wetted area as ``radial``. The native SCHISM type
name is ``radial_relheight``; the older name ``radial_rh`` is not accepted by
the current parser.

``weir_culvert``
^^^^^^^^^^^^^^^^

A compound structure combines independent weir and circular-culvert units that
share node pairs and reference water levels. SCHISM evaluates the weir and
culvert with their own geometry, coefficients, directional controls, and
numbers of duplicate units, then adds the signed discharges:

.. math::

   Q_s=Q_{weir}+Q_{culvert}.

The weir can therefore be dry while the lower culvert continues to flow. The
native input stores the weir as the main structure and the culvert as a
substructure.

.. _structure_preprocessing:

Geographic YAML input
---------------------

The schimpy preprocessor converts geographic structure definitions to global
node pairs and writes ``hydraulics.in``. Reference this file from the main
preprocessor input:

.. code-block:: yaml

   hydraulics:
     input file: hydraulics.yaml
     output file: hydraulics.in

A minimal structure file is:

.. code-block:: yaml

   nudging: 0.1
   structures:
     grantline_weir:
       type: weir
       end_points:
         - [636576.0, 4186940.0]
         - [636575.0, 4187000.0]
       configuration:
         n_duplicates: 1
         elevation: 1.0
         width: 54.9
         coefficient: 0.8
         op_downstream: 1.0
         op_upstream: 1.0
         use_time_series: true
       reference: self

``end_points`` gives the two ends of the structure span in the mesh coordinate
system. For a non-straight barrier, ``pathway`` may supply an ordered polyline;
``gate_span`` then identifies the part of that path used as the paired hydraulic
boundary. The ends of the selected span must reach land boundaries.

``reference`` is optional and defaults to ``self``. With ``self``, schimpy
selects reference nodes on opposite sides of this structure. A structure may
instead name another structure and reuse its reference pair, which is useful
for adjacent weir and culvert definitions that must use identical heads.

Configuration keys
^^^^^^^^^^^^^^^^^^

The following keys use underscores. They are values for one physical unit unless
noted otherwise.

.. list-table:: YAML configuration by type
   :header-rows: 1
   :widths: 22 55

   * - Type
     - Configuration keys
   * - ``transfer``
     - ``flow``, ``use_time_series``
   * - ``weir``
     - ``n_duplicates``, ``elevation``, ``width``, ``coefficient``,
       ``op_downstream``, ``op_upstream``, ``use_time_series``
   * - ``orifice`` or ``radial``
     - Weir keys plus ``height``
   * - ``radial_relheight``
     - Orifice/radial keys plus ``coefficient_height``
   * - ``culvert``
     - Weir keys, with ``radius`` in place of ``width``
   * - ``weir_culvert``
     - Weir keys plus ``culvert_n_duplicates``, ``culvert_elevation``,
       ``culvert_radius``, ``culvert_coefficient``,
       ``culvert_op_downstream``, and ``culvert_op_upstream``

For example, a culvert that uses the same reference nodes as the preceding weir
can be written as:

.. code-block:: yaml

   grantline_culvert:
     type: culvert
     end_points:
       - [636576.0, 4186918.0]
       - [636576.0, 4186940.0]
     configuration:
       n_duplicates: 6
       elevation: -1.28
       radius: 0.61
       coefficient: 0.6
       op_downstream: 0.0
       op_upstream: 1.0
       use_time_series: true
     reference: grantline_weir

.. _structure_input:

Native ``hydraulics.in`` input
------------------------------

Set ``ihydraulics = 1`` in the SCHISM parameter file to enable hydraulic
structures. The native file begins with two global records:

.. code-block:: text

   number_of_structures
   nudging

Each structure then has this common prefix:

.. code-block:: text

   index name
   number_of_node_pairs upstream_reference_node downstream_reference_node
   upstream_node_1 downstream_node_1
   ...
   upstream_node_N downstream_node_N
   structure_type
   n_duplicates

Indices must be sequential, names are limited to 32 non-whitespace characters,
and all node numbers are global. At least two node pairs are required. The type
specific records follow ``n_duplicates``.

.. list-table:: Type-specific native records
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Records after ``n_duplicates``
   * - ``transfer``
     - ``flow``
   * - ``weir`` or ``culvert``
     - ``elevation width_or_radius``; then
       ``coefficient op_downstream op_upstream``
   * - ``orifice`` or ``radial``
     - ``elevation width height``; then
       ``coefficient op_downstream op_upstream``
   * - ``radial_relheight``
     - ``elevation width height``; then
       ``coefficient coefficient_height``; then
       ``op_downstream op_upstream``
   * - ``weir_culvert``
     - Weir ``elevation width``; weir
       ``coefficient op_downstream op_upstream``; culvert
       ``n_duplicates``; culvert ``elevation radius``; then culvert
       ``coefficient op_downstream op_upstream``

The final record in every structure block is ``use_time_series``: zero disables
time-series control and any nonzero value enables it. Thus a complete native
weir block resembles:

.. code-block:: text

   1
   0.1
   1 grantline_weir
   3 120 126
   121 125
   122 124
   123 127
   weir
   1
   1.0 54.9
   0.8 1.0 1.0
   0

Time-series control
-------------------

When ``use_time_series`` is enabled, provide ``NAME.th`` in the SCHISM input
directory, where ``NAME`` exactly matches the structure name. The first column
is elapsed model time in seconds. The second is integer ``install``: zero
removes the structure and one installs it. Subsequent values replace selected
values from ``hydraulics.in``.

SCHISM advances these irregular time histories using stepwise changes rather
than linear interpolation. Include a record applicable at the beginning of the
run. The exact columns are:

.. list-table:: Structure time-series columns
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Columns
   * - ``transfer``
     - ``time install flow``
   * - ``weir`` or ``culvert``
     - ``time install n_duplicates op_downstream op_upstream elevation width_or_radius``
   * - ``orifice``, ``radial``, or ``radial_relheight``
     - ``time install n_duplicates op_downstream op_upstream elevation width height``
   * - ``weir_culvert``
     - ``time install weir_n_duplicates weir_op_downstream weir_op_upstream
       weir_elevation weir_width culvert_n_duplicates culvert_op_downstream
       culvert_op_upstream culvert_elevation culvert_radius``

The time series does not update ``coefficient`` or ``coefficient_height`` for
any structure type. For ``weir_culvert``, it does not update either component's
flow coefficient.

Implementation limits and modeling cautions
-------------------------------------------

The current implementation has several deliberate simplifications that should
be considered during calibration and model design:

* Reference head is water-surface elevation; velocity head is omitted.
* Culverts are circular orifices, not a full culvert inlet/outlet-control model.
* Orifices and culverts do not receive the weir's Villemonte submerged-flow
  correction.
* Radial-gate transition thresholds are fixed at :math:`2/3` and :math:`0.8`.
* Flow is distributed uniformly across the paired boundary rather than resolved
  through the physical opening.
* Time-series changes are stepwise, and flow relaxation means the imposed flow
  does not jump immediately to the newly calculated target.

These limits do not prevent calibration, but the fitted coefficient represents
all unresolved contraction, loss, and geometric effects. Validate flow direction,
reference-node wetness, crest or invert datum, duplicate count, and operation in
both directions before relying on a structure in production.

References
----------

Bodhaine, G. L. (1968). *Measurement of peak discharge at culverts by indirect
methods*. U.S. Geological Survey Techniques of Water-Resources Investigations,
Book 3, Chapter A3.

Rantz, S. E., and others (1982). *Measurement and computation of streamflow,
Volume 2: Computation of discharge*. U.S. Geological Survey Water-Supply Paper
2175.

Villemonte, J. R. (1947). Submerged-weir discharge studies. *Engineering News
Record*, 139(26), 54--56.
