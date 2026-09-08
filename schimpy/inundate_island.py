"""Generate preprocessor inputs for gradually inundating restored islands.

A breach is described by a pathway along the levee, ordered looking into the
island, and optionally a shorter ``gate_span`` naming the part of it the
hydraulic structure actually occupies. Two points are the common case and may be
given as ``left``/``right`` instead. From the gate span the generator derives one
polygon per breach and reuses it verbatim across every artifact, so the dredge,
the initial pool and the hotstart region cannot drift apart. The full pathway is
what holds the island fill in.

The pool level follows from the dredge. :func:`schimpy.ellipse.ellipse` clamps its
radial term at 1, so every node inside a breach polygon ends up at ``min_depth`` or
deeper; ``-min_depth - freeboard`` is therefore a constant water level that meets
the island's ``-z-freeboard`` exactly at the ellipse rim.

A breach that names a ``breach_date`` also gets an operating schedule written as a
``.th`` pair, dated and elapsed. The date belongs to the breach, because a run
spanning several years may combine restoration sites that come online years apart.
It may also be set once on the island or once at the top of the file, and the
narrowest setting wins: breach, then island, then file, then ``--breach-date``.
"""

import datetime
import logging
import os
from collections import deque

import click
import numpy as np
import pandas as pd

from schimpy import schism_yaml
from schimpy.schism_polygon import SchismPolygon, SchismPolygonDictConverter

__all__ = [
    "breach_axes",
    "local_edge_length",
    "breach_polygon",
    "infer_island_polygon",
    "build_inundation_artifacts",
    "write_inundation_inputs",
    "build_breach_schedule",
    "write_breach_timeseries",
    "validate_structures",
    "generate_inundation_inputs",
]

FREEBOARD = 0.01
EDGE_LENGTH_FACTOR = 3.0
POLYGON_POINTS = 24
# push the ring just outside the ellipse so its rim really reaches min_depth
POLYGON_SCALE = 1.05

EPOCH_START = "2005-01-01"
RAMP_HOURS = 24.0
RAMP_STEPS = 12
RAMP_TARGET = 0.6

# Column order per structure type, matching the read in read_struct_ts()
# (hydraulic_structures.F90). Names are the schimpy configuration keys; the
# dated header uses the short aliases below.
_TS_COLUMNS = {
    "transfer": ["install", "flow"],
    "weir": [
        "install",
        "n_duplicates",
        "op_downstream",
        "op_upstream",
        "elevation",
        "width",
    ],
    "orifice": [
        "install",
        "n_duplicates",
        "op_downstream",
        "op_upstream",
        "elevation",
        "width",
        "height",
    ],
    "weir_culvert": [
        "install",
        "n_duplicates",
        "op_downstream",
        "op_upstream",
        "elevation",
        "width",
        "culvert_n_duplicates",
        "culvert_op_downstream",
        "culvert_op_upstream",
        "culvert_elevation",
        "culvert_width",
    ],
}
_TS_COLUMNS["culvert"] = _TS_COLUMNS["weir"]
_TS_COLUMNS["radial"] = _TS_COLUMNS["orifice"]
_TS_COLUMNS["radial_relheight"] = _TS_COLUMNS["orifice"]

_TS_ALIAS = {
    "n_duplicates": "ndup",
    "op_downstream": "op_down",
    "op_upstream": "op_up",
    "elevation": "elev",
    "culvert_n_duplicates": "cul_ndup",
    "culvert_op_downstream": "cul_op_down",
    "culvert_op_upstream": "cul_op_up",
    "culvert_elevation": "cul_elev",
    "culvert_width": "cul_width",
}

# Sensible stand-ins when the structure configuration omits a scheduled column.
_TS_DEFAULTS = {
    "install": 1,
    "n_duplicates": 1,
    "op_downstream": 1.0,
    "op_upstream": 1.0,
    "culvert_n_duplicates": 1,
    "culvert_op_downstream": 1.0,
    "culvert_op_upstream": 1.0,
}

_TS_INTEGER = {
    "install",
    "n_duplicates",
    "culvert_n_duplicates",
}


def breach_axes(left, right):
    """Return the centre and unit axes of a breach.

    Parameters
    ----------
    left, right : array_like
        The two levee points, ordered looking into the island.

    Returns
    -------
    center : numpy.ndarray
        Midpoint of the two points.
    along : numpy.ndarray
        Unit vector from left to right, the ellipse minor axis.
    into : numpy.ndarray
        Unit vector perpendicular to ``along``, pointing into the island.

    Raises
    ------
    ValueError
        If the two points coincide.
    """
    left = np.asarray(left, dtype=float)[:2]
    right = np.asarray(right, dtype=float)[:2]
    span = right - left
    length = float(np.hypot(*span))
    if length == 0.0:
        raise ValueError("Breach points coincide; they must straddle the opening.")
    along = span / length
    # left-hand normal of left->right, which is the island side by the ordering rule
    into = np.array([-along[1], along[0]])
    return 0.5 * (left + right), along, into


def local_edge_length(mesh, center, radius):
    """Return the median mesh edge length near a point.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
        Mesh to measure.
    center : array_like
        Point to measure around.
    radius : float
        Search radius.

    Returns
    -------
    float

    Raises
    ------
    ValueError
        If no edge midpoint falls within ``radius``.
    """
    center = np.asarray(center, dtype=float)[:2]
    edges = mesh.edges[:, :2]
    midpoints = 0.5 * (mesh.nodes[edges[:, 0], :2] + mesh.nodes[edges[:, 1], :2])
    near = np.flatnonzero(np.linalg.norm(midpoints - center, axis=1) <= radius)
    if near.size == 0:
        raise ValueError(
            "No mesh edge within %g of (%.1f, %.1f); is the breach on the mesh?"
            % (radius, center[0], center[1])
        )
    return float(np.median(mesh.edge_len()[near]))


def breach_polygon(left, right, major_axis_len, n_points=POLYGON_POINTS):
    """Return a closed ring enclosing the dredge ellipse for a breach.

    The ring is the ellipse itself scaled out slightly, so it encloses the taper
    without sweeping in the corners a bounding box would.

    Parameters
    ----------
    left, right : array_like
        The two levee points, ordered looking into the island.
    major_axis_len : float
        Full length of the axis perpendicular to ``left``-``right``.
    n_points : int, optional
        Number of ring vertices.

    Returns
    -------
    list of list of float
        Vertices, first repeated as last.
    """
    center, along, into = breach_axes(left, right)
    semi_minor = 0.5 * float(np.linalg.norm(np.asarray(right)[:2] - np.asarray(left)[:2]))
    semi_major = 0.5 * float(major_axis_len)

    angles = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    ring = [
        center
        + POLYGON_SCALE * (semi_minor * np.cos(a) * along + semi_major * np.sin(a) * into)
        for a in angles
    ]
    ring.append(ring[0])
    return [[round(float(p[0]), 3), round(float(p[1]), 3)] for p in ring]


def _segments_cross(p1, p2, q1, q2):
    """Vectorized proper-or-touching segment intersection test."""

    def cross(o, a, b):
        return (a[:, 0] - o[:, 0]) * (b[:, 1] - o[:, 1]) - (a[:, 1] - o[:, 1]) * (
            b[:, 0] - o[:, 0]
        )

    d1 = cross(q1, q2, p1)
    d2 = cross(q1, q2, p2)
    d3 = cross(p1, p2, q1)
    d4 = cross(p1, p2, q2)
    return ((d1 * d2) <= 0.0) & ((d3 * d4) <= 0.0)


def _edges_crossing(mesh, left, right):
    """Return indices of mesh edges the breach line crosses."""
    left = np.asarray(left, dtype=float)[:2]
    right = np.asarray(right, dtype=float)[:2]
    nodes = mesh.nodes[:, :2]
    edges = mesh.edges[:, :2].astype(int)
    a = nodes[edges[:, 0]]
    b = nodes[edges[:, 1]]

    span = np.linalg.norm(right - left)
    lo = np.minimum(left, right) - span
    hi = np.maximum(left, right) + span
    mid = 0.5 * (a + b)
    near = np.flatnonzero(
        np.all((mid >= lo) & (mid <= hi), axis=1)
    )
    if near.size == 0:
        return np.array([], dtype=int)

    p1 = np.repeat(left[None, :], near.size, axis=0)
    p2 = np.repeat(right[None, :], near.size, axis=0)
    hit = _segments_cross(p1, p2, a[near], b[near])
    return near[hit]


def _elements_crossing(mesh, left, right):
    """Return indices of elements the breach line passes through."""
    adjacency = mesh.edges[:, 3:5].astype(int)
    blocked = set()
    for edge_i in _edges_crossing(mesh, left, right):
        for elem in adjacency[int(edge_i)]:
            if int(elem) >= 0:
                blocked.add(int(elem))
    return blocked


def _pathway_crossing(mesh, pathway, name=None):
    """Return indices of elements any segment of a pathway passes through."""
    blocked = set()
    for i, (start, end) in enumerate(zip(pathway[:-1], pathway[1:])):
        hit = _elements_crossing(mesh, start, end)
        if not hit:
            raise ValueError(
                "Breach %s pathway segment %d, (%.1f, %.1f) to (%.1f, %.1f), crosses "
                "no element. A gap in the barrier lets the island fill escape."
                % (name, i, start[0], start[1], end[0], end[1])
            )
        blocked.update(hit)
    return blocked


def _breach_pathway(breach, island_name, name):
    """Return the levee polyline for a breach, as a list of 2-vectors.

    ``pathway`` is the general form; ``left``/``right`` is the two-point
    shorthand. Both are ordered so the island lies to the left.
    """
    has_pathway = "pathway" in breach
    has_pair = "left" in breach or "right" in breach
    if has_pathway and has_pair:
        raise ValueError(
            "Breach %s/%s gives both pathway and left/right; use one or the other."
            % (island_name, name)
        )
    if has_pathway:
        points = [np.asarray(p, dtype=float)[:2] for p in breach["pathway"]]
        if len(points) < 2:
            raise ValueError(
                "Breach %s/%s pathway needs at least two points."
                % (island_name, name)
            )
        return points
    for key in ("left", "right"):
        if key not in breach:
            raise ValueError(
                "Breach %s/%s is missing %r; give left and right, or a pathway."
                % (island_name, name, key)
            )
    return [
        np.asarray(breach["left"], dtype=float)[:2],
        np.asarray(breach["right"], dtype=float)[:2],
    ]


def _gate_span(breach, pathway, island_name, name):
    """Return the two points the hydraulic structure itself spans.

    The gate need not cover the whole levee pathway, but it must be a straight
    chord, since that is what the preprocessor resolves into node pairs.
    """
    if "gate_span" not in breach:
        if len(pathway) > 2:
            raise ValueError(
                "Breach %s/%s has a pathway of %d points, so it needs a gate_span "
                "naming the two ends of the structure itself; a bent pathway does "
                "not resolve into node pairs."
                % (island_name, name, len(pathway))
            )
        return pathway[0], pathway[1]

    span = [np.asarray(p, dtype=float)[:2] for p in breach["gate_span"]]
    if len(span) != 2:
        raise ValueError(
            "Breach %s/%s gate_span needs exactly two points, got %d."
            % (island_name, name, len(span))
        )
    left, right = span
    along = right - left
    if float(np.hypot(*along)) == 0.0:
        raise ValueError(
            "Breach %s/%s gate_span points coincide." % (island_name, name)
        )

    # the gate inherits the pathway's sense, so the island stays on the left
    center = 0.5 * (left + right)
    seg = min(
        zip(pathway[:-1], pathway[1:]),
        key=lambda s: float(np.linalg.norm(0.5 * (s[0] + s[1]) - center)),
    )
    if float(np.dot(along, seg[1] - seg[0])) < 0.0:
        raise ValueError(
            "Breach %s/%s gate_span runs against its pathway. Order both looking "
            "into the island." % (island_name, name)
        )
    return left, right


def infer_island_polygon(mesh, breaches, logger=None):
    """Derive an island outline by filling inward from its breaches.

    Every breach line is a barrier, so a set of breaches together pens in the
    restoration area even where no single one of them closes it. The whole
    pathway blocks, not just the gated part of it, so a short gate in a long
    levee still holds. Each breach contributes the reference node the structure
    itself will use on its island side, and the fill from every one of them must
    come out identical and must never reach the far side of any breach. If it
    does, the breaches disagree about which side the island is on, or the barrier
    they form has a gap.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
        Mesh to walk.
    breaches : list of dict
        Resolved breaches carrying ``pathway``, ``left``, ``right`` and ``name``.
    logger : logging.Logger, optional

    Returns
    -------
    list of list of float
        Outline vertices, first repeated as last.

    Raises
    ------
    ValueError
        If a breach has no island-side element, or the breaches disagree.
    """
    from schimpy.schism_setup import SchismSetup

    blocked = set()
    for breach in breaches:
        blocked.update(_pathway_crossing(mesh, breach["pathway"], breach["name"]))

    setup = SchismSetup(logger)
    setup.mesh = mesh

    seeds = {}
    for breach in breaches:
        inside, outside = _breach_seeds(mesh, setup, breach, blocked)
        if inside[0] is None:
            raise ValueError(
                "Breach %s has no open element at its island-side reference node. "
                "Check that its pathway and gate_span are ordered looking into the "
                "island." % breach["name"]
            )
        seeds[breach["name"]] = (inside, outside)

    first = breaches[0]["name"]
    filled = _fill_from(mesh, seeds[first][0][0], blocked)

    # agreement, not size, is what says the barrier holds
    for name, (inside, outside) in seeds.items():
        if name != first and not np.array_equal(
            _fill_from(mesh, inside[0], blocked), filled
        ):
            x, y = mesh.nodes[inside[1]][:2]
            raise ValueError(
                "Breach %s fills a different region than breach %s from its "
                "island-side reference node at (%.1f, %.1f). Breaches that bound "
                "one island must all open onto it; check the pathway ordering."
                % (name, first, x, y)
            )
        if outside[0] is not None and filled[outside[0]]:
            x, y = mesh.nodes[outside[1]][:2]
            raise ValueError(
                "The island fill reached (%.1f, %.1f), which is on the outside of "
                "breach %s, so it went around the barrier rather than being penned "
                "by it. The breaches do not pen in the island, so its polygon must "
                "be given explicitly." % (x, y, name)
            )

    n_filled = int(filled.sum())
    ring = _fill_outline(mesh, filled)
    if ring is None:
        raise ValueError(
            "Island fill of %d of %d elements produced no closed outline."
            % (n_filled, mesh.n_elems())
        )
    return ring


def _breach_seeds(mesh, setup, breach, blocked):
    """Return open elements at a breach's island-side and channel-side references.

    Ordering the pathway looking into the island puts the interior on the
    ``down`` side of the node paths, so that side's reference node is the one the
    fill starts from. Resolving it the same way the structure will means the fill
    and the structure agree by construction. The ``up`` side is the outside, and
    a fill that reaches it has escaped.
    """
    up_path, down_path = setup.structure_node_paths(
        breach["name"], breach["pathway"], breach.get("gate_span")
    )
    name = breach["name"]
    inside = setup._find_reference_node(down_path, up_path, name)
    outside = setup._find_reference_node(up_path, down_path, name)
    return (
        (_open_element(mesh, inside, blocked), inside),
        (_open_element(mesh, outside, blocked), outside),
    )


def _open_element(mesh, node_i, blocked):
    """Return an element at a node that the barrier does not pass through."""
    for elem in mesh.get_elems_i_from_node(node_i):
        if int(elem) not in blocked:
            return int(elem)
    return None


def _fill_from(mesh, seed, blocked):
    """Flood the element graph from a seed, refusing to enter blocked elements."""
    filled = np.zeros(mesh.n_elems(), dtype=bool)
    filled[seed] = True
    adjacency = mesh.edges[:, 3:5].astype(int)
    queue = deque([seed])
    while queue:
        elem = queue.popleft()
        for edge_i in mesh.element2edges(elem):
            pair = adjacency[int(edge_i)]
            other = int(pair[1]) if int(pair[0]) == elem else int(pair[0])
            if other < 0 or filled[other] or other in blocked:
                continue
            filled[other] = True
            queue.append(other)
    return filled


def _fill_outline(mesh, filled, pad_fraction=0.25):
    """Trace the longest ring bounding a filled element set.

    The ring runs through mesh nodes, and ``contains`` excludes the boundary, so
    it is pushed out by a fraction of an edge to keep those nodes inside.
    """
    from shapely.geometry import MultiPolygon, Polygon

    adjacency = mesh.edges[:, 3:5].astype(int)
    a = adjacency[:, 0]
    b = adjacency[:, 1]
    inside_a = np.where(a >= 0, filled[np.clip(a, 0, None)], False)
    inside_b = np.where(b >= 0, filled[np.clip(b, 0, None)], False)
    boundary = np.flatnonzero(inside_a.astype(int) + inside_b.astype(int) == 1)
    if boundary.size == 0:
        return None

    links = {}
    for edge_i in boundary:
        n1, n2 = int(mesh.edges[edge_i][0]), int(mesh.edges[edge_i][1])
        links.setdefault(n1, []).append(n2)
        links.setdefault(n2, []).append(n1)

    best = None
    unused = set(links)
    while unused:
        start = unused.pop()
        walk = [start]
        prev, cur = None, start
        while True:
            options = [n for n in links[cur] if n != prev]
            if not options:
                break
            nxt = options[0]
            if nxt == start:
                break
            walk.append(nxt)
            unused.discard(nxt)
            prev, cur = cur, nxt
            if len(walk) > len(links):
                break
        if best is None or len(walk) > len(best):
            best = walk

    if best is None or len(best) < 3:
        return None

    poly = Polygon([mesh.nodes[n][:2] for n in best])
    if not poly.is_valid:
        poly = poly.buffer(0.0)
    pad = pad_fraction * float(np.median(mesh.edge_len()[boundary]))
    poly = poly.buffer(pad)
    if isinstance(poly, MultiPolygon):
        poly = max(poly.geoms, key=lambda g: g.area)
    return [
        [round(float(x), 3), round(float(y), 3)] for x, y in poly.exterior.coords
    ]


def _island_polygon(island, mesh, resolved_breaches):
    """Resolve an island outline from the spec, inferring it when asked."""
    spec = island.get("polygon")
    name = island["name"]
    if spec is None:
        raise ValueError(
            "Island %r needs a polygon; give vertices or the word 'infer'." % name
        )
    if isinstance(spec, str):
        if spec.strip().lower() != "infer":
            raise ValueError(
                "Island %r polygon must be a polygon mapping or 'infer', got %r."
                % (name, spec)
            )
        vertices = infer_island_polygon(mesh, resolved_breaches)
        return SchismPolygon(vertices, prop={"name": name, "type": "none"})

    if not isinstance(spec, dict) or "vertices" not in spec:
        raise ValueError(
            "Island %r polygon must be a mapping with 'vertices', in the same form "
            "used by the preprocessor polygon files." % name
        )
    item = dict(spec)
    item.setdefault("name", name)
    item.setdefault("type", "none")
    return SchismPolygonDictConverter().read({"polygons": [item]})[0]


def _ring_of(polygon):
    """Return the exterior ring of a SchismPolygon as plain vertices."""
    return [
        [round(float(x), 3), round(float(y), 3)]
        for x, y in polygon.exterior.coords
    ]


def _points(seq):
    """Return points as plain nested floats, fit for yaml."""
    return [[float(p[0]), float(p[1])] for p in seq]


def _dredge_attribute(left, right, min_depth, max_depth, major_axis_len):
    return (
        "ellipse(x, y, z, [{:.1f}, {:.1f}], [{:.1f}, {:.1f}], "
        "min_depth={:.3f}, max_depth={:.3f}, major_axis_len={:.1f})".format(
            left[0], left[1], right[0], right[1], min_depth, max_depth, major_axis_len
        )
    )


def _resolve_breach(breach, island_name, mesh, edge_length_factor):
    """Fill in the derived geometry for one breach."""
    name = breach.get("name")
    if not name:
        raise ValueError("Every breach in island %r needs a name." % island_name)
    for key in ("min_depth", "max_depth"):
        if key not in breach:
            raise ValueError(
                "Breach %s/%s is missing %r." % (island_name, name, key)
            )

    pathway = _breach_pathway(breach, island_name, name)
    left, right = _gate_span(breach, pathway, island_name, name)
    min_depth = float(breach["min_depth"])
    max_depth = float(breach["max_depth"])
    if max_depth <= min_depth:
        raise ValueError(
            "Breach %s/%s has max_depth %g no deeper than min_depth %g; depths are "
            "positive down." % (island_name, name, max_depth, min_depth)
        )

    center, along, _ = breach_axes(left, right)
    span = float(np.linalg.norm(right - left))
    edge = local_edge_length(mesh, center, span)

    if "major_axis_len" in breach:
        # taken as the full axis, matching the ellipse signature
        major_axis_len = float(breach["major_axis_len"])
    else:
        major_axis_len = 2.0 * edge_length_factor * edge

    # the dredge runs a little past the opening, so the structure's own end nodes
    # and its reference pair sit on dredged ground rather than on the taper
    pad = float(breach.get("dredge_pad", edge))
    dredge_left = left - along * pad
    dredge_right = right + along * pad

    return {
        "name": "%s_%s" % (island_name, name),
        "left": left,
        "right": right,
        "pathway": pathway,
        "gate_span": [left, right] if "gate_span" in breach else None,
        "dredge_left": dredge_left,
        "dredge_right": dredge_right,
        "min_depth": min_depth,
        "max_depth": max_depth,
        "major_axis_len": major_axis_len,
        "edge": edge,
        "vertices": breach_polygon(dredge_left, dredge_right, major_axis_len),
    }


def build_inundation_artifacts(
    islands,
    mesh,
    edge_length_factor=EDGE_LENGTH_FACTOR,
    freeboard=FREEBOARD,
    ambient=0.96,
    breach_date=None,
):
    """Build the preprocessor inputs for a set of restoration islands.

    Parameters
    ----------
    islands : list of dict
        Island specifications. Each needs ``name``, ``breaches`` and ``polygon``.
        The polygon is either the word ``infer``, or a mapping in the same form
        the preprocessor polygon files use, carrying at least ``vertices``.
        Each breach needs ``name``, ``min_depth``, ``max_depth`` and its geometry
        as either ``pathway`` or ``left``/``right``. A ``pathway`` of more than
        two points also needs ``gate_span``, the straight chord the structure
        occupies. ``major_axis_len``, ``structure`` and ``breach_date`` are
        optional; ``breach_date`` may also be given once per island.
    mesh : schimpy.schism_mesh.SchismMesh
        Mesh the breaches sit on, used to size the dredge footprint.
    edge_length_factor : float, optional
        Dredge penetration each way from the levee line, as a multiple of the
        local edge length. The ellipse major axis is twice this.
    freeboard : float, optional
        Depth the initial surface sits below the bed on dry ground.
    ambient : float, optional
        Water level away from the islands, used for the domain entry of the
        elevation file.

    Returns
    -------
    dict
        Maps artifact name to a structure ready for yaml serialization:
        ``structures``, ``depth_enforcement``, ``elevation`` and ``regions``.

    Raises
    ------
    ValueError
        If a specification is incomplete or geometrically inconsistent.
    """
    structures = []
    schedules = []
    dredge_polygons = []
    elev_polygons = []
    region_polygons = []

    bounds = mesh.nodes[:, :2]
    pad = 0.05 * float(np.ptp(bounds, axis=0).max())
    lo = bounds.min(axis=0) - pad
    hi = bounds.max(axis=0) + pad
    domain_ring = [
        [float(lo[0]), float(lo[1])],
        [float(hi[0]), float(lo[1])],
        [float(hi[0]), float(hi[1])],
        [float(lo[0]), float(hi[1])],
        [float(lo[0]), float(lo[1])],
    ]
    elev_polygons.append(
        {
            "name": "domain",
            "type": "none",
            "attribute": "max(%s, -z-%s)" % (ambient, freeboard),
            "vertices": domain_ring,
        }
    )
    # listed first so island regions win the overlap under allow_overlap
    region_polygons.append(
        {"name": "domain", "attribute": 0, "vertices": domain_ring}
    )

    seen = set()
    for island in islands:
        island_name = island.get("name")
        if not island_name:
            raise ValueError("Every island needs a name.")
        if island_name in seen:
            raise ValueError("Duplicate island name %r." % island_name)
        seen.add(island_name)
        if any(c.isspace() for c in island_name):
            raise ValueError("Island name %r cannot contain whitespace." % island_name)

        breaches = island.get("breaches") or []
        if not breaches:
            raise ValueError("Island %r has no breaches." % island_name)
        resolved_breaches = [
            _resolve_breach(b, island_name, mesh, edge_length_factor) for b in breaches
        ]

        island_ring = _ring_of(_island_polygon(island, mesh, resolved_breaches))

        elev_polygons.append(
            {
                "name": island_name,
                "type": "none",
                "attribute": "-z-%s" % freeboard,
                "vertices": island_ring,
            }
        )
        region_polygons.append(
            {"name": island_name, "attribute": 1, "vertices": island_ring}
        )

        for breach, resolved in zip(breaches, resolved_breaches):
            name = resolved["name"]
            if name in seen:
                raise ValueError("Duplicate breach name %r." % name)
            seen.add(name)

            dredge_polygons.append(
                {
                    "name": name,
                    "type": "none",
                    "attribute": _dredge_attribute(
                        resolved["dredge_left"],
                        resolved["dredge_right"],
                        resolved["min_depth"],
                        resolved["max_depth"],
                        resolved["major_axis_len"],
                    ),
                    "vertices": resolved["vertices"],
                }
            )
            elev_polygons.append(
                {
                    "name": name,
                    "type": "min",
                    "attribute": round(-resolved["min_depth"] - freeboard, 4),
                    "vertices": resolved["vertices"],
                }
            )

            structure = breach.get("structure") or {}
            configuration = dict(
                structure.get("configuration", {"flow": 0.0})
            )
            struct_type = structure.get("type", "transfer")
            when = breach.get("breach_date", island.get("breach_date", breach_date))
            if when is not None:
                # SCHISM looks for <struct_name>.th only when this flag is set
                configuration["use_time_series"] = 1
                schedules.append(
                    {
                        "name": name,
                        "type": struct_type,
                        "configuration": configuration,
                        "breach_date": when,
                    }
                )
            entry = {
                "name": name,
                "type": struct_type,
                "configuration": configuration,
                "reference": structure.get("reference", "self"),
            }
            pathway = resolved["pathway"]
            gate_span = resolved["gate_span"]
            if len(pathway) == 2 and gate_span is None:
                entry["end_points"] = _points(pathway)
            else:
                entry["pathway"] = _points(pathway)
                if gate_span is not None:
                    entry["gate_span"] = _points(gate_span)
            structures.append(entry)

    return {
        "structures": {"structures": structures},
        "depth_enforcement": {"polygons": dredge_polygons},
        "elevation": {"default": ambient, "polygons": elev_polygons},
        "regions": {"default": 0, "polygons": region_polygons},
        "schedules": schedules,
    }


_HEADERS = {
    "structures": (
        "# Generated by schimpy.inundate_island. Include alongside the base\n"
        "# hydraulics file; this file deliberately carries no nudging value.\n"
    ),
    "depth_enforcement": (
        "# Generated by schimpy.inundate_island.\n"
        "# Requires 'schimpy.ellipse.ellipse' among the preprocessor imports.\n"
    ),
    "elevation": (
        "# Generated by schimpy.inundate_island.\n"
        "# Relies on order. Breach depth comes after general island depth and thus\n"
        "# takes precedence. The breach entries are type 'min', a lower bound, so they\n"
        "# pool the dredged opening while leaving ambient water and dry ground alone.\n"
    ),
    "regions": (
        "# Generated by schimpy.inundate_island for hotstart patch_init.\n"
        "# This is a complete partition only after overlap resolution: use\n"
        "# allow_overlap: True and list domain first in patch_init.regions, followed\n"
        "# by the restoration regions. Later configured regions take precedence.\n"
    ),
}

_FILENAMES = {
    "structures": "hydraulic_structures_inundate.yaml",
    "depth_enforcement": "depth_enforce_inundate.yaml",
    "elevation": "elev_inundate.yaml",
    "regions": "inundate_regions.yaml",
}


def write_inundation_inputs(artifacts, out_dir, prefix=None):
    """Write generated artifacts as yaml.

    Parameters
    ----------
    artifacts : dict
        Output of :func:`build_inundation_artifacts`.
    out_dir : str
        Directory to write into. Created if absent.
    prefix : str, optional
        Prepended to each filename.

    Returns
    -------
    dict
        Maps artifact name to the path written.
    """
    os.makedirs(out_dir, exist_ok=True)
    written = {}
    for key, payload in artifacts.items():
        if key not in _FILENAMES:
            continue
        fname = _FILENAMES[key]
        if prefix:
            fname = "%s_%s" % (prefix, fname)
        path = os.path.join(out_dir, fname)
        with open(path, "w") as fh:
            fh.write(_HEADERS[key])
            schism_yaml.safe_dump(payload, fh, default_flow_style=False, sort_keys=False)
        written[key] = path
    return written


def build_breach_schedule(
    schedule,
    epoch_start=EPOCH_START,
    ramp_hours=RAMP_HOURS,
    ramp_steps=RAMP_STEPS,
    ramp_target=RAMP_TARGET,
):
    """Build the operating schedule for one breach as a dated table.

    The structure holds its configured setting from ``epoch_start`` until
    ``ramp_hours`` before the breach, opens over ``ramp_steps`` evenly spaced
    increments, and is then deinstalled at the breach itself. Deinstalling
    returns the block elements to the momentum solve; ramping first keeps that
    switch small, and ``block_nudge`` relaxes what is left.

    Parameters
    ----------
    schedule : dict
        One entry from ``artifacts['schedules']``, carrying ``name``, ``type``,
        ``configuration`` and ``breach_date``.
    epoch_start : str or datetime, optional
        First row, standing in for "infinitely in the past".
    ramp_hours : float, optional
        Length of the opening ramp before the breach.
    ramp_steps : int, optional
        Number of increments in the ramp.
    ramp_target : float, optional
        Value the ramped column reaches at the end of the ramp.

    Returns
    -------
    pandas.DataFrame
        Indexed by datetime, with one column per scheduled field.
    """
    struct_type = schedule["type"]
    try:
        columns = _TS_COLUMNS[struct_type]
    except KeyError:
        raise ValueError(
            "No time series layout known for structure type %r." % struct_type
        )

    configuration = schedule["configuration"]
    missing = [
        c
        for c in columns
        if c not in configuration and c not in _TS_DEFAULTS
    ]
    if missing:
        raise ValueError(
            "Structure %r is missing %s, needed for its time series."
            % (schedule["name"], ", ".join(missing))
        )
    base = {c: configuration.get(c, _TS_DEFAULTS.get(c)) for c in columns}
    base["install"] = 1

    ramped = "flow" if struct_type == "transfer" else "op_downstream"
    start = float(base[ramped])

    breach = pd.Timestamp(schedule["breach_date"])
    epoch = pd.Timestamp(epoch_start)
    if epoch >= breach - pd.Timedelta(hours=ramp_hours):
        raise ValueError(
            "epoch_start %s must precede the ramp into breach_date %s."
            % (epoch, breach)
        )
    if ramp_steps < 1:
        raise ValueError("ramp_steps must be at least 1.")

    times = [epoch]
    rows = [dict(base)]

    step = pd.Timedelta(hours=ramp_hours) / ramp_steps
    for k in range(ramp_steps):
        row = dict(base)
        row[ramped] = start + (float(ramp_target) - start) * (k + 1) / ramp_steps
        times.append(breach - pd.Timedelta(hours=ramp_hours) + k * step)
        rows.append(row)

    final = dict(rows[-1])
    final["install"] = 0
    times.append(breach)
    rows.append(final)

    frame = pd.DataFrame(rows, index=pd.DatetimeIndex(times), columns=columns)
    for col in columns:
        if col in _TS_INTEGER:
            frame[col] = frame[col].astype(int)
        else:
            frame[col] = frame[col].astype(float)
    return frame


def _format_schedule(frame, dated):
    """Render a schedule as whitespace separated text."""
    out = frame.copy()
    for col in out.columns:
        if col in _TS_INTEGER:
            out[col] = out[col].map("{:d}".format)
        else:
            out[col] = out[col].map(lambda v: np.format_float_positional(
                v, precision=6, unique=True, trim="0", fractional=True
            ))
    if dated:
        out.index = [t.strftime("%Y-%m-%dT%H:%M") for t in frame.index]
        header = ["datetime"] + [_TS_ALIAS.get(c, c) for c in out.columns]
        lines = ["\t".join(header)]
    else:
        # SCHISM reads these with a list-directed read, so no header
        out.index = ["{:.1f}".format(v) for v in frame.index]
        lines = []
    for label, row in out.iterrows():
        lines.append("\t".join([label] + list(row.values)))
    return "\n".join(lines) + "\n"


def write_breach_timeseries(
    schedules,
    out_dir,
    run_start,
    th_dir="th_files",
    epoch_start=EPOCH_START,
    ramp_hours=RAMP_HOURS,
    ramp_steps=RAMP_STEPS,
    ramp_target=RAMP_TARGET,
    logger=None,
):
    """Write dated and elapsed ``.th`` files for every scheduled breach.

    Parameters
    ----------
    schedules : list of dict
        ``artifacts['schedules']``.
    out_dir : str
        Directory the artifacts were written into.
    run_start : str or datetime
        Model time origin the elapsed files are measured from.
    th_dir : str, optional
        Subdirectory holding ``dated/`` and ``elapsed/``.
    epoch_start, ramp_hours, ramp_steps, ramp_target
        Passed through to :func:`build_breach_schedule`.
    logger : logging.Logger, optional

    Returns
    -------
    dict
        Maps structure name to ``(dated_path, elapsed_path)``.
    """
    from vtools.data.timeseries import datetime_elapsed

    if logger is None:
        logger = logging.getLogger(__name__)

    dated_dir = os.path.join(out_dir, th_dir, "dated")
    elapsed_dir = os.path.join(out_dir, th_dir, "elapsed")
    os.makedirs(dated_dir, exist_ok=True)
    os.makedirs(elapsed_dir, exist_ok=True)

    reference = pd.Timestamp(run_start)
    written = {}
    for schedule in schedules:
        frame = build_breach_schedule(
            schedule,
            epoch_start=epoch_start,
            ramp_hours=ramp_hours,
            ramp_steps=ramp_steps,
            ramp_target=ramp_target,
        )
        name = schedule["name"]
        dated_path = os.path.join(dated_dir, "%s.th" % name)
        with open(dated_path, "w") as fh:
            fh.write(_format_schedule(frame, dated=True))

        elapsed = datetime_elapsed(frame, reftime=reference, dtype="d")
        elapsed_path = os.path.join(elapsed_dir, "%s.th" % name)
        with open(elapsed_path, "w") as fh:
            fh.write(_format_schedule(elapsed, dated=False))

        logger.info(
            "Breach %s: %d rows, breach at %s, elapsed %.0f s from %s",
            name,
            len(frame),
            pd.Timestamp(schedule["breach_date"]),
            elapsed.index[-1],
            reference,
        )
        written[name] = (dated_path, elapsed_path)
    return written


def validate_structures(artifacts, hgrid_file, logger=None):
    """Check that each generated structure resolves inside its own breach polygon.

    Builds the structures the way the preprocessor will, then confirms the node
    pairs and the reference pair land inside the polygon that dredges the opening.
    A reference node outside it would sit on undredged ground and could go dry.

    Parameters
    ----------
    artifacts : dict
        Output of :func:`build_inundation_artifacts`.
    hgrid_file : str
        Mesh the structures are placed on.
    logger : logging.Logger, optional
        Passed to the setup so failures report through the usual channel.

    Returns
    -------
    dict
        Maps structure name to the number of nodes checked.

    Raises
    ------
    ValueError
        If a structure cannot be built, or a node falls outside its polygon.
    """
    from shapely.geometry import Point, Polygon

    from schimpy.schism_setup import create_schism_setup

    if logger is None:
        logger = logging.getLogger(__name__)
    setup = create_schism_setup(hgrid_file, logger)
    setup.create_structures(artifacts["structures"]["structures"])

    rings = {
        p["name"]: Polygon(p["vertices"])
        for p in artifacts["depth_enforcement"]["polygons"]
    }

    checked = {}
    for built in setup._input.structures:
        ring = rings[built.name]
        nodes = [n for pair in built.node_pairs for n in pair]
        nodes.extend(built.reference_pair)
        outside = [n for n in nodes if not ring.contains(Point(setup.mesh.nodes[n][:2]))]
        if outside:
            coords = ", ".join(
                "(%.1f, %.1f)" % tuple(setup.mesh.nodes[n][:2]) for n in outside[:5]
            )
            raise ValueError(
                "Structure %s has %d node(s) outside its dredge polygon: %s. Widen "
                "major_axis_len or move the breach points."
                % (built.name, len(outside), coords)
            )
        checked[built.name] = len(nodes)
        logger.info(
            "Structure %s: %d node pairs, reference pair inside the dredge polygon",
            built.name,
            len(built.node_pairs),
        )
    return checked


def generate_inundation_inputs(
    config_file,
    hgrid_file,
    out_dir,
    prefix=None,
    edge_length_factor=EDGE_LENGTH_FACTOR,
    freeboard=FREEBOARD,
    ambient=0.96,
    validate=True,
    breach_date=None,
    run_start=None,
    th_dir="th_files",
    epoch_start=EPOCH_START,
    ramp_hours=RAMP_HOURS,
    ramp_steps=RAMP_STEPS,
    ramp_target=RAMP_TARGET,
    logger=None,
):
    """Read a breach configuration and write the preprocessor inputs.

    Parameters
    ----------
    config_file : str
        Breach configuration yaml, carrying an ``islands`` list.
    hgrid_file : str
        Mesh the breaches sit on.
    out_dir : str
        Directory to write into.
    prefix : str, optional
        Prepended to each output filename.
    edge_length_factor : float, optional
        Dredge penetration each way from the levee line, as a multiple of the
        local edge length. The ellipse major axis is twice this.
    freeboard : float, optional
        Depth the initial surface sits below the bed on dry ground.
    ambient : float, optional
        Water level away from the islands.
    validate : bool, optional
        Build the structures and confirm they resolve inside their polygons.
    breach_date : str or datetime, optional
        Last-resort breach moment for breaches that name neither their own nor
        their island's. A breach with no date resolves to no time series.
    run_start : str or datetime, optional
        Model time origin. Required to write the elapsed ``.th`` files.
    th_dir : str, optional
        Subdirectory for ``dated/`` and ``elapsed/``.
    epoch_start, ramp_hours, ramp_steps, ramp_target
        Passed through to :func:`build_breach_schedule`.
    logger : logging.Logger, optional

    Returns
    -------
    dict
        Maps artifact name to the path written.
    """
    from schimpy.schism_mesh import read_mesh

    if logger is None:
        logger = logging.getLogger(__name__)

    with open(config_file) as fh:
        config = schism_yaml.load(fh)
    islands = config.get("islands")
    if not islands:
        raise ValueError("%s has no 'islands' list." % config_file)

    logger.info("Reading mesh %s", hgrid_file)
    mesh = read_mesh(hgrid_file)

    artifacts = build_inundation_artifacts(
        islands,
        mesh,
        edge_length_factor=edge_length_factor,
        freeboard=freeboard,
        ambient=ambient,
        breach_date=breach_date or config.get("breach_date"),
    )
    if validate:
        validate_structures(artifacts, hgrid_file, logger=logger)

    paths = write_inundation_inputs(artifacts, out_dir, prefix=prefix)
    for key, path in paths.items():
        logger.info("Wrote %s: %s", key, path)

    schedules = artifacts.get("schedules") or []
    if schedules:
        origin = run_start or config.get("run_start")
        if origin is None:
            raise ValueError(
                "A breach_date was given, so run_start is needed to write the "
                "elapsed time series."
            )
        written = write_breach_timeseries(
            schedules,
            out_dir,
            origin,
            th_dir=th_dir,
            epoch_start=epoch_start,
            ramp_hours=ramp_hours,
            ramp_steps=ramp_steps,
            ramp_target=ramp_target,
            logger=logger,
        )
        for name, (dated_path, elapsed_path) in written.items():
            paths["%s_dated" % name] = dated_path
            paths["%s_elapsed" % name] = elapsed_path
    return paths


@click.command(
    help=(
        "Generate preprocessor inputs for gradually inundating restored islands.\n"
        "\n"
        "Reads a breach configuration and writes four yaml files: hydraulic "
        "structures, depth enforcement for the dredged openings, an elev.ic "
        "polygon set, and regions for hotstart patch_init. Each breach is a "
        "pathway along the levee, ordered looking into the island, and "
        "optionally a shorter gate_span the structure itself occupies."
    )
)
@click.option(
    "--config",
    "config_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Breach configuration yaml.",
)
@click.option(
    "--hgrid",
    "hgrid_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Mesh the breaches sit on.",
)
@click.option(
    "--out-dir",
    default=".",
    type=click.Path(file_okay=False),
    help="Directory for the generated yaml (default: current directory).",
)
@click.option("--prefix", default=None, type=str, help="Prefix for output filenames.")
@click.option(
    "--edge-length-factor",
    default=EDGE_LENGTH_FACTOR,
    type=float,
    show_default=True,
    help="Dredge penetration each way from the levee, in local edge lengths.",
)
@click.option(
    "--freeboard",
    default=FREEBOARD,
    type=float,
    show_default=True,
    help="Depth the initial surface sits below the bed on dry ground.",
)
@click.option(
    "--ambient",
    default=0.96,
    type=float,
    show_default=True,
    help="Water level away from the islands.",
)
@click.option(
    "--validate/--no-validate",
    default=True,
    show_default=True,
    help="Build the structures and confirm they resolve inside their polygons.",
)
@click.option(
    "--breach-date",
    default=None,
    type=str,
    help=(
        "Last-resort breach moment, e.g. 2022-04-15. Prefer setting breach_date "
        "on each breach, since sites in one run can come online years apart."
    ),
)
@click.option(
    "--run-start",
    default=None,
    type=str,
    help="Model time origin for the elapsed .th, e.g. 2020-09-30.",
)
@click.option(
    "--th-dir",
    default="th_files",
    show_default=True,
    type=str,
    help="Subdirectory holding dated/ and elapsed/.",
)
@click.option(
    "--epoch-start",
    default=EPOCH_START,
    show_default=True,
    type=str,
    help="First schedule row, standing in for the infinite past.",
)
@click.option(
    "--ramp-hours",
    default=RAMP_HOURS,
    show_default=True,
    type=float,
    help="Length of the opening ramp before the breach.",
)
@click.option(
    "--ramp-steps",
    default=RAMP_STEPS,
    show_default=True,
    type=int,
    help="Number of increments in the opening ramp.",
)
@click.option(
    "--ramp-target",
    default=RAMP_TARGET,
    show_default=True,
    type=float,
    help="Value the ramped column reaches at the end of the ramp.",
)
@click.option(
    "--logdir", default=None, type=click.Path(file_okay=False), help="Log directory."
)
@click.option("--debug", is_flag=True, help="Verbose logging.")
@click.help_option("-h", "--help")
def inundate_island_cli(
    config_file,
    hgrid_file,
    out_dir,
    prefix,
    edge_length_factor,
    freeboard,
    ambient,
    validate,
    breach_date,
    run_start,
    th_dir,
    epoch_start,
    ramp_hours,
    ramp_steps,
    ramp_target,
    logdir,
    debug,
):
    """CLI wrapper for :func:`generate_inundation_inputs`."""
    handlers = [logging.StreamHandler()]
    if logdir:
        os.makedirs(logdir, exist_ok=True)
        handlers.append(
            logging.FileHandler(os.path.join(logdir, "inundate_island.log"))
        )
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(levelname)s %(message)s",
        handlers=handlers,
    )

    try:
        paths = generate_inundation_inputs(
            config_file=config_file,
            hgrid_file=hgrid_file,
            out_dir=out_dir,
            prefix=prefix,
            edge_length_factor=edge_length_factor,
            freeboard=freeboard,
            ambient=ambient,
            validate=validate,
            breach_date=breach_date,
            run_start=run_start,
            th_dir=th_dir,
            epoch_start=epoch_start,
            ramp_hours=ramp_hours,
            ramp_steps=ramp_steps,
            ramp_target=ramp_target,
        )
    except (ValueError, OSError) as exc:
        raise click.ClickException(str(exc))

    click.echo("Generated:")
    for path in paths.values():
        click.echo("  %s" % path)


if __name__ == "__main__":
    inundate_island_cli()
