"""Generate preprocessor inputs for gradually inundating restored islands.

A breach is described by two points along the levee, ordered left then right
looking into the island. From those the generator derives one polygon per breach
and reuses it verbatim across every artifact, so the dredge, the initial pool and
the hotstart region cannot drift apart.

The pool level follows from the dredge. :func:`schimpy.ellipse.ellipse` clamps its
radial term at 1, so every node inside a breach polygon ends up at ``min_depth`` or
deeper; ``-min_depth - freeboard`` is therefore a constant water level that meets
the island's ``-z-freeboard`` exactly at the ellipse rim.
"""

import logging
import os
from collections import deque

import click
import numpy as np

from schimpy import schism_yaml
from schimpy.schism_polygon import SchismPolygon, SchismPolygonDictConverter

__all__ = [
    "breach_axes",
    "local_edge_length",
    "breach_polygon",
    "infer_island_polygon",
    "build_inundation_artifacts",
    "write_inundation_inputs",
    "validate_structures",
    "generate_inundation_inputs",
]

FREEBOARD = 0.01
EDGE_LENGTH_FACTOR = 3.0
POLYGON_POINTS = 24
# push the ring just outside the ellipse so its rim really reaches min_depth
POLYGON_SCALE = 1.05
# a fill covering more than this share of the mesh has escaped the island
MAX_FILL_FRACTION = 0.5


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


def infer_island_polygon(
    mesh, breaches, max_fill_fraction=MAX_FILL_FRACTION, logger=None
):
    """Derive an island outline by filling inward from its breaches.

    Every breach line is a barrier, so a set of breaches together pens in the
    restoration area even where no single one of them closes it. Each breach
    contributes the reference node the structure itself will use on its island
    side, and all of them must land in one connected component. If they do not,
    the breaches disagree about which side the island is on, which almost always
    means one pair of points is ordered the wrong way round.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
        Mesh to walk.
    breaches : list of dict
        Resolved breaches carrying ``left``, ``right`` and ``name``.
    max_fill_fraction : float, optional
        Fail if the fill covers more than this share of the elements, which
        means it escaped past the levee.
    logger : logging.Logger, optional

    Returns
    -------
    list of list of float
        Outline vertices, first repeated as last.

    Raises
    ------
    ValueError
        If a breach has no island-side element, the breaches disagree, or the
        fill escapes.
    """
    from schimpy.schism_setup import SchismSetup

    blocked = set()
    for breach in breaches:
        blocked.update(_elements_crossing(mesh, breach["left"], breach["right"]))

    setup = SchismSetup(logger)
    setup.mesh = mesh

    seeds = {}
    for breach in breaches:
        seed, ref = _island_seed(mesh, setup, breach, blocked)
        if seed is None:
            raise ValueError(
                "Breach %s has no open element at its island-side reference node. "
                "Check that its left and right points are ordered looking into the "
                "island." % breach["name"]
            )
        seeds[breach["name"]] = (seed, ref)

    first = breaches[0]["name"]
    filled = _fill_from(mesh, seeds[first][0], blocked)

    for name, (seed, ref) in seeds.items():
        if not filled[seed]:
            x, y = mesh.nodes[ref][:2]
            raise ValueError(
                "Breach %s has its island-side reference node at (%.1f, %.1f), which "
                "is not in the same region as breach %s. Breaches that bound one "
                "island must all open onto it; check the left/right ordering."
                % (name, x, y, first)
            )

    n_elems = mesh.n_elems()
    n_filled = int(filled.sum())
    if n_filled > max_fill_fraction * n_elems:
        raise ValueError(
            "Island fill reached %d of %d elements, past the %g limit. The breaches "
            "do not pen in the island, so its polygon must be given explicitly."
            % (n_filled, n_elems, max_fill_fraction)
        )

    ring = _fill_outline(mesh, filled)
    if ring is None:
        raise ValueError(
            "Island fill of %d of %d elements produced no closed outline."
            % (n_filled, n_elems)
        )
    return ring


def _island_seed(mesh, setup, breach, blocked):
    """Return an open element at a breach's island-side reference node.

    Left to right looking into the island puts the interior on the ``down`` side
    of the node paths, so that side's reference node is the one the fill starts
    from. Reusing it means the fill and the structure agree by construction.
    """
    coords = np.array([breach["left"], breach["right"]], dtype=float)
    up_path, down_path = mesh.find_two_neighboring_node_paths(coords)
    if not up_path or not down_path:
        raise ValueError(
            "Breach %s does not cut the mesh; its points must straddle the opening."
            % breach["name"]
        )
    ref = setup._find_reference_node(down_path, up_path)
    for elem in mesh.get_elems_i_from_node(ref):
        if int(elem) not in blocked:
            return int(elem), ref
    return None, ref


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
    for key in ("left", "right", "min_depth", "max_depth"):
        if key not in breach:
            raise ValueError(
                "Breach %s/%s is missing %r." % (island_name, name, key)
            )

    left = np.asarray(breach["left"], dtype=float)[:2]
    right = np.asarray(breach["right"], dtype=float)[:2]
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
):
    """Build the preprocessor inputs for a set of restoration islands.

    Parameters
    ----------
    islands : list of dict
        Island specifications. Each needs ``name``, ``breaches`` and ``polygon``.
        The polygon is either the word ``infer``, or a mapping in the same form
        the preprocessor polygon files use, carrying at least ``vertices``.
        Each breach needs ``name``, ``left``, ``right``, ``min_depth`` and
        ``max_depth``, and may set ``major_axis_len`` and ``structure``.
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
            structures.append(
                {
                    "name": name,
                    "type": structure.get("type", "transfer"),
                    "end_points": [
                        [float(resolved["left"][0]), float(resolved["left"][1])],
                        [float(resolved["right"][0]), float(resolved["right"][1])],
                    ],
                    "configuration": structure.get("configuration", {"flow": 0.0}),
                    "reference": structure.get("reference", "self"),
                }
            )

    return {
        "structures": {"structures": structures},
        "depth_enforcement": {"polygons": dredge_polygons},
        "elevation": {"default": ambient, "polygons": elev_polygons},
        "regions": {"default": 0, "polygons": region_polygons},
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
        "# Use with allow_overlap: True; the domain entry is listed first so the\n"
        "# island regions win where they overlap it.\n"
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
        fname = _FILENAMES[key]
        if prefix:
            fname = "%s_%s" % (prefix, fname)
        path = os.path.join(out_dir, fname)
        with open(path, "w") as fh:
            fh.write(_HEADERS[key])
            schism_yaml.safe_dump(payload, fh, default_flow_style=False, sort_keys=False)
        written[key] = path
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
    spec_file,
    hgrid_file,
    out_dir,
    prefix=None,
    edge_length_factor=EDGE_LENGTH_FACTOR,
    freeboard=FREEBOARD,
    ambient=0.96,
    validate=True,
    logger=None,
):
    """Read a breach specification and write the preprocessor inputs.

    Parameters
    ----------
    spec_file : str
        Breach specification yaml, carrying an ``islands`` list.
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
    logger : logging.Logger, optional

    Returns
    -------
    dict
        Maps artifact name to the path written.
    """
    from schimpy.schism_mesh import read_mesh

    if logger is None:
        logger = logging.getLogger(__name__)

    with open(spec_file) as fh:
        spec = schism_yaml.load(fh)
    islands = spec.get("islands")
    if not islands:
        raise ValueError("%s has no 'islands' list." % spec_file)

    logger.info("Reading mesh %s", hgrid_file)
    mesh = read_mesh(hgrid_file)

    artifacts = build_inundation_artifacts(
        islands,
        mesh,
        edge_length_factor=edge_length_factor,
        freeboard=freeboard,
        ambient=ambient,
    )
    if validate:
        validate_structures(artifacts, hgrid_file, logger=logger)

    paths = write_inundation_inputs(artifacts, out_dir, prefix=prefix)
    for key, path in paths.items():
        logger.info("Wrote %s: %s", key, path)
    return paths


@click.command(
    help=(
        "Generate preprocessor inputs for gradually inundating restored islands.\n"
        "\n"
        "Reads a breach specification and writes four yaml files: hydraulic "
        "structures, depth enforcement for the dredged openings, an elev.ic "
        "polygon set, and regions for hotstart patch_init. Each breach is two "
        "points along the levee, ordered left then right looking into the island."
    )
)
@click.option(
    "--spec",
    "spec_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Breach specification yaml.",
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
    "--logdir", default=None, type=click.Path(file_okay=False), help="Log directory."
)
@click.option("--debug", is_flag=True, help="Verbose logging.")
@click.help_option("-h", "--help")
def inundate_island_cli(
    spec_file,
    hgrid_file,
    out_dir,
    prefix,
    edge_length_factor,
    freeboard,
    ambient,
    validate,
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
            spec_file=spec_file,
            hgrid_file=hgrid_file,
            out_dir=out_dir,
            prefix=prefix,
            edge_length_factor=edge_length_factor,
            freeboard=freeboard,
            ambient=ambient,
            validate=validate,
        )
    except (ValueError, OSError) as exc:
        raise click.ClickException(str(exc))

    click.echo("Generated:")
    for path in paths.values():
        click.echo("  %s" % path)


if __name__ == "__main__":
    inundate_island_cli()
