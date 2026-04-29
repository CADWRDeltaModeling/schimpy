#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Relocate co-located SCHISM source/sink point pairs.

Overview
--------
SCHISM source and sink points that share an element create two problems:

1. The element can dry, in which case SCHISM silently drops the source.
2. Co-located withdrawals and discharges short-circuit ambient mixing because
   the model treats them as cancelling at the same control volume.

This module reads an ``hgrid.gr3`` mesh together with a paired set of source
and sink point YAMLs (one prefix family at a time, ``delta`` or ``suisun``)
and offsets every still-co-located pair to nearby distinct elements that
satisfy depth criteria. The deeper of the two chosen elements is assigned to
the sink and the shallower to the source, except that when both candidates
are deep (mean depth ``>= --random-threshold``) the assignment is randomized
so we do not bias the choice on rerun.

Search strategy
---------------
The search is in **element-ring hops** (BFS over element adjacency through
shared nodes), not in meters:

* ``--rings-initial`` (default 3) ring-hops are tried first.
* If no element satisfies ``min-node-depth`` and ``mean-depth``, the search
  expands to ``--rings-fallback`` (default 10) ring-hops.
* A physical safety cap ``--max-distance`` (default 100 m) further rejects
  candidates whose centroid is farther than that radius from the seed
  centroid (which is the original point when on-grid, or the snapped
  nearest-element centroid when off-grid).

Off-grid and previously-separated inputs
----------------------------------------
Each endpoint of every pair is validated independently against the depth
criteria, regardless of whether the pair is co-located or already separated
beyond ``--separation-tol``:

* On-grid AND qualifying → coordinate is preserved verbatim.
* On-grid but too shallow (``min_node_depth`` / ``mean_depth`` not met) →
  the endpoint is relocated via element-ring search to the nearest
  qualifying element. This matters when a previously-validated source/sink
  set is re-run on a *swapped mesh* where prior separation no longer
  implies prior depth qualification.
* Off-grid (``mesh.find_elem`` returns ``None``) → snapped to the nearest
  element via ``mesh.find_closest_elems``, then the same ring search is
  applied. The result row carries ``snapped_to_mesh: true`` and a
  ``UserWarning`` is emitted. Co-located off-grid pairs follow the
  co-located path; separated pairs with one or both endpoints off-grid
  relocate each endpoint independently.

A pair therefore becomes ``status='separated'`` (true passthrough) only
when both endpoints are on-grid AND qualify with zero ring expansion.
Otherwise it becomes ``'relocated'`` (success) or ``'failed'``.

Failure handling
----------------
A *soft failure* means no element met the depth criteria within the search
radius, but at least two reachable elements were found. The pair is then
placed at the two deepest elements (sink → deepest, source → second-deepest)
and a min-depth dredge polygon is emitted for **each** location with
``attribute = max(0.5, mean_depth_at_that_element)``.

A *hard failure* is the degenerate case of only one reachable element. The
pair remains co-located there and a single polygon is emitted.

For separated pairs that fail (e.g. the swap-mesh scenario above), each
failing endpoint contributes its own min-depth polygon, so the user can
dredge each location independently. The pair-level severity is the worse
of the two endpoint severities (``hard`` if either is hard, otherwise
``soft``).

Outputs
-------
For prefix ``<P>`` and ``--out-dir <D>`` the run writes:

* ``<D>/<P>_source_relocated.yaml`` — full source set, drop-in replacement.
* ``<D>/<P>_sink_relocated.yaml`` — full sink set, drop-in replacement.
* ``<D>/<P>_relocation_report.yaml`` — *always written*; per-pair status,
  ring-hops used, final coordinates, and (for failures) severity and the
  mean depth at each chosen element. The canonical sentinel for failures is
  this file's ``failures`` list.
* ``<D>/depth_enforcement_<P>_source_sink.yaml`` — only when failures exist; one
  ``SchismPolygon`` per failed source / sink location, ready to splice into
  the project polygon yaml.
* ``<D>/<P>_skipped.yaml`` — only when there were name/pairing skips
  (entries that did not match ``(delta|suisun)_(src|sink)_<int>`` or had no
  matching partner with the same suffix).

Idempotence
-----------
Given the same mesh, parameters, and ``--seed``, the run is deterministic and
idempotent: pairs whose endpoints are already separated *and* on-grid *and*
qualify the depth criteria pass through unchanged, so a second run on the
same mesh is a no-op for all successful relocations. Reruns on a *different*
mesh may legitimately re-trigger relocation or surface new failures.

Dependency convention
---------------------
The library functions in this module deliberately avoid spatial-stack imports
(``geopandas``, ``fiona``, GDAL bindings) so that this module can be loaded
by lightweight preprocessor pipelines. Spatial imports are confined to the
optional shapefile-export branch of the CLI; see the maintainer note above
:func:`relocate_source_sink_cli`.

CLI
---
Registered as ``sch relocate_source_sink`` (and the standalone command
``relocate_source_sink``).

Author
------
Delta Modeling Section, DWR.
"""
from __future__ import annotations

import math
import os
import re
import warnings
from typing import Dict, Iterable, List, Optional, Set, Tuple

import click
import numpy as np
import pandas as pd
import yaml

from schimpy import schism_yaml
from schimpy.schism_mesh import read_mesh
from schimpy.schism_polygon import SchismPolygon, write_polygons


# -------------------------------------------------------------------------
# YAML I/O (auto-detects either the ``sources:``/``sinks:`` block layout
# read by :func:`schimpy.schism_sources_sinks.read_source_sink_yaml`, or the
# flat ``name: [x, y]`` mapping seen in ``source_suisun.yaml``/``sink_*.yaml``).
# -------------------------------------------------------------------------

_NAME_RE = re.compile(r"^(delta|suisun)_(src|sink)_(\d+)$")


def _load_yaml(fpath: str) -> dict:
    """Load a YAML file into a dictionary, tolerating empty files.

    Parameters
    ----------
    fpath : str
        Path to a YAML file.

    Returns
    -------
    dict
        Parsed mapping. An empty file returns ``{}``.

    Raises
    ------
    ValueError
        If the top-level YAML node is not a mapping.
    """
    with open(fpath, "r") as fh:
        data = yaml.safe_load(fh)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"{fpath}: expected a mapping at top level")
    return data


def read_points_yaml(fpath: str, kind: str) -> pd.DataFrame:
    """Read a points YAML in either the block or flat layout.

    Parameters
    ----------
    fpath : str
        Path to the YAML file.
    kind : {"sources", "sinks"}
        Which block to extract when the YAML uses the block layout. Ignored
        for the flat layout.

    Returns
    -------
    pandas.DataFrame
        Indexed by point name with float columns ``x`` and ``y``.
    """
    if kind not in ("sources", "sinks"):
        raise ValueError("kind must be 'sources' or 'sinks'")
    data = _load_yaml(fpath)
    if "sources" in data or "sinks" in data:
        block = data.get(kind, {}) or {}
    else:
        block = data
    rows = {}
    for name, xy in block.items():
        if xy is None or len(xy) < 2:
            warnings.warn(f"{fpath}: skipping malformed entry {name!r}")
            continue
        rows[str(name)] = (float(xy[0]), float(xy[1]))
    return pd.DataFrame.from_dict(rows, orient="index", columns=["x", "y"])


def write_points_yaml(
    fpath: str,
    df: pd.DataFrame,
    header: str = "",
) -> None:
    """Write a flat ``name: [x, y]`` YAML mirroring the input convention.

    Parameters
    ----------
    fpath : str
        Output YAML path.
    df : pandas.DataFrame
        Points indexed by name with float columns ``x`` and ``y``. Coordinates
        are rounded to 4 decimals on write.
    header : str, default ""
        Optional comment line(s) to prepend to the file verbatim. Each line
        should begin with ``#``.
    """
    payload = {name: [round(float(r.x), 4), round(float(r.y), 4)]
               for name, r in df.iterrows()}
    with open(fpath, "w") as fh:
        if header:
            fh.write(header.rstrip("\n") + "\n")
        yaml.safe_dump(payload, fh, default_flow_style=False, sort_keys=False)


# -------------------------------------------------------------------------
# Pairing
# -------------------------------------------------------------------------

def _parse_name(name: str) -> Optional[Tuple[str, str, int]]:
    """Parse a point name into ``(prefix, role, suffix)``.

    Parameters
    ----------
    name : str
        Point name such as ``"delta_src_5"`` or ``"suisun_sink_227"``.

    Returns
    -------
    tuple of (str, str, int) or None
        ``(prefix, role, suffix)`` if the name matches
        ``^(delta|suisun)_(src|sink)_(\\d+)$``; ``None`` otherwise.
    """
    m = _NAME_RE.match(name)
    if m is None:
        return None
    return m.group(1), m.group(2), int(m.group(3))


def build_pairs(
    df_source: pd.DataFrame,
    df_sink: pd.DataFrame,
    prefix: str,
) -> Tuple[List[Tuple[int, str, str]], Dict[str, str]]:
    """Pair sources and sinks by integer suffix for a given ``prefix``.

    Parameters
    ----------
    df_source : pandas.DataFrame
        Source points indexed by name (columns ``x``, ``y``).
    df_sink : pandas.DataFrame
        Sink points indexed by name (columns ``x``, ``y``).
    prefix : {"delta", "suisun"}
        Which family of names to keep. Names from the other family or with
        unrecognized layouts are reported in ``skipped``.

    Returns
    -------
    pairs : list of tuple
        ``(suffix, source_name, sink_name)`` for every name pair that shares
        the prefix and an integer suffix.
    skipped : dict
        Mapping ``name -> reason`` for entries excluded from pairing
        (bad name, wrong prefix/role, or no partner with the same suffix).
    """
    pairs: List[Tuple[int, str, str]] = []
    skipped: Dict[str, str] = {}

    src_by_n: Dict[int, str] = {}
    for name in df_source.index:
        parsed = _parse_name(name)
        if parsed is None:
            skipped[name] = "name does not match (delta|suisun)_src_<int>"
            continue
        pfx, role, n = parsed
        if pfx != prefix or role != "src":
            skipped[name] = f"prefix/role mismatch (expected {prefix}_src_*)"
            continue
        src_by_n[n] = name

    sink_by_n: Dict[int, str] = {}
    for name in df_sink.index:
        parsed = _parse_name(name)
        if parsed is None:
            skipped[name] = "name does not match (delta|suisun)_sink_<int>"
            continue
        pfx, role, n = parsed
        if pfx != prefix or role != "sink":
            skipped[name] = f"prefix/role mismatch (expected {prefix}_sink_*)"
            continue
        sink_by_n[n] = name

    common = sorted(set(src_by_n) & set(sink_by_n))
    for n in common:
        pairs.append((n, src_by_n[n], sink_by_n[n]))

    for n in sorted(set(src_by_n) - set(sink_by_n)):
        skipped[src_by_n[n]] = "no matching sink with same suffix"
    for n in sorted(set(sink_by_n) - set(src_by_n)):
        skipped[sink_by_n[n]] = "no matching source with same suffix"

    return pairs, skipped


# -------------------------------------------------------------------------
# Mesh helpers
# -------------------------------------------------------------------------

def _element_depth_stats(mesh, elem_i: int) -> Tuple[float, float]:
    """Return depth statistics for an element.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
        Mesh whose ``nodes`` array uses positive-down depth (deep > 0).
    elem_i : int
        Element index.

    Returns
    -------
    tuple of float
        ``(min_node_depth, mean_node_depth)``.
    """
    nodes = mesh.elem(elem_i)
    z = mesh.nodes[nodes, 2]
    return float(np.min(z)), float(np.mean(z))


def _element_edge_length(mesh, elem_i: int) -> float:
    """Mean edge length of an element in mesh units (typically meters).

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
    elem_i : int
        Element index.

    Returns
    -------
    float
        Arithmetic mean of the element's edge lengths.
    """
    nodes = mesh.elem(elem_i)
    xy = mesh.nodes[nodes, :2]
    n = len(nodes)
    lengths = [
        math.hypot(*(xy[(i + 1) % n] - xy[i])) for i in range(n)
    ]
    return float(np.mean(lengths))


def _seed_element(mesh, x: float, y: float) -> Tuple[int, bool]:
    """Locate the element containing ``(x, y)``, or the closest element.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
    x, y : float
        Query coordinates in mesh units.

    Returns
    -------
    elem_i : int
        Element index.
    snapped : bool
        ``True`` if the input point was outside the mesh and we fell back
        to :meth:`find_closest_elems`; ``False`` if the point was located
        inside an element.
    """
    elem = mesh.find_elem((x, y))
    if elem is not None:
        return int(elem), False
    elem = mesh.find_closest_elems((x, y), count=1)
    if isinstance(elem, list):
        elem = elem[0]
    return int(elem), True


def _ring_expand(mesh, seed: int, n_rings: int) -> List[Set[int]]:
    """BFS through element adjacency (via shared nodes) up to ``n_rings`` hops.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
    seed : int
        Seed element index.
    n_rings : int
        Maximum number of ring hops to expand.

    Returns
    -------
    list of set of int
        ``rings`` of length ``n_rings + 1`` (or shorter on early
        termination). ``rings[0] == {seed}`` and ``rings[k]`` is the set of
        elements first reached at hop ``k``.
    """
    visited: Set[int] = {seed}
    rings: List[Set[int]] = [{seed}]
    frontier: Set[int] = {seed}
    for _ in range(n_rings):
        nxt: Set[int] = set()
        for e in frontier:
            for node_i in mesh.elem(e):
                for neigh in mesh.get_elems_i_from_node(int(node_i)):
                    neigh = int(neigh)
                    if neigh not in visited:
                        nxt.add(neigh)
        if not nxt:
            rings.append(set())
            break
        visited |= nxt
        rings.append(nxt)
        frontier = nxt
    return rings


def _candidate_elements(rings: List[Set[int]]) -> List[int]:
    """Flatten BFS rings into a single ordered list (closest first).

    Parameters
    ----------
    rings : list of set of int
        Output of :func:`_ring_expand`.

    Returns
    -------
    list of int
        Element indices ordered by ring distance, then by element id within
        each ring for determinism.
    """
    out: List[int] = []
    for ring in rings:
        out.extend(sorted(ring))
    return out


def _element_centroid(mesh, elem_i: int) -> Tuple[float, float]:
    """Planar centroid of an element.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
    elem_i : int
        Element index.

    Returns
    -------
    tuple of float
        ``(x, y)`` in mesh units.
    """
    nodes = mesh.elem(elem_i)
    xy = mesh.nodes[nodes, :2]
    return float(xy[:, 0].mean()), float(xy[:, 1].mean())


# -------------------------------------------------------------------------
# Relocation core
# -------------------------------------------------------------------------

class RelocationFailure(Exception):
    """Raised internally when a pair fails the depth criteria within the
    requested ring radius."""


def _filter_qualifying(
    mesh,
    candidates: Iterable[int],
    min_node: float,
    mean_thresh: float,
    origin: Optional[Tuple[float, float]] = None,
    max_distance: Optional[float] = None,
) -> List[Tuple[int, float]]:
    """Filter candidate elements by depth criteria and optional distance cap.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
    candidates : iterable of int
        Element indices to test.
    min_node : float
        Required minimum nodal depth across all nodes of an element
        (positive-down).
    mean_thresh : float
        Required mean nodal depth.
    origin : tuple of float, optional
        ``(x, y)`` reference point for the optional distance cap.
    max_distance : float, optional
        If supplied with ``origin``, candidate elements whose centroid is
        farther than ``max_distance`` from ``origin`` are rejected.

    Returns
    -------
    list of tuple
        ``(elem_id, mean_depth)`` for every qualifying element, in the
        order they were visited (preserves caller-side determinism).
    """
    out = []
    for e in candidates:
        zmin, zmean = _element_depth_stats(mesh, e)
        if zmin < min_node or zmean < mean_thresh:
            continue
        if max_distance is not None and origin is not None:
            cx, cy = _element_centroid(mesh, e)
            if math.hypot(cx - origin[0], cy - origin[1]) > max_distance:
                continue
        out.append((e, zmean))
    return out


def _all_with_depth(
    mesh, candidates: Iterable[int]
) -> List[Tuple[int, float, float]]:
    """Return depth statistics for every candidate element.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
    candidates : iterable of int
        Element indices to score.

    Returns
    -------
    list of tuple
        ``(elem_id, min_depth, mean_depth)`` for each candidate.
    """
    return [(e, *_element_depth_stats(mesh, e)) for e in candidates]


def _assign_roles(
    mean_a: float,
    mean_b: float,
    rng: np.random.Generator,
    random_threshold: float,
) -> Tuple[int, int]:
    """Decide which of two candidate means becomes (source, sink).

    Parameters
    ----------
    mean_a, mean_b : float
        Mean nodal depths of the two candidate elements (positive-down).
    rng : numpy.random.Generator
        RNG used only when both candidates are deep enough to randomize.
    random_threshold : float
        Mean-depth threshold above which the deeper-to-sink rule is
        replaced by random assignment.

    Returns
    -------
    tuple of int
        ``(idx_for_source, idx_for_sink)`` where each idx is 0 or 1
        referring to ``(mean_a, mean_b)``.

    Notes
    -----
    Default rule: deeper -> sink, shallower -> source. When the deeper of
    the two has mean depth >= ``random_threshold`` the assignment is
    randomized so reruns do not perpetually bias one direction.
    """
    deeper = max(mean_a, mean_b)
    if deeper >= random_threshold:
        if rng.random() < 0.5:
            return 0, 1
        return 1, 0
    # deeper -> sink
    if mean_a >= mean_b:
        return 1, 0  # source=b, sink=a
    return 0, 1  # source=a, sink=b


# -------------------------------------------------------------------------
# Polygon helpers
# -------------------------------------------------------------------------

def _build_min_depth_polygon(
    mesh,
    elem_i: int,
    pair_name: str,
    deepest_mean: float,
    floor: float = 0.5,
) -> SchismPolygon:
    """Build a min-depth dredge polygon around an element.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
    elem_i : int
        Element to enclose.
    pair_name : str
        Used to derive the polygon ``name`` as ``f"{pair_name}_mindepth"``.
    deepest_mean : float
        Mean depth observed at the chosen element; floors the polygon
        attribute.
    floor : float, default 0.5
        Lower bound on the polygon attribute (m).

    Returns
    -------
    schimpy.schism_polygon.SchismPolygon
        Convex hull of the element's nodes buffered by half the local mean
        edge length, with ``type='min'`` and ``attribute`` set to
        ``max(floor, deepest_mean)`` formatted to 3 decimals.
    """
    from shapely.geometry import Polygon as _Polygon
    from shapely.geometry.polygon import orient as _orient

    nodes = mesh.elem(elem_i)
    xy = mesh.nodes[nodes, :2]
    poly = _Polygon([(float(p[0]), float(p[1])) for p in xy])
    poly = poly.convex_hull
    buf = 0.5 * _element_edge_length(mesh, elem_i)
    poly = poly.buffer(buf, quad_segs=8)
    poly = _orient(poly)
    attr = max(floor, float(deepest_mean))
    prop = {
        "name": f"{pair_name}_mindepth",
        "type": "min",
        "attribute": f"{attr:.3f}",
    }
    return SchismPolygon(shell=list(poly.exterior.coords), prop=prop)


def _relocate_single_point(
    mesh,
    x: float,
    y: float,
    rings_initial: int,
    rings_fallback: int,
    min_node_depth: float,
    mean_depth: float,
    max_distance: Optional[float],
) -> dict:
    """Project a single point onto a qualifying mesh element.

    Used for ``separated`` pairs whose individual endpoints lie outside the
    mesh (or on a dry/shallow element). Performs the same seed + ring
    expansion + depth filter as :func:`relocate_pair`, but for a lone
    point with no role-pairing constraint.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
    x, y : float
        Input coordinates (may be off-grid).
    rings_initial, rings_fallback : int
        BFS hop radii.
    min_node_depth, mean_depth : float
        Depth criteria.
    max_distance : float or None
        Distance cap measured from the seed centroid.

    Returns
    -------
    dict
        Keys: ``status`` ('relocated'|'failed'), ``x``, ``y``, ``elem``,
        ``snapped``, ``rings_used``. On ``'failed'`` also: ``severity``
        ('soft'|'hard') and ``mean`` (mean depth of the chosen element).
    """
    seed, snapped = _seed_element(mesh, x, y)
    seed_min, seed_mean = _element_depth_stats(mesh, seed)
    origin = _element_centroid(mesh, seed) if snapped else (x, y)

    if seed_min >= min_node_depth and seed_mean >= mean_depth:
        cx, cy = _element_centroid(mesh, seed)
        return {
            "status": "relocated",
            "x": cx, "y": cy,
            "elem": int(seed),
            "snapped": snapped,
            "rings_used": 0,
        }

    cands_last: List[int] = [seed]
    for n_rings in (rings_initial, rings_fallback):
        rings = _ring_expand(mesh, seed, n_rings)
        cands_last = _candidate_elements(rings)
        qualifying = _filter_qualifying(
            mesh, cands_last, min_node_depth, mean_depth,
            origin=origin, max_distance=max_distance,
        )
        if qualifying:
            # Pick deepest qualifying element (deterministic).
            qualifying.sort(key=lambda t: t[1], reverse=True)
            e = qualifying[0][0]
            cx, cy = _element_centroid(mesh, e)
            return {
                "status": "relocated",
                "x": cx, "y": cy,
                "elem": int(e),
                "snapped": snapped,
                "rings_used": n_rings,
            }

    # Failure: fall back to deepest available element in the explored set.
    scored = _all_with_depth(mesh, cands_last)
    if scored:
        scored.sort(key=lambda t: t[2], reverse=True)
        e, _zmin, zmean = scored[0]
        cx, cy = _element_centroid(mesh, e)
        return {
            "status": "failed",
            "severity": "soft",
            "x": cx, "y": cy,
            "elem": int(e),
            "snapped": snapped,
            "rings_used": rings_fallback,
            "mean": float(zmean),
        }
    cx, cy = _element_centroid(mesh, seed)
    return {
        "status": "failed",
        "severity": "hard",
        "x": cx, "y": cy,
        "elem": int(seed),
        "snapped": snapped,
        "rings_used": rings_fallback,
        "mean": float(seed_mean),
    }


# -------------------------------------------------------------------------
# Main driver
# -------------------------------------------------------------------------

def relocate_pair(
    mesh,
    xs: float,
    ys: float,
    xn: float,
    yn: float,
    rng: np.random.Generator,
    rings_initial: int = 3,
    rings_fallback: int = 10,
    min_node_depth: float = 0.0,
    mean_depth: float = 1.0,
    random_threshold: float = 2.0,
    separation_tol: float = 0.01,
    always_move_both: bool = True,
    max_distance: Optional[float] = 100.0,
) -> dict:
    """Compute the relocated coordinates for a single source/sink pair.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
        Mesh with positive-down depths in ``mesh.nodes[:, 2]``.
    xs, ys : float
        Original source coordinates (mesh units).
    xn, yn : float
        Original sink coordinates (mesh units).
    rng : numpy.random.Generator
        Source of randomness for shuffling qualifying candidates and for the
        deep-water role-assignment tie-breaker.
    rings_initial : int, default 3
        Initial element-ring search radius (BFS hops).
    rings_fallback : int, default 10
        Expanded search radius when the initial radius yields no qualifying
        elements.
    min_node_depth : float, default 0.0
        Required minimum nodal depth at every node of the chosen element.
    mean_depth : float, default 1.0
        Required mean nodal depth at the chosen element.
    random_threshold : float, default 2.0
        See :func:`_assign_roles`.
    separation_tol : float, default 0.01
        Pairs with planar distance above this threshold (m) are considered
        already separated and pass through unchanged.
    always_move_both : bool, default True
        If True, both points relocate to two distinct elements. If False and
        the seed element qualifies, the source stays at the seed and only
        the sink moves.
    max_distance : float or None, default 100.0
        Physical cap (mesh units) on how far a candidate centroid may lie
        from the original co-location point. ``None`` disables the cap.

    Returns
    -------
    dict
        Common keys:

        * ``status`` : {'separated', 'relocated', 'failed'}
        * ``source`` : tuple of float -- final ``(x, y)``
        * ``sink`` : tuple of float -- final ``(x, y)``
        * ``source_elem``, ``sink_elem`` : int or None
        * ``rings_used`` : int

        On ``'failed'`` also:

        * ``severity`` : {'soft', 'hard'} -- 'hard' indicates only one
          reachable element so source and sink remain co-located
        * ``source_mean``, ``sink_mean`` : float -- mean depths at the two
          chosen elements

        Always present (except for ``'separated'`` where it is ``False``):

        * ``snapped`` : bool -- ``True`` if the original co-location point
          lay outside the mesh and was snapped to the nearest element
          before searching.
    """
    if math.hypot(xs - xn, ys - yn) > separation_tol:
        # Pair is already separated. Validate each endpoint independently:
        # if it sits in a qualifying mesh element it passes through unchanged;
        # otherwise (off-grid OR on-grid but too shallow) it is projected
        # onto the nearest qualifying element via _relocate_single_point.
        # This matters when a previously-validated source/sink set is being
        # rerun on a swapped mesh: prior separation does not imply prior
        # depth qualification.
        def _resolve(x, y):
            elem = mesh.find_elem((x, y))
            if elem is not None:
                zmin, zmean = _element_depth_stats(mesh, int(elem))
                if zmin >= min_node_depth and zmean >= mean_depth:
                    return {
                        "kind": "ok",
                        "x": x, "y": y,
                        "elem": int(elem),
                        "snapped": False,
                        "rings_used": 0,
                    }
            return _relocate_single_point(
                mesh, x, y,
                rings_initial=rings_initial,
                rings_fallback=rings_fallback,
                min_node_depth=min_node_depth,
                mean_depth=mean_depth,
                max_distance=max_distance,
            )

        src_res = _resolve(xs, ys)
        sink_res = _resolve(xn, yn)
        # Treat the {"kind": "ok"} branch as a benign relocation that keeps
        # input coordinates and is *not* counted as snapped.
        src_status = src_res.get("status", "ok")
        sink_status = sink_res.get("status", "ok")
        src_snapped = bool(src_res.get("snapped", False))
        sink_snapped = bool(sink_res.get("snapped", False))

        # If both endpoints validated cleanly with no relocation, this is a
        # true passthrough.
        if (src_status == "ok" and sink_status == "ok"
                and src_res["rings_used"] == 0
                and sink_res["rings_used"] == 0):
            return {
                "status": "separated",
                "source": (xs, ys),
                "sink": (xn, yn),
                "source_elem": None,
                "sink_elem": None,
                "rings_used": 0,
                "snapped": False,
            }

        if src_status == "failed" or sink_status == "failed":
            # Take the worse of the two severities. ``hard`` means a
            # degenerate single-element search; ``soft`` means depth was
            # insufficient at the chosen element.
            severities = [
                s for s in (src_res.get("severity"), sink_res.get("severity"))
                if s is not None
            ]
            severity = "hard" if "hard" in severities else "soft"
            src_mean = src_res.get(
                "mean",
                _element_depth_stats(mesh, src_res["elem"])[1],
            )
            sink_mean = sink_res.get(
                "mean",
                _element_depth_stats(mesh, sink_res["elem"])[1],
            )
            return {
                "status": "failed",
                "severity": severity,
                "source": (src_res["x"], src_res["y"]),
                "sink": (sink_res["x"], sink_res["y"]),
                "source_elem": int(src_res["elem"]),
                "sink_elem": int(sink_res["elem"]),
                "rings_used": rings_fallback,
                "source_mean": float(src_mean),
                "sink_mean": float(sink_mean),
                "snapped": bool(src_snapped or sink_snapped),
            }

        return {
            "status": "relocated",
            "source": (src_res["x"], src_res["y"]),
            "sink": (sink_res["x"], sink_res["y"]),
            "source_elem": int(src_res["elem"]),
            "sink_elem": int(sink_res["elem"]),
            "rings_used": max(src_res["rings_used"], sink_res["rings_used"]),
            "snapped": bool(src_snapped or sink_snapped),
        }

    # Co-located: pick a single seed at the midpoint.
    x0, y0 = 0.5 * (xs + xn), 0.5 * (ys + yn)
    seed, snapped = _seed_element(mesh, x0, y0)
    seed_min, seed_mean = _element_depth_stats(mesh, seed)
    seed_qualifies = (seed_min >= min_node_depth) and (seed_mean >= mean_depth)

    # When the original point lay outside the mesh, the physical-distance
    # cap must be measured from the snapped seed centroid, not from the
    # off-grid input.
    if snapped:
        origin = _element_centroid(mesh, seed)
    else:
        origin = (x0, y0)

    for n_rings in (rings_initial, rings_fallback):
        rings = _ring_expand(mesh, seed, n_rings)
        cands = _candidate_elements(rings)
        qualifying = _filter_qualifying(
            mesh, cands, min_node_depth, mean_depth,
            origin=origin, max_distance=max_distance,
        )
        if not qualifying:
            continue

        # Decide whether seed is reused (minimal-move) or both move.
        move_both = always_move_both or (not seed_qualifies)

        # Sort qualifying with deterministic random pick.
        rng.shuffle(qualifying)

        # Source and sink must end up in *different* elements. Partition the
        # qualifying set into "non-seed" and "seed" so we can guarantee that.
        non_seed = [(e, m) for (e, m) in qualifying if e != seed]
        has_seed_qual = any(e == seed for (e, _m) in qualifying)

        if move_both:
            if len(non_seed) >= 2:
                cand_pair = [non_seed[0], non_seed[1]]
            elif len(non_seed) == 1 and seed_qualifies:
                # one non-seed element + the seed itself
                cand_pair = [(seed, seed_mean), non_seed[0]]
            else:
                # cannot satisfy two distinct qualifying elements at this
                # ring radius; expand search
                continue
        else:
            # Keep one at the seed; move the partner to the best non-seed.
            if not non_seed:
                continue
            e_q, m_q = non_seed[0]
            cand_pair = [(seed, seed_mean), (e_q, m_q)]

        # Final invariant check.
        assert cand_pair[0][0] != cand_pair[1][0], (
            "internal error: source and sink resolved to the same element"
        )

        idx_src, idx_sink = _assign_roles(
            cand_pair[0][1], cand_pair[1][1], rng, random_threshold
        )
        e_src = cand_pair[idx_src][0]
        e_sink = cand_pair[idx_sink][0]
        x_src, y_src = _element_centroid(mesh, e_src)
        x_sink, y_sink = _element_centroid(mesh, e_sink)
        return {
            "status": "relocated",
            "source": (x_src, y_src),
            "sink": (x_sink, y_sink),
            "source_elem": int(e_src),
            "sink_elem": int(e_sink),
            "rings_used": n_rings,
            "snapped": snapped,
        }

    # Failure: both ring radii exhausted. Pick the two deepest distinct
    # elements found in the search so source and sink are still separated.
    rings = _ring_expand(mesh, seed, rings_fallback)
    cands = _candidate_elements(rings)
    scored = _all_with_depth(mesh, cands)
    # rank by mean depth descending
    scored.sort(key=lambda t: t[2], reverse=True)
    if len(scored) >= 2:
        sink_elem, _smin, sink_mean = scored[0]
        src_elem, _smin2, src_mean = scored[1]
        severity = "soft"   # placed in best-available water; depth too shallow
    else:
        # Degenerate: only a single reachable element. Source and sink remain
        # co-located in this case — this is a *hard* failure.
        sink_elem, _smin, sink_mean = scored[0]
        src_elem, src_mean = sink_elem, sink_mean
        severity = "hard"
    x_sink, y_sink = _element_centroid(mesh, sink_elem)
    x_src, y_src = _element_centroid(mesh, src_elem)
    return {
        "status": "failed",
        "severity": severity,
        "source": (x_src, y_src),
        "sink": (x_sink, y_sink),
        "source_elem": int(src_elem),
        "sink_elem": int(sink_elem),
        "rings_used": rings_fallback,
        "source_mean": float(src_mean),
        "sink_mean": float(sink_mean),
        "snapped": snapped,
    }


def relocate_source_sink(
    hgrid: str,
    source_yaml: str,
    sink_yaml: str,
    prefix: str,
    out_dir: str = ".",
    rings_initial: int = 3,
    rings_fallback: int = 10,
    min_node_depth: float = 0.0,
    mean_depth: float = 1.0,
    random_threshold: float = 2.0,
    separation_tol: float = 0.01,
    always_move_both: bool = True,
    seed: int = 42,
    max_distance: Optional[float] = 100.0,
    out_suffix: str = "_adj",
) -> dict:
    """End-to-end driver: read inputs, relocate co-located pairs, write outputs.

    Parameters
    ----------
    hgrid : str
        Path to the SCHISM ``hgrid.gr3`` mesh (positive-down depth).
    source_yaml, sink_yaml : str
        Paths to the paired source and sink YAMLs. Either the
        ``sources:``/``sinks:`` block layout or a flat ``name: [x, y]``
        mapping is accepted.
    prefix : {"delta", "suisun"}
        Family of names processed in this run; the other family is sent to
        the skipped report.
    out_dir : str, default "."
        Output directory; created if missing.
    rings_initial, rings_fallback : int
        Element-ring search radii (BFS hops). See :func:`relocate_pair`.
    min_node_depth, mean_depth : float
        Depth criteria. See :func:`relocate_pair`.
    random_threshold : float
        Mean-depth threshold above which role assignment is randomized.
    separation_tol : float
        Distance threshold (m) above which a pair is treated as already
        separated and passes through unchanged.
    always_move_both : bool
        Whether to always move both points to distinct elements.
    seed : int, default 42
        RNG seed (``numpy.random.default_rng``) for repeatable runs.
    max_distance : float or None, default 100.0
        Physical cap (m) on relocation distance.
    out_suffix : str, default "_adj"
        Suffix appended to ``prefix`` in the output filenames, e.g.
        ``delta_source_adj.yaml``. Change to ``"_relocated"`` to keep the
        previous naming convention.

    Returns
    -------
    dict
        Result with the following keys:

        * ``counts`` -- dict with keys ``separated``, ``relocated``,
          ``failed``, ``failed_soft``, ``failed_hard``
        * ``n_pairs``, ``n_skipped`` -- ints
        * ``failures`` -- list of ``f"{prefix}_pair_{n}"`` strings
        * ``source_yaml``, ``sink_yaml``, ``report_yaml`` -- always-present
          output paths
        * ``polygon_yaml`` -- path or ``None``
        * ``skipped_yaml`` -- path or ``None``
        * ``summary`` -- ``pandas.DataFrame`` with one row per processed
          pair (status, severity, depths, final coordinates)

    Notes
    -----
    The canonical sentinel for failures is ``result['failures']`` (or
    ``result['counts']['failed'] > 0``). The
    ``<prefix>_relocation_report.yaml`` file is always written and contains
    the same per-pair summary.
    """
    if prefix not in ("delta", "suisun"):
        raise ValueError("prefix must be 'delta' or 'suisun'")
    os.makedirs(out_dir, exist_ok=True)

    mesh = read_mesh(hgrid)
    df_source = read_points_yaml(source_yaml, "sources")
    df_sink = read_points_yaml(sink_yaml, "sinks")

    pairs, skipped = build_pairs(df_source, df_sink, prefix)
    rng = np.random.default_rng(seed)

    out_source = df_source.copy()
    out_sink = df_sink.copy()
    polygons: List[SchismPolygon] = []
    counts = {"separated": 0, "relocated": 0, "failed": 0,
              "failed_soft": 0, "failed_hard": 0}
    failure_names: List[str] = []
    summary_rows = []

    for n, sname, kname in pairs:
        xs, ys = float(df_source.loc[sname, "x"]), float(df_source.loc[sname, "y"])
        xn, yn = float(df_sink.loc[kname, "x"]), float(df_sink.loc[kname, "y"])
        result = relocate_pair(
            mesh, xs, ys, xn, yn, rng,
            rings_initial=rings_initial,
            rings_fallback=rings_fallback,
            min_node_depth=min_node_depth,
            mean_depth=mean_depth,
            random_threshold=random_threshold,
            separation_tol=separation_tol,
            always_move_both=always_move_both,
            max_distance=max_distance,
        )
        counts[result["status"]] += 1
        out_source.loc[sname, ["x", "y"]] = result["source"]
        out_sink.loc[kname, ["x", "y"]] = result["sink"]
        row = {
            "suffix": n,
            "source": sname,
            "sink": kname,
            "status": result["status"],
            "rings": result["rings_used"],
            "source_xy": [float(result["source"][0]), float(result["source"][1])],
            "sink_xy": [float(result["sink"][0]), float(result["sink"][1])],
        }
        if result.get("snapped"):
            row["snapped_to_mesh"] = True
            warnings.warn(
                f"{sname}/{kname}: original location lay outside the mesh; "
                f"snapped to nearest element before relocation"
            )
        if result["status"] == "failed":
            row["severity"] = result["severity"]
            row["source_mean_depth"] = result["source_mean"]
            row["sink_mean_depth"] = result["sink_mean"]
            counts["failed_" + result["severity"]] += 1
        summary_rows.append(row)
        if result["status"] == "failed":
            failure_names.append(f"{prefix}_pair_{n}")
            # Emit one min-depth polygon per point so each location can be
            # dredged independently. Same element -> one polygon (degenerate
            # single-element search).
            poly_src = _build_min_depth_polygon(
                mesh,
                result["source_elem"],
                pair_name=sname,
                deepest_mean=result["source_mean"],
            )
            polygons.append(poly_src)
            if result["sink_elem"] != result["source_elem"]:
                poly_sink = _build_min_depth_polygon(
                    mesh,
                    result["sink_elem"],
                    pair_name=kname,
                    deepest_mean=result["sink_mean"],
                )
                polygons.append(poly_sink)

    _HEADER = (
        "# Adjusted using relocate_source_sink to separate source sink pairs "
        "and ensure placement at adequate depth."
    )

    # Derive output filenames from the input basenames so that e.g.
    # source_dcd.yaml → source_dcd_adj.yaml regardless of prefix.
    def _out_path(in_yaml: str) -> str:
        base, ext = os.path.splitext(os.path.basename(in_yaml))
        return os.path.join(out_dir, f"{base}{out_suffix}{ext}")

    src_out = _out_path(source_yaml)
    sink_out = _out_path(sink_yaml)
    write_points_yaml(src_out, out_source, header=_HEADER)
    write_points_yaml(sink_out, out_sink, header=_HEADER)

    poly_out = None
    if polygons:
        poly_out = os.path.join(out_dir, f"depth_enforcement_{prefix}_source_sink.yaml")
        write_polygons(poly_out, polygons)

    # Always write the per-pair report (sentinel for failures).
    report_out = os.path.join(out_dir, f"{prefix}_relocation_report.yaml")
    with open(report_out, "w") as fh:
        yaml.safe_dump(
            {
                "counts": counts,
                "n_pairs": len(pairs),
                "failures": failure_names,
                "pairs": summary_rows,
            },
            fh,
            default_flow_style=False,
            sort_keys=False,
        )

    skipped_out = None
    if skipped:
        skipped_out = os.path.join(out_dir, f"{prefix}_skipped.yaml")
        with open(skipped_out, "w") as fh:
            yaml.safe_dump(
                {"skipped": skipped},
                fh,
                default_flow_style=False,
                sort_keys=True,
            )
        for name, why in skipped.items():
            warnings.warn(f"{name}: {why}")

    if failure_names:
        soft = counts["failed_soft"]
        hard = counts["failed_hard"]
        warnings.warn(
            f"{len(failure_names)} pair(s) failed depth criteria "
            f"({soft} soft, {hard} hard); see {report_out}"
            + (f" and {poly_out}" if poly_out else "")
        )

    return {
        "counts": counts,
        "n_pairs": len(pairs),
        "n_skipped": len(skipped),
        "failures": failure_names,
        "source_yaml": src_out,
        "sink_yaml": sink_out,
        "polygon_yaml": poly_out,
        "skipped_yaml": skipped_out,
        "report_yaml": report_out,
        "summary": pd.DataFrame(summary_rows),
    }


# -------------------------------------------------------------------------
# CLI
#
# NOTE TO MAINTAINERS:
#   The spatial stack (geopandas, fiona, GDAL bindings) is heavy and pulls in
#   a long dependency chain. Keep its imports inside this CLI function so
#   that the library callers above remain importable in lightweight
#   environments. Do *not* hoist `import geopandas` etc. to module scope.
# -------------------------------------------------------------------------

@click.command(
    help=(
        "Relocate co-located SCHISM source/sink point pairs. Pairs are "
        "matched by integer suffix within a single prefix ('delta' or "
        "'suisun'); pairs already separated beyond the tolerance pass "
        "through unchanged. Failure cases are recorded in a min-depth "
        "polygon YAML for manual splicing into the project polygon file."
    )
)
@click.option("--hgrid", required=True, type=click.Path(exists=True),
              help="Path to hgrid.gr3 (positive-down depth convention).")
@click.option("--source", "source_yaml", required=True,
              type=click.Path(exists=True),
              help="Input source YAML (e.g. source_suisun.yaml).")
@click.option("--sink", "sink_yaml", required=True,
              type=click.Path(exists=True),
              help="Input sink YAML (e.g. sink_suisun.yaml).")
@click.option("--prefix", required=True,
              type=click.Choice(["delta", "suisun"]),
              help="Name prefix for the family of pairs to process.")
@click.option("--out-dir", type=click.Path(), default=".",
              show_default=True,
              help="Directory for output YAML files.")
@click.option("--rings-initial", type=int, default=6, show_default=True,
              help="Initial BFS ring radius for the local search.")
@click.option("--rings-fallback", type=int, default=12, show_default=True,
              help="Expanded ring radius if the initial search fails.")
@click.option("--min-node-depth", type=float, default=0.0, show_default=True,
              help="Minimum depth (m) required at every node of the chosen "
                   "element (positive-down).")
@click.option("--mean-depth", type=float, default=0.5, show_default=True,
              help="Minimum mean nodal depth (m) required at the chosen "
                   "element.")
@click.option("--random-threshold", type=float, default=2.0,
              show_default=True,
              help="If both candidates have mean depth >= this value, the "
                   "deeper-to-sink rule is replaced by random assignment.")
@click.option("--separation-tol", type=float, default=0.01, show_default=True,
              help="Pairs with planar distance above this threshold (m) are "
                   "considered already separated.")
@click.option("--max-distance", type=float, default=500.0, show_default=True,
              help="Maximum physical distance (m) a relocated point may move "
                   "from the original co-location. Candidate elements whose "
                   "centroid exceeds this radius are rejected. Set to 0 or a "
                   "negative value to disable.")
@click.option("--seed", type=int, default=42, show_default=True,
              help="RNG seed for repeatable random rejection sampling.")
@click.option("--always-move-both/--minimal-move", default=True,
              show_default=True,
              help="When relocating, move both points to distinct elements "
                   "(default), or keep one at the seed if it qualifies.")
@click.option("--out-suffix", default="_adj", show_default=True,
              help="Suffix appended to the prefix in output filenames, e.g. "
                   "delta_source_adj.yaml. Use _relocated for the old naming.")
@click.option("--shp/--no-shp", default=False, show_default=True,
              help="Also export relocated points and polygons as shapefiles.")
@click.help_option("-h", "--help")
def relocate_source_sink_cli(
    hgrid, source_yaml, sink_yaml, prefix, out_dir,
    rings_initial, rings_fallback, min_node_depth, mean_depth,
    random_threshold, separation_tol, max_distance, seed,
    always_move_both, out_suffix, shp,
):
    """Click entry point for :func:`relocate_source_sink`.

    See the module docstring for the high-level workflow and output
    contract. This function is intentionally the *only* place where
    spatial-stack imports (``geopandas``, ``shapely.geometry``) appear at
    module load time, and even then only inside the ``--shp`` branch.
    """
    md = max_distance if max_distance is not None and max_distance > 0 else None
    result = relocate_source_sink(
        hgrid=hgrid,
        source_yaml=source_yaml,
        sink_yaml=sink_yaml,
        prefix=prefix,
        out_dir=out_dir,
        rings_initial=rings_initial,
        rings_fallback=rings_fallback,
        min_node_depth=min_node_depth,
        mean_depth=mean_depth,
        random_threshold=random_threshold,
        separation_tol=separation_tol,
        always_move_both=always_move_both,
        seed=seed,
        max_distance=md,
        out_suffix=out_suffix,
    )

    click.echo(
        f"Pairs: {result['n_pairs']}  "
        f"(separated={result['counts']['separated']}, "
        f"relocated={result['counts']['relocated']}, "
        f"failed={result['counts']['failed']} "
        f"[soft={result['counts']['failed_soft']}, "
        f"hard={result['counts']['failed_hard']}])"
    )
    click.echo(f"Wrote {result['source_yaml']}")
    click.echo(f"Wrote {result['sink_yaml']}")
    click.echo(f"Wrote {result['report_yaml']}")
    if result["polygon_yaml"]:
        click.echo(f"Wrote {result['polygon_yaml']}")
    if result["skipped_yaml"]:
        click.echo(f"Wrote {result['skipped_yaml']}")
    if result["failures"]:
        click.echo(
            "\nFAILED pairs (depth criteria not met within the search "
            "radius; points were placed at the deepest reachable elements "
            "and a min-depth polygon was emitted for each location):"
        )
        for row in result["summary"].itertuples(index=False):
            if row.status != "failed":
                continue
            click.echo(
                f"  {prefix}_pair_{row.suffix}  severity={row.severity}  "
                f"src_mean={row.source_mean_depth:.2f} m  "
                f"sink_mean={row.sink_mean_depth:.2f} m"
            )
        click.echo(
            "  (severity=soft: source/sink separated to two distinct "
            "elements but neither met the mean-depth threshold; "
            "severity=hard: only one reachable element \u2014 source and "
            "sink remain co-located)"
        )

    if shp:
        # Heavy spatial-stack imports kept local to this branch on purpose.
        # See module docstring and the maintainer note above.
        from schimpy.convert_points import points_to_shp  # noqa: F401
        import geopandas as gpd  # noqa: F401
        from shapely.geometry import Point

        for label, path in (
            ("source", result["source_yaml"]),
            ("sink", result["sink_yaml"]),
        ):
            df = read_points_yaml(path, label + "s")
            df = df.reset_index().rename(columns={"index": "sites"})
            df["stype"] = label
            shp_path = path.replace(".yaml", ".shp")
            from schimpy.convert_points import points_to_shp as _to_shp
            _to_shp(shp_path, df)
            click.echo(f"Wrote {shp_path}")
        if result["polygon_yaml"]:
            shp_path = result["polygon_yaml"].replace(".yaml", ".shp")
            from schimpy.schism_polygon import write_polygons as _wp
            from schimpy.schism_polygon import read_polygons as _rp
            _wp(shp_path, _rp(result["polygon_yaml"]))
            click.echo(f"Wrote {shp_path}")

    click.echo(
        "\nReminder: convert YAML outputs to shapefiles any time with:\n"
        "  sch convert_points  --input <file>.yaml --output <file>.shp\n"
        "  sch convert_polygons --input <file>.yaml --output <file>.shp"
    )


if __name__ == "__main__":
    relocate_source_sink_cli()
