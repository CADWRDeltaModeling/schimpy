#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import annotations
import logging


"""
mesh_volume_tvd
===============

TVD-regularized volume tuning workflow wrapped for two entry points:

1) **Library API** (`refine_volume_tvd`) for the SCHISM preprocessor
   so it can be driven from a partial YAML tree (`depth_optimization: {method:
   volume_tvd, ...}`).
2) **CLI** (`python -m schimpy.mesh_volume_tvd ...`) that mirrors what
   `mesh_refine_example.py` has been doing, so you can run the workflow
   outside the preprocessor.

This module deliberately gathers logic that previously lived in the example
script (shoreline detection, floor construction, TVD run) without changing the
core algorithm implementations in :mod:`mesh_refinement_tvd` or
:mod:`shore_edge`.

Notes
-----
- We **do not** mutate boundary objects; we only update ``mesh.nodes[:, 2]``.
  Boundary information in the SCHISM mesh object is preserved.
- DEM sampling uses :func:`schimpy.stacked_dem_fill.create_dem_sampler` and is
  compatible with the point-population path that
  :func:`schimpy.stacked_dem_fill.stacked_dem_fill` provides.
"""


from dataclasses import dataclass, field
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
    Callable,
    Set,
)

import os
import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from schimpy.stacked_dem_fill import create_dem_sampler
from schimpy.schism_mesh import *
import matplotlib.pyplot as plt
from collections import defaultdict, deque
import csv
import argparse
import logging

try:
    from shapely.geometry import Point, LineString, mapping
except Exception:
    raise SystemExit("This script requires shapely (pip install shapely).")

EdgeLabel = str  # 'deep' | 'shore'

# ------------------------------- utilities -----------------------------------


def _n_nodes(mesh) -> int:
    return mesh.n_nodes() if hasattr(mesh, "n_nodes") else int(mesh.nodes.shape[0])


def _n_elems(mesh) -> int:
    return (
        mesh.n_elems()
        if hasattr(mesh, "n_elems")
        else int(getattr(mesh, "elem", np.empty((0,))).shape[0])
    )


def _edge_mean_depth(mesh, node_depth: np.ndarray) -> np.ndarray:
    e = mesh.edges
    return 0.5 * (node_depth[e[:, 0].astype(int)] + node_depth[e[:, 1].astype(int)])


def _build_elem_edges(mesh) -> List[List[int]]:
    ne = _n_elems(mesh)
    elem_edges: List[List[int]] = [[] for _ in range(ne)]
    edges = mesh.edges  # [n1, n2, type, elem_a, elem_b]
    for ei in range(edges.shape[0]):
        ea = int(edges[ei, 3])
        eb = int(edges[ei, 4])
        if ea >= 0:
            elem_edges[ea].append(ei)
        if eb >= 0:
            elem_edges[eb].append(ei)
    return elem_edges


def _seed_elems_from_xy(mesh, seeds_xy: Iterable[Tuple[float, float]]) -> List[int]:
    if getattr(mesh, "_polygons", None) is None:
        mesh._create_2d_polygons()
    if getattr(mesh, "_centroids", None) is None:
        mesh._calculate_centroids()
    out: List[int] = []
    for x, y in seeds_xy:
        p = Point(x, y)
        hit = None
        for ei, poly in enumerate(mesh._polygons):
            if poly.contains(p) or poly.touches(p):
                hit = ei
                break
        if hit is None:
            cxy = mesh._centroids
            d2 = (cxy[:, 0] - x) ** 2 + (cxy[:, 1] - y) ** 2
            hit = int(np.argmin(d2))
        out.append(int(hit))
    return out


def _neighbor_across(mesh, elem_i: int, edge_i: int) -> int:
    ea, eb = int(mesh.edges[edge_i, 3]), int(mesh.edges[edge_i, 4])
    if ea == elem_i:
        return eb
    if eb == elem_i:
        return ea
    return -1


def _perimeter_edges_of_region(mesh, in_region: np.ndarray) -> Set[int]:
    perim: Set[int] = set()
    edges = mesh.edges
    for ei in range(edges.shape[0]):
        ea, eb = int(edges[ei, 3]), int(edges[ei, 4])
        a_in = (ea >= 0) and bool(in_region[ea])
        b_in = (eb >= 0) and bool(in_region[eb])
        if a_in ^ b_in or (a_in and eb == -1) or (b_in and ea == -1):
            perim.add(ei)
    return perim


def _build_perimeter_loops(mesh, perim_edges: Set[int]) -> List[Dict]:
    """Deterministic ordering of loops and edges."""

    edges = mesh.edges
    adj = defaultdict(list)
    for eidx in sorted(perim_edges):
        n1, n2 = int(edges[eidx, 0]), int(edges[eidx, 1])
        adj[n1].append((eidx, n2))
        adj[n2].append((eidx, n1))
    for n in adj:
        adj[n].sort(key=lambda t: t[0])

    unvisited = set(perim_edges)
    loops = []
    while unvisited:
        e0 = min(unvisited)
        unvisited.remove(e0)
        n1, n2 = int(edges[e0, 0]), int(edges[e0, 1])
        path_nodes = [n1, n2]
        path_edges = [e0]
        curr, prev = n2, n1
        while True:
            nxt = None
            for eidx, other in adj[curr]:
                if eidx in unvisited and other != prev:
                    nxt = (eidx, other)
                    break
            if nxt is None:
                for eidx, other in adj[curr]:
                    if eidx in unvisited:
                        nxt = (eidx, other)
                        break
            if nxt is None:
                break
            eidx, other = nxt
            unvisited.remove(eidx)
            path_edges.append(eidx)
            path_nodes.append(other)
            prev, curr = curr, other
            if other == path_nodes[0]:
                break
        if path_nodes[-1] != path_nodes[0]:
            path_nodes.append(path_nodes[0])
        loops.append(dict(nodes=path_nodes, edges=path_edges))
    loops.sort(key=lambda L: min(L["edges"]) if L["edges"] else 10**12)
    return loops


def _turn_angle(p_prev: np.ndarray, p_curr: np.ndarray, p_next: np.ndarray) -> float:
    v1 = p_curr - p_prev
    v2 = p_next - p_curr
    cross = v1[0] * v2[1] - v1[1] * v2[0]
    dot = v1[0] * v2[0] + v1[1] * v2[1]
    return np.arctan2(cross, dot)


# ------------------------------- h0 loading ----------------------------------


def _load_h0(h0_arg: str, mesh) -> np.ndarray:
    """
    Accepts either a float string or a .gr3 path. Returns h0_node array (len = n_nodes).
    """
    s = str(h0_arg).strip()
    # .gr3 file?
    if s.lower().endswith(".gr3"):
        h0_mesh = read_mesh(s)
        if _n_nodes(h0_mesh) != _n_nodes(mesh):
            raise ValueError(
                f"h0 .gr3 node count mismatch: { _n_nodes(h0_mesh) } vs { _n_nodes(mesh) }"
            )
        # depth (positive down), no sign change
        return h0_mesh.nodes[:, 2].astype(float)
    # else try float
    try:
        h0_scalar = float(s)
    except Exception as e:
        raise ValueError(f"--h0 must be a float or .gr3 path; got {h0_arg!r}") from e
    return np.full(_n_nodes(mesh), h0_scalar, dtype=float)


def _edge_field_from_nodes(mesh, nodal: np.ndarray) -> np.ndarray:
    e = mesh.edges[:, :2].astype(int)
    return 0.5 * (nodal[e[:, 0]] + nodal[e[:, 1]])


# ----------------------------- stage 1 & 2 -----------------------------------


@dataclass
class FloodFillResult:
    wet_elems_stage1: np.ndarray
    edge_labels_stage1: Dict[int, str]
    wet_elems_stage2: np.ndarray
    edge_labels_stage2: Dict[int, str]


def floodfill_always_wet(
    mesh,
    seeds_xy: Iterable[Tuple[float, float]],
    h0_node: np.ndarray,
    h0_edge: np.ndarray,
    delta0: float,
    delta1: float,
) -> FloodFillResult:
    node_depth = mesh.nodes[:, 2].astype(float)
    edge_mean_d = _edge_mean_depth(mesh, node_depth)

    # Edge eligibility (for growth and per-element qualification)
    edge_is_wet = edge_mean_d > (h0_edge - delta1)

    elem_edges = _build_elem_edges(mesh)
    nE = _n_elems(mesh)
    qualifies = np.array(
        [
            len(elem_edges[ei]) > 0 and all(edge_is_wet[e] for e in elem_edges[ei])
            for ei in range(nE)
        ],
        dtype=bool,
    )

    # Seeds must qualify
    seeds = _seed_elems_from_xy(mesh, seeds_xy)
    for ei in seeds:
        if not qualifies[ei]:
            raise RuntimeError(
                f"Seed element {ei} does not qualify (edge depth test failed). This can sometimes be caused by insufficient DEM coverage or by a .dem_cache directory being requeried after being populated with bad values )delete directory)"
            )

    # Stage 1 BFS
    in1 = np.zeros(nE, dtype=bool)

    q = deque()
    for ei in seeds:
        in1[ei] = True
        q.append(ei)
    while q:
        ei = q.popleft()
        for eidx in elem_edges[ei]:
            if not edge_is_wet[eidx]:
                continue
            nb = _neighbor_across(mesh, ei, eidx)
            if nb < 0 or in1[nb] or not qualifies[nb]:
                continue
            in1[nb] = True
            q.append(nb)

    # Stage 1 labels (nodewise, BOTH endpoints)
    perim1 = _perimeter_edges_of_region(mesh, in1)
    labels1 = {}
    for eidx in perim1:
        n1 = int(mesh.edges[eidx, 0])
        n2 = int(mesh.edges[eidx, 1])
        deep = mesh.nodes[n1, 2] >= (h0_node[n1] + delta0) and mesh.nodes[n2, 2] >= (
            h0_node[n2] + delta0
        )
        labels1[eidx] = "deep" if deep else "shore"

    # Stage 2: one-ring soft node test (node-centered h0)
    in2 = in1.copy()
    for eidx in perim1:
        ea, eb = int(mesh.edges[eidx, 3]), int(mesh.edges[eidx, 4])
        outside = eb if (ea >= 0 and in2[ea]) else (ea if (eb >= 0 and in2[eb]) else -1)
        if outside < 0 or not qualifies[outside]:
            continue
        n1, n2 = int(mesh.edges[eidx, 0]), int(mesh.edges[eidx, 1])
        if (node_depth[n1] >= (h0_node[n1] + delta1)) or (
            node_depth[n2] >= (h0_node[n2] + delta1)
        ):
            in2[outside] = True

    # Stage 2 labels (nodewise, BOTH endpoints)
    perim2 = _perimeter_edges_of_region(mesh, in2)
    labels2 = {}
    for eidx in perim2:
        n1 = int(mesh.edges[eidx, 0])
        n2 = int(mesh.edges[eidx, 1])
        deep = mesh.nodes[n1, 2] >= (h0_node[n1] + delta0) and mesh.nodes[n2, 2] >= (
            h0_node[n2] + delta0
        )
        labels2[eidx] = "deep" if deep else "shore"

    return FloodFillResult(in1, labels1, in2, labels2)


# --------------------------- Stage 3 (single-pass) ---------------------------


def _stage3_one_ring_candidates(
    mesh, in_region: np.ndarray, elem_edges, edge_mean_d, h0_edge, relax_factor, delta1
) -> Dict[int, int]:
    """
    Perimeter edge -> outside element, relaxed all-edges test using per-edge h0.
    Condition: edge_mean_d[ed] > (h0_edge[ed] - relax_factor*delta1) for all edges of the outside element.
    """
    candidates: Dict[int, int] = {}
    for eidx in _perimeter_edges_of_region(mesh, in_region):
        ea, eb = int(mesh.edges[eidx, 3]), int(mesh.edges[eidx, 4])
        if ea >= 0 and in_region[ea]:
            outside = eb
        elif eb >= 0 and in_region[eb]:
            outside = ea
        else:
            outside = -1
        if outside < 0:
            continue
        es = elem_edges[outside]
        if es and all(
            edge_mean_d[ed] > (h0_edge[ed] - relax_factor * delta1) for ed in es
        ):
            candidates[eidx] = outside
    return candidates


def _perim_vectors_at_node_with_added(
    mesh,
    elem_idx: int,
    node: int,
    base_mask: np.ndarray,
    add_elems: Set[int],
    elem_edges: List[List[int]],
) -> List[np.ndarray]:
    """True new-perimeter directions at node for elem_idx when ONLY add_elems are included."""
    xs, ys = mesh.nodes[:, 0], mesh.nodes[:, 1]
    vecs: List[np.ndarray] = []
    for ed in elem_edges[elem_idx]:
        a, b = int(mesh.edges[ed, 0]), int(mesh.edges[ed, 1])
        if node != a and node != b:
            continue
        ea, eb = int(mesh.edges[ed, 3]), int(mesh.edges[ed, 4])
        nb = eb if ea == elem_idx else (ea if eb == elem_idx else -2)
        if nb == -2:
            continue
        if nb < 0 or (not base_mask[nb] and nb not in add_elems):
            onode = b if node == a else a
            vecs.append(np.array([xs[onode] - xs[node], ys[onode] - ys[node]], float))
    return vecs


def _angle_improves_at_node(coords, prev_idx, curr_idx, next_idx, new_vecs, rad_eps):
    """
    Compare |turn angle| before vs. after using new-perimeter directions.
    IMPORTANT: if the corner disappears from the perimeter (no new_vecs),
    treat that as an improvement (after-angle ~ 0).
    """
    p_prev, p_curr, p_next = coords[prev_idx], coords[curr_idx], coords[next_idx]
    before = abs(_turn_angle(p_prev, p_curr, p_next))
    if not new_vecs:
        return before >= rad_eps
    if before == 0.0 and rad_eps > 0:
        return False
    for v in new_vecs:
        if np.allclose(v, 0.0):
            continue
        after = abs(_turn_angle(p_prev, p_curr, p_curr + v))
        if after <= max(0.0, before - rad_eps):
            return True
    return False


def stage3_single_pass_entry_exit(
    mesh,
    in_region_before: np.ndarray,
    elem_edges: List[List[int]],
    edge_mean_d: np.ndarray,
    h0_node: np.ndarray,
    h0_edge: np.ndarray,
    delta0: float,
    delta1: float,
    relax_factor: float,
    strip_nmax: int,
    angle_eps_deg: float,
    debug_nodes: Optional[List[int]] = None,
):
    xs, ys = mesh.nodes[:, 0], mesh.nodes[:, 1]
    coords = np.column_stack([xs, ys])

    cand_map = _stage3_one_ring_candidates(
        mesh, in_region_before, elem_edges, edge_mean_d, h0_edge, relax_factor, delta1
    )
    perim_before = _perimeter_edges_of_region(mesh, in_region_before)
    loops = _build_perimeter_loops(mesh, perim_before)

    in_after = in_region_before.copy()
    accepted_strips = 0
    elems_added = 0
    consumed_edges: Set[int] = set()
    rad_eps = np.deg2rad(max(0.0, angle_eps_deg))

    # Optional focused perimeter-edge listing near chosen nodes
    if debug_nodes:
        dbgset = set(int(n) for n in debug_nodes)
        print("Stage3 DEBUG: nodes =", sorted(dbgset))
        node2edges: Dict[int, List[int]] = {}
        for eidx in sorted(perim_before):
            a, b = int(mesh.edges[eidx, 0]), int(mesh.edges[eidx, 1])
            if a in dbgset or b in dbgset:
                node2edges.setdefault(a, []).append(eidx)
                node2edges.setdefault(b, []).append(eidx)
        for n in sorted(dbgset):
            eids = sorted(node2edges.get(n, []))
            print(f"  node {n}: incident perimeter edges = {eids}")

    for loop in loops:
        nodes = loop["nodes"]
        edges = loop["edges"]
        m = len(edges)
        if m == 0:
            continue

        def node_id(k):
            return nodes[k % len(nodes)]

        def edge_id(k):
            return edges[k % m]

        i = 0
        while i < m:
            e_start = edge_id(i)
            if e_start in consumed_edges or e_start not in cand_map:
                i += 1
                continue

            elem_out = cand_map[e_start]
            A, B = node_id(i), node_id(i + 1)
            A_prev, B_next = node_id(i - 1), node_id(i + 2)

            # local start-only simulation
            start_add = {elem_out}
            vecs_A = _perim_vectors_at_node_with_added(
                mesh, elem_out, A, in_region_before, start_add, elem_edges
            )
            vecs_B = _perim_vectors_at_node_with_added(
                mesh, elem_out, B, in_region_before, start_add, elem_edges
            )
            improves_A = _angle_improves_at_node(coords, A_prev, A, B, vecs_A, rad_eps)
            improves_B = _angle_improves_at_node(coords, A, B, B_next, vecs_B, rad_eps)

            if debug_nodes and (A in set(debug_nodes) or B in set(debug_nodes)):
                print(
                    f"  DEBUG start edge {e_start} at nodes A={A}, B={B}: "
                    f"cand_elem={elem_out}, improves_A={improves_A}, improves_B={improves_B}"
                )

            if not improves_A and not improves_B:
                i += 1
                continue

            direction = +1 if improves_A else -1

            strip_edges: List[int] = []
            strip_elems: List[int] = []
            j = i
            steps = 0
            terminal_found = False

            while steps < strip_nmax:
                ej = edge_id(j)
                if ej in consumed_edges or ej not in cand_map:
                    break
                elemj = cand_map[ej]
                strip_edges.append(ej)
                strip_elems.append(elemj)

                added = set(strip_elems)
                if direction == +1:
                    A_j, B_j, Bn_j = node_id(j), node_id(j + 1), node_id(j + 2)
                    vecs_far = _perim_vectors_at_node_with_added(
                        mesh, elemj, B_j, in_region_before, added, elem_edges
                    )
                    improves_far = _angle_improves_at_node(
                        coords, A_j, B_j, Bn_j, vecs_far, rad_eps
                    )
                else:
                    Ap_j, A_j, B_j = node_id(j - 1), node_id(j), node_id(j + 1)
                    vecs_far = _perim_vectors_at_node_with_added(
                        mesh, elemj, A_j, in_region_before, added, elem_edges
                    )
                    improves_far = _angle_improves_at_node(
                        coords, Ap_j, A_j, B_j, vecs_far, rad_eps
                    )

                if improves_far:
                    terminal_found = True
                    break

                j += direction
                steps += 1

            if terminal_found and strip_elems:
                # Focused debug: list accepted edges/nodes if touching debug nodes
                if debug_nodes:
                    dbgset = set(int(n) for n in debug_nodes)
                    touches_dbg = any(
                        int(mesh.edges[e, 0]) in dbgset
                        or int(mesh.edges[e, 1]) in dbgset
                        for e in strip_edges
                    )
                    if touches_dbg:
                        print(
                            f"ACCEPT (near debug nodes): start_edge={e_start}, "
                            f"dir={'forward' if direction==+1 else 'backward'}, "
                            f"len={len(strip_edges)}"
                        )
                        uniq_nodes = []
                        for e in strip_edges:
                            n1, n2 = int(mesh.edges[e, 0]), int(mesh.edges[e, 1])
                            print(f"  edge {e}: nodes ({n1}, {n2})")
                            uniq_nodes.extend([n1, n2])
                        print(f"  unique_nodes_in_strip: {sorted(set(uniq_nodes))}")

                for ee in strip_elems:
                    if not in_after[ee]:
                        in_after[ee] = True
                        elems_added += 1
                for eused in strip_edges:
                    consumed_edges.add(eused)
                accepted_strips += 1
                i = (j + 1) if direction == +1 else (i + 1)
            else:
                i += 1

    # Final perimeter + labels (node-centered h0 applied per-edge)
    perim_final = _perimeter_edges_of_region(mesh, in_after)
    node_depth = mesh.nodes[:, 2].astype(float)
    edge_mean_d_final = _edge_mean_depth(mesh, node_depth)
    labels_final = {}
    for eidx in perim_final:
        n1 = int(mesh.edges[eidx, 0])
        n2 = int(mesh.edges[eidx, 1])
        deep_nodewise = (node_depth[n1] >= (h0_node[n1] + delta0)) and (
            node_depth[n2] >= (h0_node[n2] + delta0)
        )
        labels_final[eidx] = "deep" if deep_nodewise else "shore"

    return in_after, labels_final, accepted_strips, elems_added


# ------------------------------ I/O helpers ----------------------------------


def _write_perimeter_csv(mesh, labels: Dict[int, str], path: str) -> None:

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["edge_id", "n1", "n2", "label"])
        for eidx, lab in sorted(labels.items()):
            n1, n2 = int(mesh.edges[eidx, 0]), int(mesh.edges[eidx, 1])
            w.writerow([eidx, n1, n2, lab])


def _write_wet_elems_csv(mask: np.ndarray, path: str) -> None:

    idxs = np.flatnonzero(mask.astype(bool))
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["elem_id"])
        for i in idxs:
            w.writerow([int(i)])


def _write_perimeter_shapefile(
    mesh, labels: Dict[int, str], path: str, epsg: int = None
) -> None:
    """
    Write 2D polyline Shapefile (.shp + .shx + .dbf, and .prj if EPSG given) for perimeter edges.
    Uses geopandas + shapely only.
    """

    xs, ys = mesh.nodes[:, 0], mesh.nodes[:, 1]
    node_depth = mesh.nodes[:, 2].astype(float)
    edge_mean_d = _edge_mean_depth(mesh, node_depth)
    if not path.lower().endswith(".shp"):
        path += ".shp"

    polylines = []
    for eidx, lab in labels.items():
        n1, n2 = int(mesh.edges[eidx, 0]), int(mesh.edges[eidx, 1])
        geom = LineString(
            [[(float(xs[n1]), float(ys[n1])), (float(xs[n2]), float(ys[n2]))]]
        )
        polylines.append(
            {
                "edge_id": int(eidx),
                "n1": n1,
                "n2": n2,
                "label": str(label),
                "geometry": geom,
                "depth_mean": float(edge_mean_d[eidx]),
            }
        )

    gdf = gpd.GeoDataFrame(polylines)

    if epsg is not None:
        gdf.set_crs(epsg=epsg, inplace=True)

    gdf.to_file(path, driver="ESRI Shapefile")


# ------------------------ Merging + Visualization -----------------------------


def _edge_pairs(mesh, edge_ids: Iterable[int]) -> List[Tuple[int, int]]:
    pairs = []
    for e in edge_ids:
        n1, n2 = int(mesh.edges[int(e), 0]), int(mesh.edges[int(e), 1])
        pairs.append((n1, n2))
    return pairs


def _build_chains_from_edges(pairs: List[Tuple[int, int]]) -> List[List[int]]:
    """
    Given undirected edge (n1,n2) pairs, split into maximal paths or cycles and
    return ordered node lists (one per chain).
    """

    adj = defaultdict(list)
    edge_set = set()
    for a, b in pairs:
        adj[a].append(b)
        adj[b].append(a)
        edge_set.add((min(a, b), max(a, b)))
    used = set()
    chains = []

    def consume(start):
        path = [start]
        cur = start
        prev = None
        while True:
            nexts = [
                n
                for n in adj[cur]
                if (min(cur, n), max(cur, n)) in edge_set
                and (min(cur, n), max(cur, n)) not in used
            ]
            if prev is not None and prev in nexts and len(nexts) > 1:
                nexts.remove(prev)
            if not nexts:
                break
            nxt = nexts[0]
            used.add((min(cur, nxt), max(cur, nxt)))
            path.append(nxt)
            prev, cur = cur, nxt
        return path

    # endpoints first (deg != 2)
    endpoints = [n for n, neis in adj.items() if len(neis) != 2]
    for ep in endpoints:
        if any((min(ep, k), max(ep, k)) not in used for k in adj[ep]):
            chains.append(consume(ep))
    # cycles
    for n in list(adj.keys()):
        if any((min(n, k), max(n, k)) not in used for k in adj[n]):
            chains.append(consume(n))
    return chains


def _merged_polylines_and_node_rows(mesh, labels_all: Dict[int, str]):
    xs, ys = mesh.nodes[:, 0], mesh.nodes[:, 1]
    zs = mesh.nodes[:, 2].astype(float)
    by_label: Dict[str, List[int]] = {}
    for eidx, lab in labels_all.items():
        by_label.setdefault(str(lab), []).append(int(eidx))

    # NEW: nodes that touch any DEEP edge (used to trim shore/smoothed endpoints)
    deep_nodes = set()
    if "deep" in by_label:
        for e in by_label["deep"]:
            n1, n2 = int(mesh.edges[e, 0]), int(mesh.edges[e, 1])
            deep_nodes.add(n1)
            deep_nodes.add(n2)

    polylines = []
    node_rows = []
    poly_id = 0
    for lab, eids in by_label.items():
        pairs = _edge_pairs(mesh, eids)
        chains = _build_chains_from_edges(pairs)
        for chain in chains:
            # NEW: trim shore/smoothed endpoints that sit on deep nodes
            if lab in ("shore", "smoothed"):
                # trim from the front
                while len(chain) > 2 and chain[0] in deep_nodes:
                    chain = chain[1:]
                # trim from the back
                while len(chain) > 2 and chain[-1] in deep_nodes:
                    chain = chain[:-1]
                if len(chain) < 2:
                    continue  # nothing left to draw / export

            coords = np.column_stack([xs[chain], ys[chain]]).tolist()
            polylines.append(
                dict(poly_id=poly_id, label=lab, nodes=chain, coords=coords)
            )
            for n in chain:
                node_rows.append(
                    dict(
                        node=int(n) + 1,
                        x=float(xs[n]),
                        y=float(ys[n]),
                        z=float(zs[n]),
                        poly_id=int(poly_id),
                        label=lab,
                    )
                )
            poly_id += 1
    return polylines, node_rows


def _write_merged_shapefile(polylines, path: str, epsg: int = None):
    """
    Write 2D polyline Shapefile (.shp + .shx + .dbf, and .prj if EPSG given) for perimeter edges.
    Uses geopandas + shapely only.
    """

    if not path.lower().endswith(".shp"):
        path += ".shp"

    # Build a GeoDataFrame from your polylines
    gdf = gpd.GeoDataFrame(
        [
            {
                "poly_id": int(p["poly_id"]),
                "label": str(p["label"]),
                "geometry": LineString(p["coords"]),
            }
            for p in polylines
        ]
    )

    # Set CRS if provided
    if epsg is not None:
        gdf.set_crs(epsg=epsg, inplace=True)

    # Write to shapefile
    gdf.to_file(path, driver="ESRI Shapefile")


def _write_nodes_csv(node_rows, path: str, epsg: int = None):

    with open(path, "w", newline="") as f:
        if epsg:
            f.write(
                f"# Node ids are 1-based, xy coordinates are based on EPSG:{int(epsg)}\n"
            )
        else:
            f.write("# Node ids are 1-based\n")
        w = csv.writer(f)
        w.writerow(["node", "x", "y", "z", "poly_id", "label"])
        for r in node_rows:
            w.writerow(
                [
                    r["node"],
                    f"{r['x']:.6f}",
                    f"{r['y']:.6f}",
                    f"{r['z']:.6f}",
                    r["poly_id"],
                    r["label"],
                ]
            )


def _plot_profiles(
    polylines,
    node_rows,
    which_ids: Set[int],
    out_dir: str,
    filt: str = "none",
    window: int = 5,
):

    rows_by_node = {int(r["node"]): r for r in node_rows}
    for p in polylines:
        pid = int(p["poly_id"])
        if which_ids and (pid not in which_ids):
            continue
        seq = []
        for n0 in p["nodes"]:
            n1 = int(n0) + 1
            r = rows_by_node.get(n1)
            if r is None:  # may be filtered out of CSV
                continue
            seq.append((r["x"], r["y"], r["z"], r["label"]))
        if len(seq) < 2:
            continue
        xs = np.array([s[0] for s in seq])
        ys = np.array([s[1] for s in seq])
        zs = np.array([s[2] for s in seq])
        labels = [s[3] for s in seq]
        dxy = np.sqrt(np.diff(xs) ** 2 + np.diff(ys) ** 2)
        s = np.concatenate([[0.0], np.cumsum(dxy)])

        def _roll(arr, k, mode):
            k = max(1, int(k) // 2 * 2 + 1)  # ensure odd
            if mode == "none" or k <= 1:
                return None, k
            pad = k // 2
            apad = np.pad(arr, (pad, pad), mode="edge")
            out = np.empty_like(arr)
            for i in range(len(arr)):
                win = apad[i : i + k]
                if mode == "median":
                    out[i] = np.median(win)
                elif mode == "max":
                    out[i] = np.max(win)
                else:
                    return None, k
            return out, k

        zf, k = _roll(zs, window, filt.lower())
        fig, ax = plt.subplots(figsize=(8, 3.0))
        (base_line,) = ax.plot(s, zs, label=f"z (depth) — {p['label']}")

        if zf is not None:
            ax.plot(
                s,
                zf,
                linestyle="--",
                label=f"{filt.lower()} (k={k})",
                color=base_line.get_color(),
            )
        ax.set_xlabel("distance (m)")
        ax.set_ylabel("depth (m, +down)")
        ax.set_xlabel("distance (m)")
        ax.invert_yaxis()  # draw positive depth downward (negative up)
        ax.set_title(f"poly_id={pid}, label={p['label']}")
        ax.legend()
        fig.tight_layout()
        out_png = os.path.join(out_dir or ".", f"profile_poly_{pid}.png")
        fig.savefig(out_png, dpi=140)
        plt.close(fig)


# ------------------- Variational/TVD refinement -------------------

# ---- Lightweight geometry/quadrature and TV machinery (adapted & trimmed) ----


class AdaptiveElementQuadrature:
    """Adaptive quadrature for tri/quad elements.

    Triangles use a fixed symmetric 6-point rule. Quads use tensor
    Gauss–Legendre nodes (3×3, 4×4, or 5×5) chosen by polygon area.

    Parameters
    ----------
    A9 : float, optional
        Area threshold (m^2) below which quads use 3×3 rule, by default 60.0.
    A16 : float, optional
        Area threshold for 4×4 rule (5×5 above), by default 200.0.
    """

    def __init__(self, A9: float = 60.0, A16: float = 200.0):
        self.A9 = float(A9)
        self.A16 = float(A16)

    @staticmethod
    def _poly_area(xy: np.ndarray) -> float:
        a = 0.0
        n = xy.shape[0]
        for i in range(n):
            x1, y1 = xy[i]
            x2, y2 = xy[(i + 1) % n]
            a += x1 * y2 - x2 * y1
        return 0.5 * abs(a)

    @staticmethod
    def _tri_rule6() -> Tuple[np.ndarray, np.ndarray]:
        a = 0.445948490915965
        b = 0.091576213509771
        w1 = 0.111690794839005
        w2 = 0.054975871827661
        pts = np.array(
            [
                [a, a, 1 - 2 * a],
                [a, 1 - 2 * a, a],
                [1 - 2 * a, a, a],
                [b, b, 1 - 2 * b],
                [b, 1 - 2 * b, b],
                [1 - 2 * b, b, b],
            ]
        )
        wts = np.array([w1, w1, w1, w2, w2, w2])
        return pts, wts

    @staticmethod
    def _map_tri(xy: np.ndarray, bary: np.ndarray) -> np.ndarray:
        return bary @ xy

    @staticmethod
    def _gauss_legendre_1d(n: int) -> Tuple[np.ndarray, np.ndarray]:

        return leggauss(n)

    @staticmethod
    def _map_quad(xy: np.ndarray, xi_eta: np.ndarray):
        x0, y0 = xy[0, 0], xy[0, 1]
        x1, y1 = xy[1, 0], xy[1, 1]
        x2, y2 = xy[2, 0], xy[2, 1]
        x3, y3 = xy[3, 0], xy[3, 1]

        xi = xi_eta[:, 0]
        eta = xi_eta[:, 1]

        N0 = 0.25 * (1 - xi) * (1 - eta)
        N1 = 0.25 * (1 + xi) * (1 - eta)
        N2 = 0.25 * (1 + xi) * (1 + eta)
        N3 = 0.25 * (1 - xi) * (1 + eta)

        x = N0 * x0 + N1 * x1 + N2 * x2 + N3 * x3
        y = N0 * y0 + N1 * y1 + N2 * y2 + N3 * y3
        pts = np.column_stack([x, y])

        dN0_dxi = -0.25 * (1 - eta)
        dN1_dxi = 0.25 * (1 - eta)
        dN2_dxi = 0.25 * (1 + eta)
        dN3_dxi = -0.25 * (1 + eta)

        dN0_deta = -0.25 * (1 - xi)
        dN1_deta = -0.25 * (1 + xi)
        dN2_deta = 0.25 * (1 + xi)
        dN3_deta = 0.25 * (1 - xi)

        dx_dxi = dN0_dxi * x0 + dN1_dxi * x1 + dN2_dxi * x2 + dN3_dxi * x3
        dy_dxi = dN0_dxi * y0 + dN1_dxi * y1 + dN2_dxi * y2 + dN3_dxi * y3
        dx_deta = dN0_deta * x0 + dN1_deta * x1 + dN2_deta * x2 + dN3_deta * x3
        dy_deta = dN0_deta * y0 + dN1_deta * y1 + dN2_deta * y2 + dN3_deta * y3

        J = np.abs(dx_dxi * dy_deta - dy_dxi * dx_deta)
        return pts, J

    def quad_points(self, xy: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        area = self._poly_area(xy)
        if xy.shape[0] == 3:
            pts_b, w_b = self._tri_rule6()
            pts = self._map_tri(xy, pts_b)
            wts = 2.0 * w_b * area
            return pts, wts
        elif xy.shape[0] == 4:
            if area < self.A9:
                n = 3
            elif area < self.A16:
                n = 4
            else:
                n = 5
            xi, wi = self._gauss_legendre_1d(n)
            XI, ETA = np.meshgrid(xi, xi)
            WI, WJ = np.meshgrid(wi, wi)
            pts, J = self._map_quad(xy, np.c_[XI.ravel(), ETA.ravel()])
            wts = (WI * WJ).ravel() * J
            # normalize weights to polygon area (robust to warping)
            A_poly = AdaptiveElementQuadrature._poly_area(xy)
            Wsum = float(np.sum(wts))
            if not np.isfinite(Wsum) or Wsum <= 0.0:
                wts = np.full_like(wts, A_poly / len(wts))
            else:
                wts *= A_poly / Wsum
            return pts, wts
        else:
            raise ValueError("Only tri/quad supported")


class TVSubgradientOperator:
    """Edge-TV subgradient step on an unstructured mesh.

    TV(z) = sum_{(i,j)} w_ij * |z_i - z_j|, with w_ij = |edge|.

    Parameters
    ----------
    nodes_xy : ndarray, shape (n,2)
        Node x,y coordinates.
    elems : list of array_like
        Element connectivities (3 or 4 node indices).
    lump_mass : ndarray, shape (n,)
        Lumped node areas (for explicit scaling).
    eps : float, optional
        Smoothing parameter in the |·| derivative, by default 1e-8.
    enforce_cfl : bool, optional
        If True, scales the step to keep dt*max_rate <= cfl_target.
    cfl_target : float, optional
        Target explicit rate cap used when ``enforce_cfl`` is True.
    clipping_eps : Optional[float], optional
        If not None, clip the updated node values to each node's neighbor
        range ± ``clipping_eps``. Use None to disable clipping entirely.
    """

    def __init__(
        self,
        nodes_xy: np.ndarray,
        elems: List[np.ndarray],
        lump_mass: np.ndarray,
        eps: float = 1e-8,
        enforce_cfl: bool = True,
        cfl_target: float = 1.0,
        clipping_eps: Optional[float] = None,
    ):
        self.nodes_xy = np.asarray(nodes_xy)
        self.elems = elems
        self.M = np.asarray(lump_mass, dtype=float)
        self.eps = float(eps)
        self.enforce_cfl = bool(enforce_cfl)
        self.cfl_target = float(cfl_target)
        self.clipping_eps = clipping_eps

        # Undirected edge list (i<j)
        edge_set = set()
        for conn in elems:
            m = len(conn)
            for a in range(m):
                i = int(conn[a])
                j = int(conn[(a + 1) % m])
                if i > j:
                    i, j = j, i
                edge_set.add((i, j))
        edges = np.array(sorted(edge_set), dtype=int)
        self.edges = edges

        p0 = self.nodes_xy[edges[:, 0], :]
        p1 = self.nodes_xy[edges[:, 1], :]
        self.w = np.linalg.norm(p1 - p0, axis=1)  # edge length

        # node incident edges
        n = self.nodes_xy.shape[0]
        inc = [[] for _ in range(n)]
        for e, (i, j) in enumerate(edges):
            inc[i].append(e)
            inc[j].append(e)
        self.incident = [np.asarray(ix, dtype=int) for ix in inc]

        # precompute rate_i ≈ (sum_j w_ij) / M_i
        wsum = np.zeros(n, dtype=float)
        for e, (i, j) in enumerate(edges):
            wij = self.w[e]
            wsum[i] += wij
            wsum[j] += wij
        self.rate = wsum / np.maximum(self.M, 1e-12)
        self.max_rate = float(np.max(self.rate)) if self.rate.size else 0.0

    def step(
        self, z: np.ndarray, mu: float = 1.0, dt: float = 1.0
    ) -> Tuple[np.ndarray, Dict]:
        z = np.asarray(z, dtype=float)

        # CFL-cap alpha so that alpha*max_rate <= cfl_target
        scale = 1.0
        if self.enforce_cfl and self.max_rate > 0.0 and dt > 0.0:
            scale = min(1.0, self.cfl_target / (self.max_rate * dt))
        alpha = dt * scale * float(mu)

        # subgradient
        g = np.zeros_like(z)
        i_idx = self.edges[:, 0]
        j_idx = self.edges[:, 1]
        dz = z[i_idx] - z[j_idx]
        sgn = dz / np.sqrt(dz * dz + self.eps * self.eps)
        wsgn = self.w * sgn
        np.add.at(g, i_idx, wsgn)
        np.add.at(g, j_idx, -wsgn)

        z_trial = z - alpha * (g / np.maximum(self.M, 1e-12))

        if self.clipping_eps is not None:
            # neighbor range per node
            n = z.shape[0]
            zmin = np.full(n, +np.inf)
            zmax = np.full(n, -np.inf)
            for node, inc in enumerate(self.incident):
                if inc.size == 0:
                    zmin[node] = z[node]
                    zmax[node] = z[node]
                    continue
                # gather neighbors from edges incident to this node
                nbrs = set()
                for e in inc:
                    i, j = self.edges[e]
                    nbrs.add(i)
                    nbrs.add(j)
                nbrs.discard(node)
                if nbrs:
                    zi = z[list(nbrs)]
                    zmin[node] = zi.min()
                    zmax[node] = zi.max()
                else:
                    zmin[node] = z[node]
                    zmax[node] = z[node]
            epsc = float(self.clipping_eps)
            z_trial = np.minimum(np.maximum(z_trial, zmin - epsc), zmax + epsc)

        info = {"scale": float(scale), "max_rate": float(self.max_rate)}
        return z_trial, info


class TVL2VolumeEvaluator:
    """Evaluator for TV energy and DEM-integrated L2 energy with cached quadrature.

    Parameters
    ----------
    mesh : SCHISM_mesh
        Mesh with ``nodes`` (n,3), ``edges`` (m,2), and ``elems`` (list of 3/4-tuples).
    dem_sampler : callable, optional
        Function mapping ``(k,2)`` XY points → DEM depth (+down) array.
    quadrature : AdaptiveElementQuadrature, optional
        Quadrature rule provider.
    """

    def __init__(
        self,
        mesh: SCHISM_mesh,
        dem_sampler: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        quadrature: Optional[AdaptiveElementQuadrature] = None,
    ):
        self.mesh = mesh
        self.dem = dem_sampler
        self.quad = quadrature or AdaptiveElementQuadrature()

        self.elem_quad_pts: List[np.ndarray] = []
        self.elem_quad_wts: List[np.ndarray] = []
        self.elem_kind: List[Tuple[str, int]] = []  # [("tri",6) or ("quad", n)]
        self.elem_areas = self._elem_areas()

        E = self.mesh.edges[:, :2].astype(int, copy=False)
        pi = self.mesh.nodes[E[:, 0], :2]
        pj = self.mesh.nodes[E[:, 1], :2]
        dij = np.linalg.norm(pi - pj, axis=1)
        self._tv_edges = E
        self._tv_w = 1.0 / (dij + 1e-12)

        self.dem_vals: Optional[List[np.ndarray]] = None
        self._S_tri6: Optional[np.ndarray] = None
        self._S_quad: Dict[int, np.ndarray] = {}
        self.dem_abs_max: float = 140.0  # meters

        if self.dem is not None:
            self._prep_dem()

    def _elem_areas(self) -> np.ndarray:
        A = np.empty(self.mesh.n_elems(), dtype=float)
        for ei, conn in enumerate(self.mesh.elems):
            xy = self.mesh.nodes[np.asarray(conn), :2]
            A[ei] = AdaptiveElementQuadrature._poly_area(xy)
        return A

    @staticmethod
    def _build_S_tri6() -> np.ndarray:
        a = 0.445948490915965
        b = 0.091576213509771
        return np.array(
            [
                [a, a, 1 - 2 * a],
                [a, 1 - 2 * a, a],
                [1 - 2 * a, a, a],
                [b, b, 1 - 2 * b],
                [b, 1 - 2 * b, b],
                [1 - 2 * b, b, b],
            ],
            dtype=float,
        )

    @staticmethod
    def _build_S_quad(n: int) -> np.ndarray:

        xi1d, _ = leggauss(n)
        XI, ETA = np.meshgrid(xi1d, xi1d)
        xi, eta = XI.ravel(), ETA.ravel()
        N0 = 0.25 * (1 - xi) * (1 - eta)
        N1 = 0.25 * (1 + xi) * (1 - eta)
        N2 = 0.25 * (1 + xi) * (1 + eta)
        N3 = 0.25 * (1 - xi) * (1 + eta)
        return np.stack([N0, N1, N2, N3], axis=1)

    def _prep_dem(self):
        pts_to_query = []
        self.elem_quad_pts.clear()
        self.elem_quad_wts.clear()
        self.elem_kind.clear()

        for conn in self.mesh.elems:
            ids = np.asarray(conn)
            xy = self.mesh.nodes[ids, :2]
            pts, wts = self.quad.quad_points(xy)
            self.elem_quad_pts.append(pts)
            self.elem_quad_wts.append(wts)
            if ids.size == 3:
                self.elem_kind.append(("tri", 6))
            elif ids.size == 4:
                n = int(round(np.sqrt(len(wts))))
                self.elem_kind.append(("quad", n))
            else:
                raise ValueError("Only tri/quad elements supported")
            pts_to_query.append(pts)

        all_pts = np.vstack(pts_to_query)
        all_dem = self.dem(all_pts).astype(float)

        bad = ~np.isfinite(all_dem) | (np.abs(all_dem) > self.dem_abs_max)
        if int(bad.sum()) > 0:
            idx = np.flatnonzero(bad)[:10]
            raise ValueError(
                "DEM contains invalid/out-of-range samples at prep "
                f"(examples indices: {idx.tolist()})."
            )

        self.dem_vals = []
        row = 0
        for wts in self.elem_quad_wts:
            n = len(wts)
            self.dem_vals.append(all_dem[row : row + n])
            row += n

    def tv_energy(self, z: np.ndarray) -> float:
        E = self._tv_edges
        return float(np.sum(self._tv_w * np.abs(z[E[:, 0]] - z[E[:, 1]])))

    def l2_dem(self, z: np.ndarray) -> float:
        if self.dem is None or self.dem_vals is None:
            return 0.0
        if self._S_tri6 is None:
            self._S_tri6 = self._build_S_tri6()
        needed_ns = {n for kind, n in self.elem_kind if kind == "quad"}
        for n in needed_ns:
            if n not in self._S_quad:
                self._S_quad[n] = self._build_S_quad(n)

        total = 0.0
        total_area = float(np.sum(self.elem_areas)) + 1e-16
        for ei, conn in enumerate(self.mesh.elems):
            ids = np.asarray(conn, dtype=int)
            wts = self.elem_quad_wts[ei]
            dem = np.asarray(self.dem_vals[ei], dtype=float)
            kind, tag = self.elem_kind[ei]
            if kind == "tri":
                zq = self._S_tri6 @ z[ids]
            else:
                S = self._S_quad[tag]
                zq = S @ z[ids]
            bad = (
                (~np.isfinite(dem))
                | (~np.isfinite(zq))
                | (np.abs(dem) > self.dem_abs_max)
                | (np.abs(zq) > self.dem_abs_max)
            )
            if np.any(bad):
                raise ValueError("Invalid DEM/mesh values encountered during L2.")
            total += float(np.sum(wts * (zq - dem) ** 2))
        return float(total / total_area)

    def global_volume(self, z: np.ndarray) -> float:
        V = 0.0
        for ei, conn in enumerate(self.mesh.elems):
            zbar = float(np.mean(z[np.asarray(conn, dtype=int)]))
            V += self.elem_areas[ei] * zbar
        return float(V)


class VolumeL2Variational:
    """Element-mean variational L2 penalty to DEM.

    J_vol(z) = 0.5 * sum_e A_e * ( z̄_e - d̄_e )^2,   r_e = z̄_e - d̄_e

    Gradient w.r.t. node i:
    dJ/dz_i = sum_{e∋i} r_e * ∫_e N_i dA
    """

    def __init__(self, evaluator: TVL2VolumeEvaluator):
        self.ev = evaluator
        if self.ev.dem is None or self.ev.dem_vals is None:
            raise ValueError("Variational L2 requires a DEM sampler.")
        if self.ev._S_tri6 is None:
            self.ev._S_tri6 = self.ev._build_S_tri6()
        needed_ns = {n for kind, n in self.ev.elem_kind if kind == "quad"}
        for n in needed_ns:
            if n not in self.ev._S_quad:
                self.ev._S_quad[n] = self.ev._build_S_quad(n)

        self.Ni_int: List[np.ndarray] = []
        for ei, conn in enumerate(self.ev.mesh.elems):
            ids = np.asarray(conn, dtype=int)
            wts = self.ev.elem_quad_wts[ei]
            kind, tag = self.ev.elem_kind[ei]
            if kind == "tri":
                S = self.ev._S_tri6  # (6,3)
                Ni = (wts[:, None] * S).sum(axis=0)  # (3,)
            else:
                n = tag
                S = self.ev._S_quad[n]  # (n*n,4)
                Ni = (wts[:, None] * S).sum(axis=0)  # (4,)
            self.Ni_int.append(Ni.astype(float))

    def residuals_and_energy(self, z: np.ndarray) -> Tuple[np.ndarray, float]:
        r = np.zeros(self.ev.mesh.n_elems(), dtype=float)
        total = 0.0
        for ei, conn in enumerate(self.ev.mesh.elems):
            ids = np.asarray(conn, dtype=int)
            wts = self.ev.elem_quad_wts[ei]
            dem = np.asarray(self.ev.dem_vals[ei], dtype=float)
            kind, tag = self.ev.elem_kind[ei]
            A = float(np.sum(wts))
            if kind == "tri":
                zq = self.ev._S_tri6 @ z[ids]
            else:
                n = tag
                zq = self.ev._S_quad[n] @ z[ids]
            zbar = float(np.sum(wts * zq) / (A + 1e-16))
            dbar = float(np.sum(wts * dem) / (A + 1e-16))
            ri = zbar - dbar
            r[ei] = ri
            total += 0.5 * A * (ri * ri)
        return r, total

    def gradient(self, r: np.ndarray) -> np.ndarray:
        g = np.zeros(self.ev.mesh.n_nodes(), dtype=float)
        for ei, conn in enumerate(self.ev.mesh.elems):
            ids = np.asarray(conn, dtype=int)
            Ni = self.Ni_int[ei]
            g[ids] += r[ei] * Ni
        return g


# ----------------------------- Options & Driver -------------------------------


@dataclass
class FlowOptions:
    """Options for :func:`run_tvd_variational`.

    Attributes
    ----------
    n_steps : int
        Max number of iterations.
    dt : Optional[float]
        Base time step (meters). If None, a CFL-derived dt is used.
    cfl : float
        CFL-like factor if ``dt`` is None.
    mu : float
        TV step strength (prox gain). Default 1.0.
    lambda_l2 : float
        Scale on the L2 descent increment (separate from ``l2_weight`` in the
        objective). Default 1.0.
    cfl_target : float
        Target explicit rate cap for TV step. Default 1.0.
    clipping_eps : Optional[float]
        If not None, clip neighbor-wise with ± ``clipping_eps`` (meters).
    rel_eta : float
        Required relative drop in the weighted objective per accept.
    max_backtracks : int
        Maximum backtracks per iteration.
    tv_weight : float
        Weight on TV in the objective.
    l2_weight : float
        Weight on L2(volume) in the objective.
    """

    n_steps: int = 25
    dt: Optional[float] = None
    cfl: float = 0.25
    mu: float = 1.0
    lambda_l2: float = 1.0
    cfl_target: float = 1.0
    clipping_eps: Optional[float] = None
    rel_eta: float = 1e-3
    max_backtracks: int = 4
    tv_weight: float = 1.0
    l2_weight: float = 0.2


def run_tvd_variational(
    mesh: SCHISM_mesh,
    z0: np.ndarray,
    dem_sampler: Callable[[np.ndarray], np.ndarray],
    options: FlowOptions = FlowOptions(),
    z_floor: np.ndarray = None,
    floor_mask: np.ndarray = None,
    reg_weight: float = 0.0,
    enforce_floor: bool = False,
    logger: logging.Logger = None,
) -> Dict[str, np.ndarray]:
    """TV-prox + variational L2 descent with backtracking.

    The iteration performs an explicit L2 descent step followed by a TV
    subgradient "prox" step. Acceptance is based on a relative drop in the
    weighted objective :math:`w_{tv} TV + w_{l2} J_{vol}`.

    Parameters
    ----------
    mesh : SCHISM_mesh
        Mesh object providing ``nodes``, ``edges``, and ``elems``.
    z0 : ndarray, shape (n_nodes,)
        Initial node elevations (meters). Use DEM samples if your incoming mesh
        has been manipulated away from the DEM and you want a safer start.
    dem_sampler : callable
        Function mapping XY→depth(+down). See :func:`create_dem_sampler`.
    options : FlowOptions, optional
        Algorithm options (see :class:`FlowOptions`).

    Returns
    -------
    dict
        Keys:
        - ``"z"``: final node elevations.
        - ``"history"``: (k+1, n_nodes) stack of accepted iterates.
    """
    evaluator = TVL2VolumeEvaluator(
        mesh, dem_sampler=dem_sampler, quadrature=AdaptiveElementQuadrature()
    )
    evaluator.dem_abs_max = 140.0

    # lumped "mass" per node for explicit scaling (area sum of incident elements / nen)
    a_node = np.zeros(mesh.n_nodes(), dtype=float)
    all_areas = evaluator.elem_areas
    for ei, conn in enumerate(mesh.elems):
        area = all_areas[ei]
        nen = len(conn)
        contrib = area / float(nen)
        for ni in conn:
            a_node[int(ni)] += contrib
    a_node = np.maximum(a_node, 1e-12)

    tvd = TVSubgradientOperator(
        mesh.nodes[:, :2],
        mesh.elems,
        lump_mass=a_node,
        eps=1e-8,
        enforce_cfl=True,
        cfl_target=float(options.cfl_target),
        clipping_eps=options.clipping_eps,
    )

    # dt setup
    hmin = np.nanmin(
        np.linalg.norm(
            mesh.nodes[mesh.edges[:, 0], :2] - mesh.nodes[mesh.edges[:, 1], :2], axis=1
        )
    )
    dt_base = options.dt if options.dt is not None else options.cfl * (hmin**2)

    # Initial metrics
    z = z0.copy()
    tv_last = evaluator.tv_energy(z)
    vol = VolumeL2Variational(evaluator)
    _, Jvol_last = vol.residuals_and_energy(z)

    w_tv = float(options.tv_weight)
    w_l2 = float(options.l2_weight)
    obj_last = w_tv * tv_last + w_l2 * Jvol_last

    # Optional shoreline-floor regularizer controls
    w_reg = float(reg_weight)
    has_floor = (z_floor is not None) and (floor_mask is not None) and (w_reg > 0.0)

    max_backtracks = int(options.max_backtracks)
    rel_eta = float(options.rel_eta)

    history = [z.copy()]

    for k in range(options.n_steps):
        accepted = False
        bt = 0
        dt_try = float(dt_base)
        logger.debug(
            f"[iter {k:03d}] dt={dt_try:.3e} tv_last={tv_last:.3e} Jvol_last={Jvol_last:.3e}"
        )

        while True:
            # --- L2 gradient step ---
            r_before, J_before = vol.residuals_and_energy(z)
            g_vol = vol.gradient(r_before)
            # One-sided hinge gradient only where z < z_floor
            if has_floor:
                act = (floor_mask) & (z < z_floor)
                g_reg = np.zeros_like(z)
                g_reg[act] = z[act] - z_floor[act]
            else:
                g_reg = 0.0
            dz_l2 = -dt_try * (w_l2 * g_vol + w_reg * g_reg) / a_node
            z_lin = z + float(options.lambda_l2) * dz_l2

            # --- TV prox step ---
            z_trial, info = tvd.step(z_lin, mu=float(options.mu), dt=dt_try)
            if enforce_floor and has_floor:
                z_trial[floor_mask] = np.maximum(
                    z_trial[floor_mask], z_floor[floor_mask]
                )
            cfl_scale = float(info.get("scale", 1.0))
            max_rate = float(info.get("max_rate", 0.0))
            logger.debug(
                f"[iter {k:03d}] TV step: cfl_scale={cfl_scale:.3f} max_rate={max_rate:.3e} mu={options.mu:.2f}"
            )

            # --- Diagnostics: TV progression
            tv_before = evaluator.tv_energy(z)
            tv_after_l2 = evaluator.tv_energy(z_lin)
            tv_after_prox = evaluator.tv_energy(z_trial)
            dTV_L2 = tv_after_l2 - tv_before
            dTV_prox = tv_after_prox - tv_after_l2
            dTV_total = tv_after_prox - tv_before
            logger.debug(
                f"[iter {k:03d}] TV: before={tv_before:.6e} -> afterL2={tv_after_l2:.6e} "
                f"(ΔL2={dTV_L2:.3e}) -> afterTV={tv_after_prox:.6e} (ΔTV={dTV_prox:.3e}, Δtot={dTV_total:.3e})"
            )

            # --- L2 energy progression
            _, J_after_l2 = vol.residuals_and_energy(z_lin)
            _, J_trial = vol.residuals_and_energy(z_trial)
            logger.debug(
                f"[iter {k:03d}] J_vol: before={J_before:.6e} -> afterL2={J_after_l2:.6e} -> afterTV={J_trial:.6e}"
            )

            # --- Objective
            tv = tv_after_prox
            Jvol = J_trial
            if has_floor:
                act_t = (floor_mask) & (z_trial < z_floor)
                R = 0.5 * float(np.sum((z_floor[act_t] - z_trial[act_t]) ** 2))
            else:
                R = 0.0
            obj = w_tv * tv + w_l2 * Jvol + w_reg * R
            rel_drop = (obj_last - obj) / (abs(obj_last) + 1e-16)

            # --- Step magnitudes (clear labels)
            step_L2 = float(np.max(np.abs(z_lin - z)))
            step_TV = float(np.max(np.abs(z_trial - z_lin)))
            step_total = float(np.max(np.abs(z_trial - z)))
            p90_total = float(np.percentile(np.abs(z_trial - z), 90.0))
            logger.debug(
                f"[iter {k:03d}] steps: L2={step_L2:.3f}m, TV={step_TV:.3f}m, total={step_total:.3f}m, p90_total={p90_total:.3f}m"
            )

            logger.info(
                f"[iter {k:03d}] Obj={obj:.6e} (w*TV={w_tv*tv:.6e}, w*L2={w_l2*Jvol:.6e})  rel_drop={rel_drop:.3e} (target {rel_eta:.3e})"
            )

            if rel_drop >= rel_eta:
                z = z_trial
                tv_last = tv
                Jvol_last = Jvol
                obj_last = obj
                accepted = True
                logger.debug(f"[iter {k:03d}] ACCEPT")
                dz_cum = z - z0
                p50, p90c, p99, p100 = np.percentile(np.abs(dz_cum), [50, 90, 99, 100])
                logger.info(
                    f"[Δz cumulative] P50={p50:.3f} P90={p90c:.3f} P99={p99:.3f} max={p100:.3f} m"
                )
                break

            if bt < max_backtracks:
                bt += 1
                dt_try *= 0.5
                logger.debug(
                    f"[iter {k:03d}] backtrack {bt}/{max_backtracks}: dt -> {dt_try:.3e}"
                )
                continue
            else:
                logger.info(
                    f"[iter {k:03d}] STOP: no improvement after {bt} backtrack(s)"
                )
                break

        history.append(z.copy())
        if not accepted:
            break

    return {"z": z, "history": np.stack(history, axis=0)}


# ------------------------------ Floor utilities ------------------------------
def build_shore_floor_from_df(
    df: pd.DataFrame,
    n_nodes: int,
    *,
    filter_kind: str = "max",  # or "median"
    window: int = 5,
    trim_endpoints: bool = True,
):
    """
    Construct a nodewise shoreline floor from an in-memory DataFrame.

    This mirrors the CSV-based implementation but avoids any file I/O.

    Parameters
    ----------
    df : pandas.DataFrame
        Columns: ``node, x, y, z, poly_id, label`` (``node`` is 1-based).
    n_nodes : int
        Total number of mesh nodes.
    filter_kind : {'max','median'}, optional
        Rolling aggregation applied along each polyline, by default ``'max'``.
    window : int, optional
        Odd window size for the rolling aggregation, by default ``7``.
    trim_endpoints : bool, optional
        Drop endpoints that touch any ``'deep'`` edge, by default ``True``.

    Returns
    -------
    z_floor : numpy.ndarray
        Floor elevation per node (length = ``n_nodes``), ``-inf`` where unset.
    mask : numpy.ndarray of bool
        True where a floor value is defined.
    """

    df = df.copy()
    df["node0"] = df["node"].astype(int) - 1
    deep_nodes = set(df.loc[df["label"] == "deep", "node0"].tolist())
    df_ss = df[df["label"] == "shore"]

    z_floor = np.full(int(n_nodes), -np.inf, dtype=float)
    k = max(1, (int(window) // 2) * 2 + 1)
    half = k // 2

    for _, g in df_ss.groupby("poly_id"):
        nodes = g["node0"].to_numpy()
        zc = g["z"].to_numpy()
        if trim_endpoints:
            while len(nodes) >= 2 and nodes[0] in deep_nodes:
                nodes, zc = nodes[1:], zc[1:]
            while len(nodes) >= 2 and nodes[-1] in deep_nodes:
                nodes, zc = nodes[:-1], zc[:-1]
        if len(nodes) == 0:
            continue
        if len(nodes) == 1:
            zf = zc.copy()
        else:
            zf = np.empty_like(zc)
            for i in range(len(zc)):
                a, b = max(0, i - half), min(len(zc), i + half + 1)
                win = zc[a:b]
                zf[i] = np.median(win) if filter_kind == "median" else np.max(win)
        for n, v in zip(nodes, zf):
            if v > z_floor[n]:
                z_floor[n] = v
    mask = np.isfinite(z_floor)
    return z_floor, mask


# -------------------------- Shoreline detection API --------------------------


@dataclass
class ShorelineOptions:
    # href can be a float or a path to a gr3 with elevations in nodes[:,2]
    href: Union[float, str] = 0.0
    deep_delta: float = 1.0
    shore_delta: float = 3.0
    # Either provide explicit seeds or use a built-in default set
    seeds: Optional[List[Tuple[float, float]]] = None
    use_default_seeds: bool = True
    # Stage 3 smoothing controls
    smooth_relax_factor: float = 0.4
    smooth_strip_max: int = 24
    smooth_eps_deg: float = 55.0
    # Output / filtering
    filter_deep: bool = True
    epsg: int = 26910  # default UTM10N NAD83
    shore_csv: Optional[str] = None
    shore_shp: Optional[str] = None


def _load_href_to_nodes(mesh: "SCHISM_mesh", href_opt: Union[float, str]) -> np.ndarray:
    """Load href (float or gr3 path) into a nodewise array."""
    try:
        v = float(href_opt)
        return np.full(mesh.n_nodes(), v, dtype=float)
    except Exception:
        m2 = read_mesh(str(href_opt))
        h = m2.nodes[:, 2].astype(float)
        if h.shape[0] != mesh.n_nodes():
            raise ValueError("Href mesh node count != target mesh node count")
        return h


@dataclass
class ShorelineData:
    """
    Container for shoreline detection outputs (pure data; no filenames).

    Attributes
    ----------
    polylines : list of dict
        One dict per merged polyline with keys:
        ``{'poly_id','label','nodes','coords'}``.
    df : pandas.DataFrame
        Per-node table equivalent to the shoreline CSV with columns
        ``['node','x','y','z','poly_id','label']`` (``node`` is 1-based).
    labels_perimeter : dict[int, str]
        Perimeter edge labels (``'shore'`` or ``'deep'``) keyed by edge id.
    wet_mask : numpy.ndarray (bool)
        Mask of elements in the final always-wet region after stage 2/3.
    """

    polylines: List[dict]
    df: pd.DataFrame
    labels_perimeter: Dict[int, str]
    wet_mask: np.ndarray


def detect_shorelines(
    mesh,
    *,
    href: float,
    deep_delta: float,
    shore_delta: float,
    seeds=None,
    use_default_seeds: bool = True,
    smooth_relax_factor: float = 1.25,
    smooth_strip_max: int = 12,
    smooth_eps_deg: float = 5.0,
    filter_deep: bool = True,
    epsg: int = 26910,
    shore_csv: Optional[str] = None,
    shore_shp: Optional[str] = None,
    logger=None,
) -> ShorelineData:
    """
    Detect perimeter shore/deep edges (stages 1–3) and return **data only**.

    Optionally writes diagnostic artifacts if ``shore_csv`` / ``shore_shp`` are
    provided. No filenames are returned.

    Parameters
    ----------
    mesh : schimpy.schism_mesh.SchismMesh
        Target mesh (``nodes`` depths are read but not mutated here).
    href : float
        Reference depth for nodewise tests (positive downward).
    deep_delta, shore_delta : float
        Threshold deltas (m) for deep vs. shore classification.
    seeds : sequence of (x, y), optional
        Seed points for the always-wet floodfill.
    use_default_seeds : bool, optional
        If True, augment with built-in seeds.
    smooth_relax_factor : float, optional
        Stage-3 relaxed edge criterion factor.
    smooth_strip_max : int, optional
        Maximum length of accepted stage-3 strips.
    smooth_eps_deg : float, optional
        Corner improvement tolerance in degrees.
    filter_deep : bool, optional
        Exclude rows labeled ``'deep'`` from CSV export.
    epsg : int, optional
        EPSG code for shapefile/CSV header.
    shore_csv, shore_shp : str or None, optional
        Optional paths to write diagnostics. If None, nothing is written.
    logger : logging.Logger, optional
        Logger for INFO/DEBUG emissions.

    Returns
    -------
    ShorelineData
        Polylines, per-node DataFrame, perimeter labels, and wet mask.
    """

    h0_node = _load_href_to_nodes(mesh, href)
    h0_edge = _edge_field_from_nodes(mesh, h0_node)

    # Seeds
    if (seeds is None) and (not use_default_seeds):
        raise ValueError("No seed coordinates and use_default_seeds is False.")
    seeds_xy = (
        list(seeds)
        if (seeds is not None)
        else [
            (548110.0, 4186370.0),
            (606258.0, 4213690.0),
            (626358.0, 4213690.0),
            (634475.0, 4186958.0),
        ]
    )

    # Stage 1-2 flood
    res = floodfill_always_wet(
        mesh, seeds_xy, h0_node, h0_edge, float(deep_delta), float(shore_delta)
    )
    in_after_stage2 = res.wet_elems_stage2
    perim_before = _perimeter_edges_of_region(mesh, in_after_stage2)

    # Stage 3 smoothing
    in_after, labels_final, n_strips, n_elems_added = stage3_single_pass_entry_exit(
        mesh,
        in_after_stage2,
        _build_elem_edges(mesh),
        _edge_mean_depth(mesh, mesh.nodes[:, 2].astype(float)),
        h0_node,
        h0_edge,
        float(deep_delta),
        float(shore_delta),
        float(smooth_relax_factor),
        int(smooth_strip_max),
        float(smooth_eps_deg),
        debug_nodes=None,
    )
    perim_final = set(_perimeter_edges_of_region(mesh, in_after))

    # Final labels + 'smoothed' where removed
    labels_all = labels_final  # from stage 3 (or 2), as you do today
    polylines, node_rows = _merged_polylines_and_node_rows(mesh, labels_all)

    # Build DataFrame for downstream use

    polylines, node_rows = _merged_polylines_and_node_rows(mesh, labels_all)
    if shore_csv:
        rows_for_csv = (
            node_rows
            if not filter_deep
            else [r for r in node_rows if str(r.get("label")) != "deep"]
        )
        _write_nodes_csv(rows_for_csv, shore_csv, epsg=epsg)

    if shore_csv:
        _write_nodes_csv(rows_for_csv, shore_csv, epsg=int(epsg))

    if shore_shp:
        _write_merged_shapefile(polylines, shore_shp, epsg=epsg)

    return ShorelineData(
        polylines=polylines,
        df=pd.DataFrame(node_rows),
        labels_perimeter=labels_all,
        wet_mask=in_after,
    )


# ----------------------------- Main driver (API) -----------------------------


@dataclass
class FloorOptions:
    filter: str = "max"  # "max" or "median"
    window: int = 5  # odd integer
    enforce_floor: bool = False  # project after TV prox
    reg_weight: float = 0.0  # hinge weight; 0 disables


@dataclass
class TVDOptions:
    # Mirrors mesh_refine_example defaults (sensible starting points)
    steps: int = 32
    dt: float = 1.0
    mu: float = 1.0
    lambda_l2: float = 500.0
    tv_weight: float = 1.0
    l2_weight: float = 5e-5
    cfl_target: float = 1.0
    clip_eps: Optional[float] = None
    rel_eta: float = 2e-3
    max_backtracks: int = 4


def refine_volume_tvd(
    mesh: "SCHISM_mesh",
    dem_spec: Union[str, Sequence[str], Mapping[str, Any]],
    out_dir: str,
    shoreline: Optional[ShorelineOptions] = None,
    floor: Optional[FloorOptions] = None,
    tvd: Optional[TVDOptions] = None,
    cache_dir: str = "./.dem_cache",
    profiles: Optional[Sequence[int]] = None,
    profiles_dir: Optional[str] = None,
    logger: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Execute the end-to-end shoreline+TVD workflow on an already loaded mesh.

    Parameters
    ----------
    mesh : SCHISM_mesh
        Mesh to refine. Boundary info is preserved.
    dem_spec : str | sequence | mapping
        DEM YAML path or object (anything :func:`create_dem_sampler` accepts).
    out_dir : str
        Where to write diagnostics (shore CSV/shape, profiles).
    shoreline : ShorelineOptions, optional
        If provided, shorelines are detected and used to build a floor.
    floor : FloorOptions, optional
        Controls floor regularization and enforcement.
    tvd : TVDOptions, optional
        Algorithm controls; sensible defaults used if None.
    cache_dir : str
        Disk cache for DEM sampler.
    profiles : sequence of int, optional
        poly_id list for which to dump triple profiles.
    profiles_dir : str, optional
        Directory where profile PNGs go (defaults to ``out_dir``).
    logger : logging.Logger, optional
        Logger instance for logging progress and messages.

    Returns
    -------
    dict
        Keys include:
        - "z": final node elevations
        - "shore_csv": CSV path if shorelines detected
        - "z0": initial DEM point samples
    """
    if logger is None:
        logger = logging.getLogger("mesh_volume_tvd")
        if not logger.handlers:
            h = logging.StreamHandler()
            h.setFormatter(logging.Formatter("%(message)s"))
            logger.addHandler(h)
        logger.setLevel(logging.INFO)

    logger.info("Starting refine_volume_tvd...")

    os.makedirs(out_dir, exist_ok=True)

    logger.info("Loading points from DEM sampler")

    # 1) DEM sampler and z0 from point samples (safer start)
    dem_sampler = create_dem_sampler(dem_spec, cache_dir=cache_dir)

    # 2) Optional shoreline detection → floor
    shore_csv = None
    z_floor = None
    floor_mask = None
    shore = None
    shore = None
    if shoreline is not None:
        shore = detect_shorelines(
            mesh,
            href=shoreline.href,
            deep_delta=float(shoreline.deep_delta),
            shore_delta=float(shoreline.shore_delta),
            seeds=shoreline.seeds,
            use_default_seeds=bool(shoreline.use_default_seeds),
            smooth_relax_factor=float(shoreline.smooth_relax_factor),
            smooth_strip_max=int(shoreline.smooth_strip_max),
            smooth_eps_deg=float(shoreline.smooth_eps_deg),
            filter_deep=bool(shoreline.filter_deep),
            epsg=int(shoreline.epsg),
            # Only if you want diagnostics written by the shoreline stage:
            shore_csv=shoreline.shore_csv,  # may be None
            shore_shp=shoreline.shore_shp,  # may be None
            logger=logger,
        )

        # Floor straight from memory (df), no CSV round-trip
        fopts = floor or FloorOptions()
        z_floor, floor_mask = build_shore_floor_from_df(
            shore.df,
            mesh.n_nodes(),
            filter_kind=fopts.filter,
            window=int(fopts.window),
            trim_endpoints=True,
        )

    # Optional profiles: plot from in-memory data (no CSV read)
    if (shore is not None) and profiles:
        out_dir_profiles = profiles_dir or out_dir
        os.makedirs(out_dir_profiles, exist_ok=True)
        _plot_profiles(
            shore.polylines,
            shore.df.to_dict("records"),
            set(int(p) for p in profiles),
            out_dir_profiles,
            filt=(floor.filter if floor else "none"),
            window=(floor.window if floor else 5),
        )

    # Build shoreline “floor” directly from in-memory data
    if shore is not None and floor is not None:
        z_floor, floor_mask = build_shore_floor_from_df(
            shore.df,
            mesh.n_nodes(),
            filter_kind=fopts.filter,
            window=int(fopts.window),
            trim_endpoints=True,
        )

    # 3) Build options and run TVD
    top = tvd or TVDOptions()
    opts = FlowOptions(
        n_steps=int(top.steps),
        dt=float(top.dt),
        mu=float(top.mu),
        lambda_l2=float(top.lambda_l2),
        tv_weight=float(top.tv_weight),
        l2_weight=float(top.l2_weight),
        cfl_target=float(top.cfl_target),
        clipping_eps=None if top.clip_eps is None else float(top.clip_eps),
        rel_eta=float(top.rel_eta),
        max_backtracks=int(top.max_backtracks),
    )

    reg_weight = 0.0
    enforce_floor = False
    if floor is not None:
        reg_weight = float(floor.reg_weight)
        enforce_floor = bool(floor.enforce_floor)

    logger.info(
        "Running Volume-TVD variational problem to improve volume fidelity, reduce noise and enforce shoreline constraints"
    )
    res = run_tvd_variational(
        mesh,
        z0=mesh.nodes[:, 2],
        dem_sampler=dem_sampler,
        options=opts,
        z_floor=z_floor,
        floor_mask=floor_mask,
        reg_weight=float(reg_weight),
        enforce_floor=bool(enforce_floor),
        logger=logger,
    )

    z = np.asarray(res["z"], dtype=float)
    mesh.nodes[:, 2] = z

    # Optional triple profiles
    if (shoreline is not None) and (profiles is not None and len(profiles) > 0):
        try:

            # Reuse the function from the example without importing to avoid circularity
            # Compute distances and dump a quick triple for selected ids
            df = pd.read_csv(shore_csv, comment="#")
            df["node0"] = df["node"].astype(int) - 1
            df_ss = df[df["label"].isin(["shore", "smoothed"])]
            out_dir_profiles = profiles_dir or out_dir
            os.makedirs(out_dir_profiles, exist_ok=True)
            for pid in profiles:
                g = df_ss[df_ss["poly_id"] == int(pid)]
                if g.empty:
                    continue
                nodes = g["node0"].to_numpy()
                xs = g["x"].to_numpy()
                ys = g["y"].to_numpy()
                if len(nodes) < 1:
                    continue
                if len(nodes) == 1:
                    s = np.array([0.0])
                else:
                    dxy = np.sqrt(np.diff(xs) ** 2 + np.diff(ys) ** 2)
                    s = np.concatenate([[0.0], np.cumsum(dxy)])
                z_orig = z0[nodes]
                z_fin = z[nodes]
                # Try to plot a filtered floor if available
                z_flt = z_floor[nodes] if (z_floor is not None) else None

                fig, ax = plt.subplots(figsize=(8, 3.0))
                (base_line,) = ax.plot(s, z_orig, label=f"original")
                if z_flt is not None:
                    ax.plot(
                        s,
                        z_flt,
                        linestyle="--",
                        label="filtered",
                        color=base_line.get_color(),
                    )
                ax.plot(s, z_fin, label="final")
                ax.set_xlabel("distance (m)")
                ax.set_ylabel("depth (m, +down)")
                ax.invert_yaxis()
                ax.set_title(f"poly_id={int(pid)}")
                ax.legend()
                fig.tight_layout()
                out_png = os.path.join(out_dir_profiles, f"profile_poly_{int(pid)}.png")
                fig.savefig(out_png, dpi=140)
                plt.close(fig)
        except Exception:
            pass
    logger.info("Mesh optimization complete")
    return res


# ------------------------------------ CLI ------------------------------------


def _cli():

    ap = argparse.ArgumentParser(
        description="TVD-regularized volume tuning (shoreline + floor + TVD)."
    )
    ap.add_argument("--mesh", required=True, help="Input mesh file (hgrid.gr3)")
    ap.add_argument("--dem", required=True, help="DEM YAML (file or object in YAML)")
    ap.add_argument("--out", required=True, help="Output mesh file")
    ap.add_argument("--cache", default="./.dem_cache", help="Disk cache directory")
    ap.add_argument(
        "--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )

    # TVD controls
    ap.add_argument("--steps", type=int, default=16)
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--mu", type=float, default=1.0)
    ap.add_argument("--lambda_l2", type=float, default=1000.0)
    ap.add_argument("--tv_weight", type=float, default=1.0)
    ap.add_argument("--l2_weight", type=float, default=2.5e-5)
    ap.add_argument("--cfl_target", type=float, default=1.0)
    ap.add_argument("--clip_eps", type=float, default=None)

    # Shoreline detection
    ap.add_argument(
        "--href",
        default=-0.35,
        help="Reference *depth* (not elev) for determining the 'always wet' contour (float or gr3 path).",
    )
    ap.add_argument(
        "--deep_delta",
        type=float,
        default=0.35,
        help="Additional depth from href for a shore to be classified 'deep' for purposes of regularization",
    )
    ap.add_argument(
        "--shore_delta",
        type=float,
        default=0.35,
        help="Additional distance above href for a shore to still be classified 'shallow' for purposes of regularization",
    )
    ap.add_argument("--use_default_seeds", action="store_true")
    ap.add_argument("--seed", nargs=2, type=float, action="append", default=None)
    ap.add_argument("--smooth_relax_factor", type=float, default=1.2)
    ap.add_argument("--smooth_strip_max", type=int, default=24)
    ap.add_argument("--smooth_eps_deg", type=float, default=55.0)
    ap.add_argument("--filter_deep", action="store_true")  # todo: necessary?
    ap.add_argument(
        "--epsg", type=int, default=26910
    )  # todo: shouldn't be unique to this process
    ap.add_argument(
        "--shore_csv", default=None, help="Output CSV file for node diagnostics"
    )
    ap.add_argument(
        "--shore_shp", default=None, help="Output Shapefile for shoreline diagnostics"
    )

    # Floor regularization
    ap.add_argument("--shore_floor_filter", choices=["max", "median"], default="median")
    ap.add_argument("--shore_floor_window", type=int, default=5)
    ap.add_argument("--reg_weight", type=float, default=8.0e-2)
    ap.add_argument("--enforce_floor", action="store_true")

    # Profiles (optional)
    ap.add_argument(
        "--profiles", nargs="*", type=int, default=None, help="poly_id(s) to plot"
    )
    ap.add_argument("--profiles_dir", default=None, help="Directory for profile PNGs")

    args = ap.parse_args()

    # Logger for CLI
    logger = logging.getLogger("mesh_volume_tvd")
    logger.handlers.clear()
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(h)

    logger.setLevel(getattr(logging, args.log_level.upper(), logging.INFO))
    mesh = read_mesh(args.mesh, nodestring_option="land")

    # Prepopulate the mesh based on point samples at nodes
    dem_sampler = create_dem_sampler(args.dem, cache_dir=args.cache)
    z0 = dem_sampler(mesh.nodes[:, :2]).astype(float)
    mesh.nodes[:, 2] = z0

    tvd = TVDOptions(
        steps=args.steps,
        dt=args.dt,
        mu=args.mu,
        lambda_l2=args.lambda_l2,
        tv_weight=args.tv_weight,
        l2_weight=args.l2_weight,
        cfl_target=args.cfl_target,
        clip_eps=args.clip_eps,
    )
    shoreline = ShorelineOptions(
        href=args.href,
        deep_delta=args.deep_delta,
        shore_delta=args.shore_delta,
        seeds=args.seed,
        use_default_seeds=bool(args.use_default_seeds),
        smooth_relax_factor=args.smooth_relax_factor,
        smooth_strip_max=args.smooth_strip_max,
        smooth_eps_deg=args.smooth_eps_deg,
        filter_deep=bool(args.filter_deep),
        epsg=int(args.epsg),
        shore_csv=args.shore_csv or None,
        shore_shp=args.shore_shp or None,
    )
    floor = FloorOptions(
        filter=args.shore_floor_filter,
        window=int(args.shore_floor_window),
        enforce_floor=bool(args.enforce_floor),
        reg_weight=float(args.reg_weight),
    )
    res = refine_volume_tvd(
        mesh,
        dem_spec=args.dem,
        out_dir=os.path.dirname(args.out) or ".",
        shoreline=shoreline,
        floor=floor,
        tvd=tvd,
        cache_dir=args.cache,
        profiles=args.profiles,
        profiles_dir=args.profiles_dir,
        logger=logger,
    )

    mesh.nodes[:, 2] = res["z"]

    write_mesh(mesh, args.out)
    logger.info(f"Wrote {args.out}")


if __name__ == "__main__":  # pragma: no cover
    _cli()
