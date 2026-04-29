#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for schimpy.relocate_source_sink."""
import os
import sys
import importlib

import numpy as np
import pandas as pd
import pytest
import yaml

import schimpy.schism_mesh as schism_mesh
from schimpy import relocate_source_sink as rss


HERE = os.path.dirname(__file__)
MESH_PATH = os.path.join(HERE, "testdata", "testmesh.gr3")


def _load_mesh_with_depths(depth_fn):
    """Load testmesh.gr3 and override node depths via a callable depth_fn(x, y)."""
    reader = schism_mesh.SchismMeshIoFactory().get_reader("gr3")
    mesh = reader.read(MESH_PATH)
    nodes = mesh.nodes.copy()
    for i in range(nodes.shape[0]):
        nodes[i, 2] = depth_fn(nodes[i, 0], nodes[i, 1])
    mesh.nodes = nodes
    return mesh


def _write_pair(tmp_path, prefix, src_pts, sink_pts):
    src_fp = tmp_path / f"source_{prefix}.yaml"
    sink_fp = tmp_path / f"sink_{prefix}.yaml"
    with open(src_fp, "w") as fh:
        yaml.safe_dump({k: list(v) for k, v in src_pts.items()}, fh)
    with open(sink_fp, "w") as fh:
        yaml.safe_dump({k: list(v) for k, v in sink_pts.items()}, fh)
    return str(src_fp), str(sink_fp)


def test_no_geopandas_at_module_load():
    """The library module must not pull in the spatial stack."""
    # Remove from sys.modules to ensure a clean check on re-import.
    for key in list(sys.modules):
        if key.startswith("geopandas") or key.startswith("fiona"):
            sys.modules.pop(key, None)
    importlib.reload(rss)
    assert "geopandas" not in sys.modules
    assert "fiona" not in sys.modules


def test_parse_name():
    assert rss._parse_name("delta_src_5") == ("delta", "src", 5)
    assert rss._parse_name("suisun_sink_227") == ("suisun", "sink", 227)
    assert rss._parse_name("foo_src_1") is None
    assert rss._parse_name("delta_other_5") is None


def test_build_pairs_filters_other_prefix():
    df_src = pd.DataFrame(
        {"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]},
        index=["delta_src_1", "suisun_src_2", "delta_src_3"],
    )
    df_sink = pd.DataFrame(
        {"x": [1.0, 3.0, 4.0], "y": [1.0, 3.0, 4.0]},
        index=["delta_sink_1", "delta_sink_3", "delta_sink_99"],
    )
    pairs, skipped = rss.build_pairs(df_src, df_sink, "delta")
    suffixes = sorted(p[0] for p in pairs)
    assert suffixes == [1, 3]
    # suisun_src_2 should be skipped (prefix mismatch); delta_sink_99 unpaired.
    assert "suisun_src_2" in skipped
    assert "delta_sink_99" in skipped


def test_separated_pair_passthrough(tmp_path):
    mesh = _load_mesh_with_depths(lambda x, y: 5.0)
    rng = np.random.default_rng(0)
    res = rss.relocate_pair(mesh, 10.0, 50.0, 30.0, 50.0, rng)
    assert res["status"] == "separated"
    assert res["source"] == (10.0, 50.0)
    assert res["sink"] == (30.0, 50.0)
    assert res["snapped"] is False


def test_off_grid_point_snaps_to_mesh():
    """A co-located pair lying outside the mesh must snap to the nearest
    element rather than failing, and the result must record the snap."""
    mesh = _load_mesh_with_depths(lambda x, y: 5.0)
    # The simple_triquad testmesh covers roughly (0,0)-(100,100); use a point
    # well outside the mesh.
    rng = np.random.default_rng(0)
    res = rss.relocate_pair(mesh, 500.0, 500.0, 500.0, 500.0, rng)
    assert res["status"] == "relocated"
    assert res["snapped"] is True
    assert res["source_elem"] != res["sink_elem"]


def test_separated_pair_with_off_grid_endpoint_relocates():
    """A pair already separated beyond ``separation_tol`` but with a sink
    lying outside the mesh must project the sink onto the nearest qualifying
    element. The on-grid endpoint must pass through unchanged."""
    mesh = _load_mesh_with_depths(lambda x, y: 5.0)
    rng = np.random.default_rng(0)
    # source on-grid; sink off-grid (mesh covers ~(0,0)-(100,100)).
    res = rss.relocate_pair(mesh, 10.0, 50.0, 500.0, 500.0, rng)
    assert res["status"] == "relocated"
    assert res["snapped"] is True
    # On-grid source coordinate is preserved.
    assert res["source"] == (10.0, 50.0)
    # Off-grid sink has been moved onto the mesh.
    sx, sy = res["sink"]
    assert 0.0 <= sx <= 100.0 and 0.0 <= sy <= 100.0
    assert res["sink_elem"] is not None


def test_separated_pair_on_dry_mesh_fails_and_emits_polygons():
    """Regression for the swap-mesh scenario: a pair that was previously
    separated and on-grid in the *old* mesh may sit on dry/shallow elements
    in a *new* mesh. The driver must detect insufficient depth at each
    on-grid endpoint, fail the pair, and emit min-depth polygons so the
    user can dredge — not silently pass the points through."""
    # All-dry mesh: every element fails depth criteria.
    mesh = _load_mesh_with_depths(lambda x, y: -10.0)
    rng = np.random.default_rng(0)
    # Both endpoints on-grid and well separated.
    res = rss.relocate_pair(mesh, 10.0, 10.0, 90.0, 90.0, rng)
    assert res["status"] == "failed"
    assert res["snapped"] is False
    assert res["source_elem"] is not None
    assert res["sink_elem"] is not None
    # Polygons must be buildable for both endpoints.
    poly_src = rss._build_min_depth_polygon(
        mesh, res["source_elem"], "delta_src_1", res["source_mean"]
    )
    poly_sink = rss._build_min_depth_polygon(
        mesh, res["sink_elem"], "delta_sink_1", res["sink_mean"]
    )
    assert float(poly_src.attribute) >= 0.5
    assert float(poly_sink.attribute) >= 0.5


def test_relocated_when_co_located_and_deep(tmp_path):
    # Uniform deep water -> qualifying everywhere; should relocate, not fail.
    mesh = _load_mesh_with_depths(lambda x, y: 5.0)
    rng = np.random.default_rng(0)
    res = rss.relocate_pair(mesh, 50.0, 50.0, 50.0, 50.0, rng)
    assert res["status"] == "relocated"
    # Both must be moved to distinct elements.
    assert res["source_elem"] != res["sink_elem"]
    # Output must be different from input (both moved).
    dx_s = res["source"][0] - 50.0
    dy_s = res["source"][1] - 50.0
    dx_n = res["sink"][0] - 50.0
    dy_n = res["sink"][1] - 50.0
    assert (dx_s, dy_s) != (0.0, 0.0) or (dx_n, dy_n) != (0.0, 0.0)


def test_failure_emits_polygon(tmp_path):
    # All-dry mesh: depth = -10 everywhere -> no element qualifies.
    mesh = _load_mesh_with_depths(lambda x, y: -10.0)
    rng = np.random.default_rng(0)
    res = rss.relocate_pair(mesh, 50.0, 50.0, 50.0, 50.0, rng)
    assert res["status"] == "failed"
    # Source and sink must end up in different elements even on failure.
    assert res["source_elem"] != res["sink_elem"]
    assert res["source"] != res["sink"]
    # Sink should be at least as deep as source (deeper -> sink rule).
    assert res["sink_mean"] >= res["source_mean"]
    poly_src = rss._build_min_depth_polygon(
        mesh, res["source_elem"], "delta_src_1", res["source_mean"]
    )
    poly_sink = rss._build_min_depth_polygon(
        mesh, res["sink_elem"], "delta_sink_1", res["sink_mean"]
    )
    assert float(poly_src.attribute) >= 0.5
    assert float(poly_sink.attribute) >= 0.5


def test_relocated_pair_distinct_elements():
    """Relocated source and sink must always live in different elements."""
    mesh = _load_mesh_with_depths(lambda x, y: 5.0)
    for s in range(20):
        rng = np.random.default_rng(s)
        res = rss.relocate_pair(mesh, 50.0, 50.0, 50.0, 50.0, rng)
        assert res["status"] == "relocated"
        assert res["source_elem"] != res["sink_elem"]
        assert res["source"] != res["sink"]


def test_distinct_when_only_seed_qualifies():
    """If only the seed element qualifies, the pair must escalate to the next
    ring radius rather than collapsing source and sink onto the same element.
    """
    # Depth: only the central element qualifies; everything else is dry.
    # Strategy: deep at exactly the seed centroid area, dry elsewhere.
    reader = schism_mesh.SchismMeshIoFactory().get_reader("gr3")
    mesh = reader.read(MESH_PATH)
    seed, _snapped = rss._seed_element(mesh, 50.0, 50.0)
    seed_nodes = set(mesh.elem(seed).tolist())
    nodes = mesh.nodes.copy()
    nodes[:, 2] = -10.0
    for i in seed_nodes:
        nodes[i, 2] = 5.0
    mesh.nodes = nodes
    rng = np.random.default_rng(0)
    res = rss.relocate_pair(mesh, 50.0, 50.0, 50.0, 50.0, rng)
    if res["status"] == "relocated":
        assert res["source_elem"] != res["sink_elem"]


def test_end_to_end_writes_outputs(tmp_path):
    """Smoke test the full driver: relocated + separated + skipped."""
    mesh = _load_mesh_with_depths(lambda x, y: 5.0)
    import schimpy.relocate_source_sink as mod
    orig = mod.read_mesh
    mod.read_mesh = lambda *a, **k: mesh
    try:
        src_fp, sink_fp = _write_pair(
            tmp_path,
            "delta",
            {
                "delta_src_1": [70.0, 50.0],   # co-located -> relocate
                "delta_src_3": [70.0, 30.0],   # already separated
                "extra_src_4": [70.0, 30.0],   # bad name -> skipped
            },
            {
                "delta_sink_1": [70.0, 50.0],
                "delta_sink_3": [80.0, 30.0],
                "extra_sink_4": [70.0, 30.0],
            },
        )
        result = rss.relocate_source_sink(
            hgrid="ignored",
            source_yaml=src_fp,
            sink_yaml=sink_fp,
            prefix="delta",
            out_dir=str(tmp_path),
            seed=1,
        )
    finally:
        mod.read_mesh = orig

    assert result["counts"]["relocated"] == 1
    assert result["counts"]["separated"] == 1
    assert result["counts"]["failed"] == 0
    assert result["polygon_yaml"] is None
    assert os.path.exists(result["source_yaml"])
    assert os.path.exists(result["sink_yaml"])
    assert os.path.exists(result["skipped_yaml"])

    # Validate that the relocated source for pair 1 differs from the sink
    out_src = rss.read_points_yaml(result["source_yaml"], "sources")
    out_sink = rss.read_points_yaml(result["sink_yaml"], "sinks")
    assert (out_src.loc["delta_src_1", "x"], out_src.loc["delta_src_1", "y"]) != \
           (out_sink.loc["delta_sink_1", "x"], out_sink.loc["delta_sink_1", "y"])

    # Output filenames are derived from the input basenames with _adj appended.
    assert result["source_yaml"].endswith("_delta_adj.yaml")
    assert result["sink_yaml"].endswith("_delta_adj.yaml")
    with open(result["source_yaml"]) as fh:
        first_line = fh.readline()
    assert first_line.startswith("# Adjusted using relocate_source_sink")
