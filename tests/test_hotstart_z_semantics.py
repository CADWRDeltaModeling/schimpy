"""Tests pinning the meaning of ``z`` in attribute expressions.

The gr3/polygon route and the hotstart ``simple_trend`` route both accept an
expression containing ``z``. For ``elevation`` they agree: ``z`` is ``dp``,
positive down. They diverge only for genuinely 3D variables, where
``simple_trend`` binds ``z`` to the ``build_z`` layer elevation instead.
"""

import numpy as np
import pandas as pd
import pytest

from schimpy import geo_tools
from schimpy.schism_setup import create_schism_setup
from schimpy.schism_hotstart import VariableField


def test_build_z_bottom_is_negated_and_clamped(triangle_mesh, triangle_dp):
    """build_z returns elevations, not depths, and clamps the column to 0.1 m."""
    z = triangle_mesh.build_z()
    assert z.shape == (3, 3)
    np.testing.assert_allclose(z[:, 0], -np.maximum(0.1, triangle_dp))
    # the clamp erases terrain: the levee node and the near-datum node collapse
    assert z[1, 0] == z[2, 0]


def test_gr3_polygon_z_is_depth_positive_down(triangle_gr3, triangle_dp):
    """The time-tested route: z is dp, so -z-0.01 puts the surface just under the bed."""
    s = create_schism_setup(str(triangle_gr3))
    polygons = [
        {
            "name": "all",
            "type": "none",
            "attribute": "-z-0.01",
            "vertices": [[-10.0, -10.0], [200.0, -10.0], [200.0, 200.0], [-10.0, 200.0]],
        }
    ]
    eta = s.apply_polygons(polygons=polygons, default=None)
    np.testing.assert_allclose(eta, -triangle_dp - 0.01)
    # water column is a uniform -1 cm across 5 m of relief
    np.testing.assert_allclose(triangle_dp + eta, -0.01)


def _field(triangle_mesh, vname, expr):
    v_meta = {"initializer": {"simple_trend": {"values": expr}}}
    return VariableField(
        v_meta,
        vname,
        triangle_mesh,
        triangle_mesh.build_z(),
        pd.Timestamp("2021-01-01"),
        "EPSG:26910",
        [],
        None,
    )


def test_simple_trend_elevation_z_agrees_with_gr3(triangle_mesh, triangle_dp):
    """elevation is node2D, whose vgrid is node_z, so z is dp exactly as in gr3."""
    eta = np.asarray(_field(triangle_mesh, "elevation", "-z-0.01").simple_trend())
    np.testing.assert_allclose(eta, -triangle_dp - 0.01)
    np.testing.assert_allclose(triangle_dp + eta, -0.01)


def test_simple_trend_3d_variable_z_is_layer_elevation(triangle_mesh, triangle_dp):
    """A 3D variable binds vgrid to build_z, so z there is clamped bed elevation."""
    vals = np.asarray(_field(triangle_mesh, "salinity", "-z-0.01").simple_trend())
    np.testing.assert_allclose(vals[:, 0], np.maximum(0.1, triangle_dp) - 0.01)
    # the 0.1 clamp makes the two shallow nodes indistinguishable
    assert vals[1, 0] == vals[2, 0]


def test_patch_init_accepts_yaml_region_polygons(triangle_mesh, monkeypatch):
    """Generated inundation regions are YAML polygons consumed by partition_check."""
    field = _field(triangle_mesh, "elevation", "0.0")
    field.ini_meta = {
        "regions_filename": "inundate_regions.yaml",
        "regions": [],
        "smoothing": False,
    }

    def accepted(*args, **kwargs):
        raise RuntimeError("region YAML accepted")

    monkeypatch.setattr(geo_tools, "partition_check", accepted)
    with pytest.raises(RuntimeError, match="region YAML accepted"):
        field.patch_init()


def test_patch_init_preserves_domain_dtype(triangle_mesh, monkeypatch):
    field = _field(triangle_mesh, "velocity_u", "0.0")
    field.ini_meta = {
        "regions_filename": "regions.yaml",
        "smoothing": False,
        "regions": [
            {"region": "domain", "initializer": {"hotstart_nc": {}}},
            {"region": "restoration", "initializer": {"simple_trend": 0.0}},
        ],
    }
    mapping = np.full(field.n_hgrid, "domain", dtype=object)
    mapping[-1] = "restoration"
    monkeypatch.setattr(geo_tools, "partition_check", lambda *args, **kwargs: mapping)
    monkeypatch.setattr(
        field,
        "hotstart_nc",
        lambda ini_meta, inpoly: np.ones(
            (len(inpoly), field.n_vgrid), dtype=np.float32
        ),
    )

    values = field.patch_init()

    assert values.dtype == np.float32
    np.testing.assert_array_equal(values[:-1], 1.0)
    np.testing.assert_array_equal(values[-1], 0.0)
