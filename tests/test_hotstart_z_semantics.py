"""Tests pinning the meaning of ``z`` in attribute expressions.

The gr3/polygon route and the hotstart ``simple_trend`` route both accept an
expression containing ``z``. For ``elevation`` they agree: ``z`` is ``dp``,
positive down. They diverge only for genuinely 3D variables, where
``simple_trend`` binds ``z`` to the ``build_z`` layer elevation instead.
"""

import numpy as np
import pandas as pd

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
