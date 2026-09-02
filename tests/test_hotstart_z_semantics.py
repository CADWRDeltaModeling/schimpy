"""Tests pinning the meaning of ``z`` in attribute expressions.

The gr3/polygon route and the hotstart ``simple_trend`` route both accept an
expression containing ``z``. For ``elevation`` they agree: ``z`` is ``dp``,
positive down. They diverge only for genuinely 3D variables, where
``simple_trend`` binds ``z`` to the ``build_z`` layer elevation instead.
"""

import numpy as np
import pandas as pd
import pytest

from schimpy.schism_mesh import read_mesh
from schimpy.schism_setup import create_schism_setup
from schimpy.schism_vertical_mesh import SchismLocalVerticalMesh
from schimpy.schism_hotstart import VariableField


# Three nodes is the minimum that forms an element, and is enough to cover the
# three regimes: submerged, dry land, and straddling the max(0.1, .) clamp.
DP = np.array([2.0, -3.0, 0.05])  # depth positive down, so bed elev = -DP
NVRT = 3


@pytest.fixture
def gr3_path(tmp_path):
    path = tmp_path / "hgrid.gr3"
    xy = [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0)]
    lines = ["one triangle", "1 3 ! # of elements and nodes"]
    for i, ((x, y), dp) in enumerate(zip(xy, DP), start=1):
        lines.append(f"{i} {x:.8f} {y:.8f} {dp:.8f}")
    lines.append("1 3 1 2 3")
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def mesh(gr3_path):
    m = read_mesh(str(gr3_path))
    vmesh = SchismLocalVerticalMesh()
    # sigma runs -1 at the bottom to 0 at the surface, so column index 0 is the bed
    vmesh.init(np.tile(np.array([-1.0, -0.5, 0.0]), (3, 1)))
    vmesh.param["nvrt"] = NVRT
    m._vmesh = vmesh
    return m


def test_build_z_bottom_is_negated_and_clamped(mesh):
    """build_z returns elevations, not depths, and clamps the column to 0.1 m."""
    z = mesh.build_z()
    assert z.shape == (3, NVRT)
    np.testing.assert_allclose(z[:, 0], -np.maximum(0.1, DP))
    # the clamp erases terrain: the levee node and the near-datum node collapse
    assert z[1, 0] == z[2, 0]


def test_gr3_polygon_z_is_depth_positive_down(gr3_path):
    """The time-tested route: z is dp, so -z-0.01 puts the surface just under the bed."""
    s = create_schism_setup(str(gr3_path))
    polygons = [
        {
            "name": "all",
            "type": "none",
            "attribute": "-z-0.01",
            "vertices": [[-10.0, -10.0], [200.0, -10.0], [200.0, 200.0], [-10.0, 200.0]],
        }
    ]
    eta = s.apply_polygons(polygons=polygons, default=None)
    np.testing.assert_allclose(eta, -DP - 0.01)
    # water column is a uniform -1 cm across 5 m of relief
    np.testing.assert_allclose(DP + eta, -0.01)


def _field(mesh, vname, expr):
    v_meta = {"initializer": {"simple_trend": {"values": expr}}}
    return VariableField(
        v_meta,
        vname,
        mesh,
        mesh.build_z(),
        pd.Timestamp("2021-01-01"),
        "EPSG:26910",
        [],
        None,
    )


def test_simple_trend_elevation_z_agrees_with_gr3(mesh):
    """elevation is node2D, whose vgrid is node_z, so z is dp exactly as in gr3."""
    eta = np.asarray(_field(mesh, "elevation", "-z-0.01").simple_trend())
    np.testing.assert_allclose(eta, -DP - 0.01)
    np.testing.assert_allclose(DP + eta, -0.01)


def test_simple_trend_3d_variable_z_is_layer_elevation(mesh):
    """A 3D variable binds vgrid to build_z, so z there is clamped bed elevation."""
    vals = np.asarray(_field(mesh, "salinity", "-z-0.01").simple_trend())
    np.testing.assert_allclose(vals[:, 0], np.maximum(0.1, DP) - 0.01)
    # the 0.1 clamp makes the two shallow nodes indistinguishable
    assert vals[1, 0] == vals[2, 0]
