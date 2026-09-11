"""Tests for wet/dry flagging in the hotstart.

A node is dry when its water column ``H = dp + eta2`` falls to or below ``h0``.
Sides and elements are dry when any of their nodes is.
"""

import numpy as np
import xarray as xr

from schimpy.schism_hotstart import hotstart, prepare_schism_hotstart_output

H0 = 0.01

# levee node put 1 cm under its own bed, the other two left wet.
# dp = -3 means the bed is 3 m ABOVE datum, so "just below bed" is eta = -dp - 0.01 = +2.99
ETA = np.array([0.0, 2.99, 0.0])
EXPECTED_IDRY = np.array([0, 1, 0])


def _hotstart_with_eta(mesh, eta):
    h = hotstart(input=None)
    h.mesh = mesh
    h.depths = mesh.build_z()
    h.h0 = H0
    h.nc_dataset = xr.Dataset(
        {
            "elevation": ("node", np.asarray(eta, dtype=float)),
            "idry": ("node", np.zeros(mesh.n_nodes(), dtype=int)),
            "idry_s": ("side", np.zeros(mesh.n_edges(), dtype=int)),
            "idry_e": ("elem", np.zeros(mesh.n_elems(), dtype=int)),
        }
    )
    return h


def test_dry_node_flagged_from_eta(triangle_mesh, triangle_dp):
    np.testing.assert_allclose(triangle_dp + ETA, [2.0, -0.01, 0.05])

    h = _hotstart_with_eta(triangle_mesh, ETA)
    h.wet_dry_check()

    np.testing.assert_array_equal(h.nc_dataset["idry"].values, EXPECTED_IDRY)


def test_dry_propagates_to_side_and_element(triangle_mesh):
    h = _hotstart_with_eta(triangle_mesh, ETA)
    h.wet_dry_check()

    edges = triangle_mesh.edges[:, :2]
    expected_s = np.array([int(EXPECTED_IDRY[a] or EXPECTED_IDRY[b]) for a, b in edges])
    np.testing.assert_array_equal(h.nc_dataset["idry_s"].values, expected_s)
    # the single element touches the dry node
    np.testing.assert_array_equal(h.nc_dataset["idry_e"].values, [1])


def test_all_wet_leaves_everything_wet(triangle_mesh):
    h = _hotstart_with_eta(triangle_mesh, np.array([0.0, 5.0, 0.0]))
    h.wet_dry_check()

    np.testing.assert_array_equal(h.nc_dataset["idry"].values, [0, 0, 0])
    np.testing.assert_array_equal(h.nc_dataset["idry_e"].values, [0])


def test_elevation_initializer_flags_preserve_dry_but_not_impossible_wet(
    triangle_mesh,
):
    h = _hotstart_with_eta(triangle_mesh, ETA)
    h.elevation_idry = np.array([1, 0, -1])

    h.wet_dry_check()

    expected_idry = np.array([1, 1, 0])
    np.testing.assert_array_equal(h.nc_dataset["idry"].values, expected_idry)
    edges = triangle_mesh.edges[:, :2]
    expected_s = np.array([int(expected_idry[a] or expected_idry[b]) for a, b in edges])
    np.testing.assert_array_equal(h.nc_dataset["idry_s"].values, expected_s)
    np.testing.assert_array_equal(h.nc_dataset["idry_e"].values, [1])


def test_schism_output_uses_model_dtypes_and_omits_helpers(tmp_path):
    dataset = xr.Dataset(
        {
            "time": ("one", np.array([90], dtype=np.int64)),
            "iths": ("one", np.array([1], dtype=np.int64)),
            "ifile": ("one", np.array([1], dtype=np.int64)),
            "nsteps_from_cold": ("one", np.array([1], dtype=np.int64)),
            "idry": ("node", np.array([0, 1], dtype=np.int64)),
            "idry_s": ("side", np.array([0], dtype=np.int64)),
            "idry_e": ("elem", np.array([0], dtype=np.int64)),
            "su2": (("side", "nVert"), np.ones((1, 2), dtype=np.float32)),
            "sv2": (("side", "nVert"), np.ones((1, 2), dtype=np.float32)),
            "we": (("elem", "nVert"), np.ones((1, 2), dtype=np.float32)),
            "tr_nd": (
                ("node", "nVert", "ntracers"),
                np.ones((2, 2, 1), dtype=np.float32),
            ),
            "tr_nd0": (
                ("node", "nVert", "ntracers"),
                np.ones((2, 2, 1), dtype=np.float32),
            ),
            "z": (("node", "nVert"), np.zeros((2, 2), dtype=np.float32)),
        },
        coords={"tracer_list": ("ntracers", np.array(["TEM"], dtype="S10"))},
    )

    output = prepare_schism_hotstart_output(dataset)
    output_path = tmp_path / "hotstart.nc"
    output.to_netcdf(output_path)

    with xr.open_dataset(output_path) as written:
        for name in (
            "iths",
            "ifile",
            "nsteps_from_cold",
            "idry",
            "idry_s",
            "idry_e",
        ):
            assert written[name].dtype == np.int32
        for name in ("time", "su2", "sv2", "we", "tr_nd", "tr_nd0"):
            assert written[name].dtype == np.float64
        assert "z" not in written
        assert "tracer_list" not in written
        assert "string10" not in written.dims

    assert dataset["su2"].dtype == np.float32
    assert "z" in dataset
