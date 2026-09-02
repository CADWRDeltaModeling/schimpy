"""Tests for donor selection when interpolating from a source hotstart.

A target node with no counterpart in the source grid must not inherit values from
a dry source node, whose state is vestigial. Donor reselection is scoped to the
nodes handed to the initializer, so a region initialized by other means is never
touched.
"""

import numpy as np
import pytest

from schimpy.schism_hotstart import redirect_to_wet_donors, VariableField

SRC_XY = np.array([[0.0, 0.0], [10.0, 0.0], [20.0, 0.0], [30.0, 0.0]])
SRC_WET = np.array([True, False, False, True])
TOL = 0.5


def test_novel_node_takes_nearest_wet_donor():
    tgt_xy = np.array([[0.0, 0.0], [11.0, 0.0], [30.0, 0.0]])
    indices = np.array([0, 1, 3])
    dist = np.array([0.0, 1.0, 0.0])

    out, novel = redirect_to_wet_donors(SRC_XY, SRC_WET, tgt_xy, indices, dist, TOL)

    # only the middle node is novel; its raw donor was dry node 1, nearest wet is node 0
    np.testing.assert_array_equal(out, [0, 0, 3])
    np.testing.assert_array_equal(novel, [False, True, False])


def test_no_novel_nodes_leaves_donors_untouched():
    tgt_xy = np.array([[0.0, 0.0], [10.0, 0.0]])
    indices = np.array([0, 1])
    dist = np.array([0.0, 0.0])

    out, novel = redirect_to_wet_donors(SRC_XY, SRC_WET, tgt_xy, indices, dist, TOL)

    np.testing.assert_array_equal(out, [0, 1])
    assert not novel.any()


def test_donor_reselection_does_not_touch_non_novel_nodes():
    """A node with a counterpart keeps its donor even if that donor is dry."""
    tgt_xy = np.array([[10.0, 0.0], [11.0, 0.0]])
    indices = np.array([1, 1])
    dist = np.array([0.0, 1.0])

    out, _ = redirect_to_wet_donors(SRC_XY, SRC_WET, tgt_xy, indices, dist, TOL)

    assert out[0] == 1  # exact match on a dry source node is preserved
    assert out[1] == 0  # novel node moved to a wet donor


def test_no_wet_source_nodes_raises():
    tgt_xy = np.array([[11.0, 0.0]])
    with pytest.raises(ValueError, match="no wet nodes"):
        redirect_to_wet_donors(
            SRC_XY, np.zeros(4, dtype=bool), tgt_xy, np.array([1]), np.array([1.0]), TOL
        )


class _Field:
    """Minimal stand-in exposing only what the validation reads."""

    def __init__(self, variable_name):
        self.variable_name = variable_name

    _validated_max_blw_bed = VariableField._validated_max_blw_bed


def test_elevation_requires_max_blw_bed():
    with pytest.raises(ValueError, match="must set max_blw_bed"):
        _Field("elevation")._validated_max_blw_bed({"data_source": "x.nc"})


def test_elevation_rejects_negative_max_blw_bed():
    with pytest.raises(ValueError, match="non-negative"):
        _Field("elevation")._validated_max_blw_bed({"max_blw_bed": -1.0})


def test_max_blw_bed_rejected_on_other_variables():
    with pytest.raises(ValueError, match="cannot be applied to 'salinity'"):
        _Field("salinity")._validated_max_blw_bed({"max_blw_bed": 0.5})


def test_other_variables_need_no_floor():
    assert _Field("salinity")._validated_max_blw_bed({"data_source": "x.nc"}) is None


def test_elevation_floor_value_is_returned():
    assert _Field("elevation")._validated_max_blw_bed({"max_blw_bed": 0.5}) == 0.5
