"""Tests for hydraulic structure node path validation."""

import numpy as np
import pytest

from schimpy.schism_mesh import read_mesh
from schimpy.schism_setup import create_schism_setup


# 3x3 node grid; nodes are numbered row by row from the origin
NODES = [
    (0.0, 0.0), (100.0, 0.0), (200.0, 0.0),
    (0.0, 100.0), (100.0, 100.0), (200.0, 100.0),
    (0.0, 200.0), (100.0, 200.0), (200.0, 200.0),
]

QUADS = [(1, 2, 5, 4), (2, 3, 6, 5), (4, 5, 8, 7), (5, 6, 9, 8)]

# same cells split on the up[i]-down[i+1] diagonal, which touches the lower-left node
TRIS_DIAG_A = [
    (1, 2, 5), (1, 5, 4), (2, 3, 6), (2, 6, 5),
    (4, 5, 8), (4, 8, 7), (5, 6, 9), (5, 9, 8),
]

# split on the other diagonal, so one triangle of each cell misses the lower-left node
TRIS_DIAG_B = [
    (1, 2, 4), (2, 5, 4), (2, 3, 5), (3, 6, 5),
    (4, 5, 7), (5, 8, 7), (5, 6, 8), (6, 9, 8),
]

# spans the full width at mid-height, so it cuts the three vertical edges
CROSSING = np.array([[-10.0, 50.0], [210.0, 50.0]])


def _write_gr3(path, elems):
    lines = ["test mesh", "{} {} ! # of elements and nodes".format(len(elems), len(NODES))]
    for i, (x, y) in enumerate(NODES, start=1):
        lines.append("{} {:.8f} {:.8f} {:.8f}".format(i, x, y, 5.0))
    for i, e in enumerate(elems, start=1):
        lines.append("{} {} {}".format(i, len(e), " ".join(str(n) for n in e)))
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def quad_setup(tmp_path):
    return create_schism_setup(str(_write_gr3(tmp_path / "quad.gr3", QUADS)))


@pytest.fixture
def tri_a_setup(tmp_path):
    return create_schism_setup(str(_write_gr3(tmp_path / "tri_a.gr3", TRIS_DIAG_A)))


@pytest.fixture
def tri_b_setup(tmp_path):
    return create_schism_setup(str(_write_gr3(tmp_path / "tri_b.gr3", TRIS_DIAG_B)))


def _paths(setup):
    return setup.mesh.find_two_neighboring_node_paths(CROSSING)


def test_quad_cells_validate(quad_setup):
    up, down = _paths(quad_setup)
    assert len(up) == len(down) == 3
    quad_setup._validate_node_paths("gate", up, down)


def test_triangle_pair_validates_either_diagonal(tri_a_setup, tri_b_setup):
    for setup in (tri_a_setup, tri_b_setup):
        up, down = _paths(setup)
        setup._validate_node_paths("gate", up, down)


def test_empty_paths_rejected(quad_setup):
    with pytest.raises(ValueError, match="did not cross any element edges"):
        quad_setup._validate_node_paths("gate", [], [])


def test_unequal_path_lengths_rejected(quad_setup):
    up, down = _paths(quad_setup)
    with pytest.raises(ValueError, match="must match one-to-one"):
        quad_setup._validate_node_paths("gate", up, down[:-1])


def test_unconnected_node_pair_rejected(quad_setup):
    # node 3 is not adjacent to node 4, so the pair is not an edge
    with pytest.raises(ValueError, match="is not an edge"):
        quad_setup._validate_node_paths("gate", [2], [3])


def test_broken_path_rejected(quad_setup):
    # nodes 1 and 3 are on the same row but not adjacent
    with pytest.raises(ValueError, match="path is broken"):
        quad_setup._validate_node_paths("gate", [0, 2], [3, 5])
