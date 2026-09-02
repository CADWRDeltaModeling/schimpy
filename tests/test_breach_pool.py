"""Tests for pooling a gated breach without needing pre-dredge depths.

``ellipse`` clamps its radial term at 1, so every node inside a breach polygon
comes out at ``min_depth`` or deeper. That makes ``-min_depth - 0.01`` a constant
pool level that meets the island's ``-z-0.01`` exactly at the ellipse rim, and a
``type: min`` polygon applies it as a lower bound so ambient water and dry ground
are left alone.
"""

import numpy as np
import pytest

from schimpy.schism_setup import create_schism_setup

# levee runs along x = 0; island is x < 0, channel is x > 0
BREACH_LEFT = [0.0, -10.0]
BREACH_RIGHT = [0.0, 10.0]
MIN_DEPTH = -0.20
MAX_DEPTH = 4.09
MAJOR_AXIS_LEN = 30.0
POOL = -MIN_DEPTH - 0.01  # 0.19
AMBIENT = 0.96

# hole (ellipse centre), rim, channel, island interior
XY = [(0.0, 0.0), (0.0, 12.0), (15.0, 0.0), (-40.0, 5.0)]
DP0 = np.array([-1.0, -0.5, 3.0, 0.5])
HOLE, RIM, CHANNEL, ISLAND = 0, 1, 2, 3

BREACH_VERTICES = [[-20.0, -15.0], [20.0, -15.0], [20.0, 15.0], [-20.0, 15.0]]
ISLAND_VERTICES = [[-50.0, -20.0], [5.0, -20.0], [5.0, 20.0], [-50.0, 20.0]]
DOMAIN_VERTICES = [[-100.0, -100.0], [100.0, -100.0], [100.0, 100.0], [-100.0, 100.0]]

IMPORTS = ["schimpy.ellipse.ellipse"]
DREDGE_ATTR = (
    f"ellipse(x, y, z, {BREACH_LEFT}, {BREACH_RIGHT}, "
    f"min_depth={MIN_DEPTH}, max_depth={MAX_DEPTH}, major_axis_len={MAJOR_AXIS_LEN})"
)


@pytest.fixture
def breach_setup(tmp_path):
    path = tmp_path / "hgrid.gr3"
    lines = ["breach test", "2 4 ! # of elements and nodes"]
    for i, ((x, y), dp) in enumerate(zip(XY, DP0), start=1):
        lines.append(f"{i} {x:.8f} {y:.8f} {dp:.8f}")
    lines.append("1 3 1 3 2")
    lines.append("2 3 1 2 4")
    path.write_text("\n".join(lines) + "\n")
    return create_schism_setup(str(path))


def _dredge(setup):
    polygons = [
        {
            "name": "breach",
            "type": "none",
            "attribute": DREDGE_ATTR,
            "vertices": BREACH_VERTICES,
        }
    ]
    setup.mesh.nodes[:, 2] = setup.apply_polygons(
        polygons=polygons, default=None, global_imports=IMPORTS
    )
    return setup.mesh.nodes[:, 2]


def test_ellipse_floors_whole_polygon_at_min_depth(breach_setup):
    z = _dredge(breach_setup)

    assert z[HOLE] == pytest.approx(MAX_DEPTH)
    assert z[RIM] == pytest.approx(MIN_DEPTH)
    # already deeper than the rim, so the dredge leaves the channel alone
    assert z[CHANNEL] == pytest.approx(DP0[CHANNEL])
    # outside the polygon entirely
    assert z[ISLAND] == pytest.approx(DP0[ISLAND])

    in_polygon = [HOLE, RIM, CHANNEL]
    assert np.all(z[in_polygon] >= MIN_DEPTH - 1e-9)


def test_breach_pool_is_a_lower_bound_over_island_and_ambient(breach_setup):
    z = _dredge(breach_setup)

    polygons = [
        {
            "name": "domain",
            "type": "none",
            "attribute": f"max({AMBIENT}, -z-0.01)",
            "vertices": DOMAIN_VERTICES,
        },
        {
            "name": "island",
            "type": "none",
            "attribute": "-z-0.01",
            "vertices": ISLAND_VERTICES,
        },
        # order matters: the breach lower bound must follow the island
        {
            "name": "breach",
            "type": "min",
            "attribute": str(POOL),
            "vertices": BREACH_VERTICES,
        },
    ]
    eta = breach_setup.apply_polygons(polygons=polygons, default=AMBIENT)

    np.testing.assert_allclose(eta, [POOL, POOL, AMBIENT, -DP0[ISLAND] - 0.01])

    h = z + eta
    assert h[HOLE] == pytest.approx(MAX_DEPTH + POOL)  # pooled at the gate
    assert h[RIM] == pytest.approx(-0.01)  # dry, and continuous with the island
    assert h[CHANNEL] == pytest.approx(DP0[CHANNEL] + AMBIENT)  # ambient preserved
    assert h[ISLAND] == pytest.approx(-0.01)  # island still dry


def test_pool_level_meets_island_formula_at_the_rim(breach_setup):
    """The constant equals -z-0.01 where the ground sits at min_depth."""
    z = _dredge(breach_setup)
    assert POOL == pytest.approx(-z[RIM] - 0.01)
