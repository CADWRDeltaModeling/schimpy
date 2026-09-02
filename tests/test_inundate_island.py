"""Tests for the restoration island inundation generator."""

import numpy as np
import pytest
from shapely.geometry import Point, Polygon

from schimpy.schism_mesh import read_mesh
from schimpy.schism_setup import create_schism_setup
from schimpy import schism_yaml
from schimpy.inundate_island import (
    breach_axes,
    breach_polygon,
    build_inundation_artifacts,
    write_inundation_inputs,
    validate_structures,
    inundate_island_cli,
    local_edge_length,
)

FREEBOARD = 0.01
MIN_DEPTH = -0.20
MAX_DEPTH = 4.09


@pytest.fixture
def grid_mesh(tmp_path):
    """A 6x6 node grid of quads, 10 m spacing, all bed 1 m above datum."""
    path = tmp_path / "hgrid.gr3"
    n = 6
    nodes, elems = [], []
    for j in range(n):
        for i in range(n):
            nodes.append((i * 10.0, j * 10.0, -1.0))
    for j in range(n - 1):
        for i in range(n - 1):
            a = j * n + i + 1
            elems.append((a, a + 1, a + n + 1, a + n))
    lines = ["grid", "%d %d ! # of elements and nodes" % (len(elems), len(nodes))]
    for k, (x, y, dp) in enumerate(nodes, start=1):
        lines.append("%d %.8f %.8f %.8f" % (k, x, y, dp))
    for k, e in enumerate(elems, start=1):
        lines.append("%d 4 %d %d %d %d" % ((k,) + e))
    path.write_text("\n".join(lines) + "\n")
    return read_mesh(str(path))


def _spec(polygon=None):
    if polygon is None:
        polygon = {
            "name": "test_island",
            "type": "none",
            "vertices": [[0.0, 0.0], [50.0, 0.0], [50.0, 20.0], [0.0, 20.0]],
        }
    return [
        {
            "name": "test_island",
            "polygon": polygon,
            "breaches": [
                {
                    "name": "gate",
                    "left": [15.0, 25.0],
                    "right": [35.0, 25.0],
                    "min_depth": MIN_DEPTH,
                    "max_depth": MAX_DEPTH,
                }
            ],
        }
    ]


def test_axes_point_into_the_island():
    """Left to right looking into the island puts the interior on the +normal side."""
    _, along, into = breach_axes([0.0, 0.0], [10.0, 0.0])
    np.testing.assert_allclose(along, [1.0, 0.0])
    np.testing.assert_allclose(into, [0.0, 1.0])


def test_coincident_breach_points_rejected():
    with pytest.raises(ValueError, match="coincide"):
        breach_axes([5.0, 5.0], [5.0, 5.0])


def test_local_edge_length_matches_grid_spacing(grid_mesh):
    assert local_edge_length(grid_mesh, [25.0, 25.0], 15.0) == pytest.approx(10.0)


def test_polygon_encloses_the_dredge_ellipse():
    ring = breach_polygon([20.0, 25.0], [30.0, 25.0], 30.0)
    poly = Polygon(ring)
    # ellipse centre and both foci sit inside
    assert poly.contains(Point(25.0, 25.0))
    assert poly.contains(Point(20.0, 25.0))
    assert poly.contains(Point(30.0, 25.0))
    # and it reaches the full penetration depth through the levee
    assert poly.contains(Point(25.0, 25.0 + 14.0))


def test_pool_level_is_derived_from_min_depth(grid_mesh):
    art = build_inundation_artifacts(_spec(), grid_mesh, freeboard=FREEBOARD)
    breach = [p for p in art["elevation"]["polygons"] if p["name"].endswith("gate")][0]

    assert breach["type"] == "min"
    assert breach["attribute"] == pytest.approx(-MIN_DEPTH - FREEBOARD)


def test_breach_entry_follows_the_island_entry(grid_mesh):
    """Order is load bearing: the lower bound must be applied last."""
    names = [p["name"] for p in build_inundation_artifacts(_spec(), grid_mesh)["elevation"]["polygons"]]
    assert names.index("domain") < names.index("test_island")
    assert names.index("test_island") < names.index("test_island_gate")


def test_dredge_and_elevation_share_one_polygon(grid_mesh):
    art = build_inundation_artifacts(_spec(), grid_mesh)
    dredge = art["depth_enforcement"]["polygons"][0]
    breach = [p for p in art["elevation"]["polygons"] if p["name"] == dredge["name"]][0]
    assert dredge["vertices"] == breach["vertices"]


def test_structure_end_points_are_the_breach_points(grid_mesh):
    struct = build_inundation_artifacts(_spec(), grid_mesh)["structures"]["structures"][0]
    assert struct["end_points"] == [[15.0, 25.0], [35.0, 25.0]]
    assert struct["type"] == "transfer"
    assert struct["reference"] == "self"


def test_structures_file_carries_no_nudging(grid_mesh):
    """A nudging scalar here would be summed with the base file on include."""
    assert set(build_inundation_artifacts(_spec(), grid_mesh)["structures"]) == {"structures"}


def test_major_axis_defaults_to_twice_the_penetration(grid_mesh):
    """The factor is the reach each way, so the full axis is twice it."""
    art = build_inundation_artifacts(_spec(), grid_mesh)
    assert "major_axis_len=60.0" in art["depth_enforcement"]["polygons"][0]["attribute"]


def test_max_depth_must_be_deeper_than_min_depth(grid_mesh):
    spec = _spec()
    spec[0]["breaches"][0]["max_depth"] = MIN_DEPTH - 1.0
    with pytest.raises(ValueError, match="no deeper than"):
        build_inundation_artifacts(spec, grid_mesh)


@pytest.fixture
def gated_mesh(tmp_path):
    """Two equal blocks of quads joined only by one gate element in a levee row.

    Rows of elements span y 0-30 and y 40-70. The row between them is absent
    except at x 30-40, so either block is reachable from the other only through
    that one element. The blocks are the same size so neither trips the fill
    guard on its own.
    """
    path = tmp_path / "gated.gr3"
    n = 8
    nodes = [(i * 10.0, j * 10.0, -1.0) for j in range(n) for i in range(n)]
    elems = []
    for j in range(n - 1):
        for i in range(n - 1):
            if j == 3 and i not in (2, 3, 4):
                continue
            a = j * n + i + 1
            elems.append((a, a + 1, a + n + 1, a + n))
    lines = ["gated", "%d %d ! # of elements and nodes" % (len(elems), len(nodes))]
    for k, (x, y, dp) in enumerate(nodes, start=1):
        lines.append("%d %.8f %.8f %.8f" % (k, x, y, dp))
    for k, e in enumerate(elems, start=1):
        lines.append("%d 4 %d %d %d %d" % ((k,) + e))
    path.write_text("\n".join(lines) + "\n")
    return read_mesh(str(path))


def _gated_spec():
    return [
        {
            "name": "gated_island",
            "polygon": "infer",
            "breaches": [
                {
                    "name": "gate",
                    "left": [20.0, 35.0],
                    "right": [50.0, 35.0],
                    "min_depth": MIN_DEPTH,
                    "max_depth": MAX_DEPTH,
                }
            ],
        }
    ]


def test_infer_recovers_the_enclosed_side(gated_mesh):
    """Left to right looking into the island selects the upper block."""
    art = build_inundation_artifacts(_gated_spec(), gated_mesh)
    island = [p for p in art["elevation"]["polygons"] if p["name"] == "gated_island"][0]
    poly = Polygon(island["vertices"])

    assert poly.contains(Point(25.0, 55.0))  # upper block, the island side
    assert not poly.contains(Point(25.0, 5.0))  # lower block, across the gate


def test_infer_follows_the_left_right_ordering(gated_mesh):
    """Swapping the two points selects the other side of the levee."""
    spec = _gated_spec()
    spec[0]["breaches"][0]["left"] = [40.0, 35.0]
    spec[0]["breaches"][0]["right"] = [30.0, 35.0]

    art = build_inundation_artifacts(spec, gated_mesh)
    island = [p for p in art["elevation"]["polygons"] if p["name"] == "gated_island"][0]
    poly = Polygon(island["vertices"])

    assert poly.contains(Point(25.0, 5.0))
    assert not poly.contains(Point(25.0, 55.0))


def test_inferred_outline_excludes_the_gate_element(gated_mesh):
    """The gate is a barrier, so the breach polygon covers it instead."""
    art = build_inundation_artifacts(_gated_spec(), gated_mesh)
    island = [p for p in art["elevation"]["polygons"] if p["name"] == "gated_island"][0]
    assert not Polygon(island["vertices"]).contains(Point(35.0, 35.0))

@pytest.fixture
def two_gate_mesh(tmp_path):
    """A block penned in by two openings, neither of which closes it alone.

    Element rows 3 and 8 are absent except at x 20-50, so the middle band is
    reachable only through two gates, one at each end. Blocking either one on
    its own still leaves the other open.
    """
    path = tmp_path / "two_gate.gr3"
    n = 13
    nodes = [(i * 10.0, j * 10.0, -1.0) for j in range(n) for i in range(n)]
    elems = []
    for j in range(n - 1):
        for i in range(n - 1):
            if j in (3, 8) and i not in (2, 3, 4):
                continue
            a = j * n + i + 1
            elems.append((a, a + 1, a + n + 1, a + n))
    lines = ["two gate", "%d %d ! # of elements and nodes" % (len(elems), len(nodes))]
    for k, (x, y, dp) in enumerate(nodes, start=1):
        lines.append("%d %.8f %.8f %.8f" % (k, x, y, dp))
    for k, e in enumerate(elems, start=1):
        lines.append("%d 4 %d %d %d %d" % ((k,) + e))
    path.write_text("\n".join(lines) + "\n")
    return read_mesh(str(path))


def _two_gate_spec(flip_second=False):
    """Both breaches look into the middle band: south gate north, north gate south."""
    south = {"left": [20.0, 35.0], "right": [50.0, 35.0]}
    north = {"left": [50.0, 85.0], "right": [20.0, 85.0]}
    if flip_second:
        north = {"left": [20.0, 85.0], "right": [50.0, 85.0]}
    breaches = []
    for label, pts in (("south", south), ("north", north)):
        breaches.append(
            dict(
                name=label,
                min_depth=MIN_DEPTH,
                max_depth=MAX_DEPTH,
                major_axis_len=40.0,
                **pts,
            )
        )
    return [{"name": "penned", "polygon": "infer", "breaches": breaches}]


def test_two_breaches_together_pen_in_one_island(two_gate_mesh):
    """Neither gate closes the band alone; together they do."""
    art = build_inundation_artifacts(_two_gate_spec(), two_gate_mesh)
    island = [p for p in art["elevation"]["polygons"] if p["name"] == "penned"][0]
    poly = Polygon(island["vertices"])

    assert poly.contains(Point(65.0, 60.0))  # the penned band
    assert not poly.contains(Point(65.0, 15.0))  # south of the first gate
    assert not poly.contains(Point(65.0, 105.0))  # north of the second


def test_disagreeing_breaches_are_rejected(two_gate_mesh):
    """A reversed pair points away from the island and must be caught."""
    with pytest.raises(ValueError, match="not in the same region"):
        build_inundation_artifacts(_two_gate_spec(flip_second=True), two_gate_mesh)


def test_validation_accepts_a_breach_that_spans_the_gate(gated_mesh, tmp_path):
    """Node pairs and the reference pair must land inside the dredge polygon."""
    art = build_inundation_artifacts(_gated_spec(), gated_mesh)
    checked = validate_structures(art, str(tmp_path / "gated.gr3"))
    assert checked["gated_island_gate"] > 0


def test_validation_rejects_a_dredge_too_narrow_for_its_structure(gated_mesh, tmp_path):
    """A short penetration leaves the reference nodes on undredged ground."""
    spec = _gated_spec()
    spec[0]["breaches"][0]["major_axis_len"] = 6.0
    spec[0]["breaches"][0]["dredge_pad"] = 0.0
    art = build_inundation_artifacts(spec, gated_mesh)
    with pytest.raises(ValueError, match="outside its dredge polygon"):
        validate_structures(art, str(tmp_path / "gated.gr3"))


def test_cli_writes_all_four_artifacts(gated_mesh, tmp_path):
    from click.testing import CliRunner

    spec_path = tmp_path / "breaches.yaml"
    schism_yaml.safe_dump({"islands": _gated_spec()}, open(spec_path, "w"))

    result = CliRunner().invoke(
        inundate_island_cli,
        [
            "--spec", str(spec_path),
            "--hgrid", str(tmp_path / "gated.gr3"),
            "--out-dir", str(tmp_path / "cli_out"),
        ],
    )
    assert result.exit_code == 0, result.output
    written = sorted(p.name for p in (tmp_path / "cli_out").iterdir())
    assert written == [
        "depth_enforce_inundate.yaml",
        "elev_inundate.yaml",
        "hydraulic_structures_inundate.yaml",
        "inundate_regions.yaml",
    ]


def test_cli_reports_a_bad_spec_without_a_traceback(gated_mesh, tmp_path):
    from click.testing import CliRunner

    spec_path = tmp_path / "empty.yaml"
    spec_path.write_text("islands: []\n")

    result = CliRunner().invoke(
        inundate_island_cli,
        ["--spec", str(spec_path), "--hgrid", str(tmp_path / "gated.gr3")],
    )
    assert result.exit_code != 0
    assert "no 'islands' list" in result.output


def test_island_without_polygon_rejected(grid_mesh):
    spec = _spec()
    del spec[0]["polygon"]
    with pytest.raises(ValueError, match="needs a polygon"):
        build_inundation_artifacts(spec, grid_mesh)


def test_bare_vertex_list_rejected(grid_mesh):
    """The polygon must use the same form as the preprocessor polygon files."""
    with pytest.raises(ValueError, match="mapping with 'vertices'"):
        build_inundation_artifacts(_spec(polygon=[[0.0, 0.0], [1.0, 0.0]]), grid_mesh)


def test_unknown_polygon_keyword_rejected(grid_mesh):
    with pytest.raises(ValueError, match="or 'infer'"):
        build_inundation_artifacts(_spec(polygon="guess"), grid_mesh)


def test_infer_refuses_when_the_fill_escapes(grid_mesh):
    """A breach that encloses nothing fills the whole grid rather than an island."""
    with pytest.raises(ValueError, match="do not pen in the island"):
        build_inundation_artifacts(_spec(polygon="infer"), grid_mesh)


def test_written_yaml_has_no_anchors(grid_mesh, tmp_path):
    """Closing a ring with the same object would serialize as a yaml alias."""
    art = build_inundation_artifacts(_spec(), grid_mesh)
    paths = write_inundation_inputs(art, str(tmp_path / "out"))
    for path in paths.values():
        text = open(path).read()
        assert "&id" not in text and "*id" not in text


def test_written_files_reload_and_dredge_applies(grid_mesh, tmp_path):
    art = build_inundation_artifacts(_spec(), grid_mesh)
    paths = write_inundation_inputs(art, str(tmp_path / "out"))

    from schimpy import schism_yaml

    with open(paths["depth_enforcement"]) as fh:
        reloaded = schism_yaml.load(fh)

    setup = create_schism_setup(str(tmp_path / "hgrid.gr3"))
    z = setup.apply_polygons(
        polygons=reloaded["polygons"],
        default=None,
        global_imports=["schimpy.ellipse.ellipse"],
    )
    # the opening is dredged and nothing in the footprint is left above the rim
    assert z.max() > 1.0
    ring = Polygon(art["depth_enforcement"]["polygons"][0]["vertices"])
    inside = [
        i for i in range(grid_mesh.n_nodes()) if ring.contains(Point(grid_mesh.nodes[i][:2]))
    ]
    assert inside
    assert np.all(z[inside] >= MIN_DEPTH - 1e-9)
