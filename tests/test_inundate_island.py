"""Tests for the restoration island inundation generator."""

import numpy as np
import pytest
from shapely.geometry import Point, Polygon

from schimpy.schism_mesh import read_mesh
from schimpy.schism_setup import create_schism_setup
from schimpy import schism_yaml
from schimpy.inundate_island import (
    build_breach_schedule,
    write_breach_timeseries,
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


def test_hotstart_regions_preserve_distinct_island_names(grid_mesh):
    regions = build_inundation_artifacts(_spec(), grid_mesh)["regions"]["polygons"]
    assert [region["name"] for region in regions] == ["domain", "test_island"]


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
    with pytest.raises(ValueError, match="fills a different region|went around"):
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

    config_path = tmp_path / "breaches.yaml"
    schism_yaml.safe_dump({"islands": _gated_spec()}, open(config_path, "w"))

    result = CliRunner().invoke(
        inundate_island_cli,
        [
            "--config", str(config_path),
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
    config_path = tmp_path / "empty.yaml"
    config_path.write_text("islands: []\n")

    result = CliRunner().invoke(
        inundate_island_cli,
        ["--config", str(config_path), "--hgrid", str(tmp_path / "gated.gr3")],
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
    """A breach that encloses nothing lets the fill reach its own far side."""
    with pytest.raises(ValueError, match="went around the barrier"):
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


def _weir_schedule(**overrides):
    schedule = {
        "name": "breach",
        "type": "weir",
        "configuration": {
            "n_duplicates": 1,
            "elevation": 1.5,
            "width": 20.0,
            "coefficient": 0.6,
            "op_downstream": 1.0e-5,
            "op_upstream": 0.0,
        },
        "breach_date": "2022-04-15",
    }
    schedule.update(overrides)
    return schedule


def test_schedule_holds_then_ramps_then_deinstalls():
    frame = build_breach_schedule(_weir_schedule())
    assert len(frame) == 14
    assert str(frame.index[0].date()) == "2005-01-01"
    # the ramp ends one increment before the breach, which then deinstalls
    assert frame["install"].tolist() == [1] * 13 + [0]
    assert frame["op_downstream"].iloc[0] == pytest.approx(1.0e-5)
    assert frame["op_downstream"].iloc[-2] == pytest.approx(0.6)
    assert frame["op_downstream"].is_monotonic_increasing
    assert str(frame.index[-1]) == "2022-04-15 00:00:00"
    assert str(frame.index[1]) == "2022-04-14 00:00:00"


def test_schedule_columns_match_the_fortran_read_order():
    """read_struct_ts reads ttt, install, nduplicate, op_down, op_up, elev, width."""
    frame = build_breach_schedule(_weir_schedule())
    assert list(frame.columns) == [
        "install",
        "n_duplicates",
        "op_downstream",
        "op_upstream",
        "elevation",
        "width",
    ]


def test_transfer_schedule_ramps_the_flow():
    schedule = _weir_schedule(type="transfer", configuration={"flow": 0.0})
    frame = build_breach_schedule(schedule, ramp_target=2.0)
    assert list(frame.columns) == ["install", "flow"]
    assert frame["flow"].iloc[0] == pytest.approx(0.0)
    assert frame["flow"].iloc[-2] == pytest.approx(2.0)


def test_epoch_must_precede_the_ramp():
    with pytest.raises(ValueError, match="must precede the ramp"):
        build_breach_schedule(_weir_schedule(), epoch_start="2022-04-14T12:00")


def test_elapsed_th_has_no_header_and_dated_does(tmp_path):
    """SCHISM reads the .th with a list-directed read, so it cannot carry a header."""
    written = write_breach_timeseries(
        [_weir_schedule()], str(tmp_path), run_start="2020-09-30"
    )
    dated, elapsed = written["breach"]

    dated_lines = open(dated).read().splitlines()
    assert dated_lines[0].split() == [
        "datetime", "install", "ndup", "op_down", "op_up", "elev", "width",
    ]
    assert dated_lines[1].startswith("2005-01-01T00:00")

    elapsed_lines = open(elapsed).read().splitlines()
    assert len(elapsed_lines) == len(dated_lines) - 1
    for line in elapsed_lines:
        float(line.split()[0])
    # 2005-01-01 sits before the run origin, and the breach lands after it
    assert float(elapsed_lines[0].split()[0]) < 0.0
    assert float(elapsed_lines[-1].split()[0]) == pytest.approx(48556800.0)


def test_scheduled_structures_are_flagged_for_a_time_series(gated_mesh):
    spec = _gated_spec()
    spec[0]["breaches"][0]["structure"] = {
        "type": "weir",
        "configuration": {
            "n_duplicates": 1,
            "elevation": 1.5,
            "width": 20.0,
            "coefficient": 0.6,
            "op_downstream": 1.0e-5,
            "op_upstream": 0.0,
        },
    }
    art = build_inundation_artifacts(spec, gated_mesh, breach_date="2022-04-15")
    entry = art["structures"]["structures"][0]
    assert entry["type"] == "weir"
    assert entry["configuration"]["use_time_series"] == 1
    assert len(art["schedules"]) == 1


def test_no_breach_date_leaves_the_transfer_default_alone(gated_mesh):
    art = build_inundation_artifacts(_gated_spec(), gated_mesh)
    entry = art["structures"]["structures"][0]
    assert entry["type"] == "transfer"
    assert entry["configuration"] == {"flow": 0.0}
    assert art["schedules"] == []


def test_breach_date_narrowest_setting_wins(two_gate_mesh):
    """Breach beats island beats file, so sites in one run can differ."""
    spec = _two_gate_spec()
    spec[0]["breach_date"] = "2023-01-01"
    spec[0]["breaches"][0]["breach_date"] = "2022-04-15"

    art = build_inundation_artifacts(spec, two_gate_mesh, breach_date="2030-12-31")
    dates = {s["name"]: s["breach_date"] for s in art["schedules"]}
    assert len(dates) == 2
    assert sorted(dates.values()) == ["2022-04-15", "2023-01-01"]


def test_sites_years_apart_get_independent_schedules(two_gate_mesh, tmp_path):
    spec = _two_gate_spec()
    spec[0]["breaches"][0]["breach_date"] = "2022-04-15"
    spec[0]["breaches"][1]["breach_date"] = "2024-09-01"

    art = build_inundation_artifacts(spec, two_gate_mesh)
    written = write_breach_timeseries(
        art["schedules"], str(tmp_path / "multi"), run_start="2020-09-30"
    )
    assert len(written) == 2
    last_rows = {}
    for name, (_dated, elapsed) in written.items():
        last_rows[name] = float(open(elapsed).read().splitlines()[-1].split()[0])
    values = sorted(last_rows.values())
    # roughly two and a half years apart, and both after the origin
    assert values[0] > 0.0
    assert values[1] - values[0] == pytest.approx(870 * 86400, rel=0.05)


def test_a_breach_without_a_date_gets_no_schedule(two_gate_mesh):
    spec = _two_gate_spec()
    spec[0]["breaches"][0]["breach_date"] = "2022-04-15"
    art = build_inundation_artifacts(spec, two_gate_mesh)
    assert [s["name"] for s in art["schedules"]] == [
        "%s_%s" % (spec[0]["name"], spec[0]["breaches"][0]["name"])
    ]
    flagged = [
        s["name"]
        for s in art["structures"]["structures"]
        if "use_time_series" in s["configuration"]
    ]
    assert len(flagged) == 1


def test_cli_writes_th_files(gated_mesh, tmp_path):
    from click.testing import CliRunner

    spec = _gated_spec()
    spec[0]["breaches"][0]["structure"] = {
        "type": "weir",
        "configuration": {
            "n_duplicates": 1,
            "elevation": 1.5,
            "width": 20.0,
            "coefficient": 0.6,
            "op_downstream": 1.0e-5,
            "op_upstream": 0.0,
        },
    }
    config_path = tmp_path / "breaches_th.yaml"
    schism_yaml.safe_dump(
        {"islands": spec, "breach_date": "2022-04-15", "run_start": "2020-09-30"},
        open(config_path, "w"),
    )

    out_dir = tmp_path / "th_out"
    result = CliRunner().invoke(
        inundate_island_cli,
        [
            "--config", str(config_path),
            "--hgrid", str(tmp_path / "gated.gr3"),
            "--out-dir", str(out_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    name = spec[0]["breaches"][0]["name"]
    struct_name = "%s_%s" % (spec[0]["name"], name)
    assert (out_dir / "th_files" / "dated" / ("%s.th" % struct_name)).exists()
    assert (out_dir / "th_files" / "elapsed" / ("%s.th" % struct_name)).exists()


def test_breach_date_without_run_start_is_rejected(gated_mesh, tmp_path):
    from click.testing import CliRunner

    config_path = tmp_path / "no_start.yaml"
    schism_yaml.safe_dump(
        {"islands": _gated_spec(), "breach_date": "2022-04-15"},
        open(config_path, "w"),
    )
    result = CliRunner().invoke(
        inundate_island_cli,
        [
            "--config", str(config_path),
            "--hgrid", str(tmp_path / "gated.gr3"),
            "--out-dir", str(tmp_path / "no_start_out"),
        ],
    )
    assert result.exit_code != 0
    assert "run_start is needed" in result.output
