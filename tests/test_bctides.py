# test_bctides.py

import pytest
import numpy as np

from schimpy.schism_mesh import BoundaryType
from schimpy.bctide import boundary

class DummyBoundary:
    def __init__(self, nodes, comment=None, btype=BoundaryType.OPEN):
        self.nodes = nodes
        self.comment = comment
        self.btype = btype

class DummyMesh:
    def __init__(self, boundaries, nodes=None):
        self.boundaries = boundaries
        if nodes is None:
            self.nodes = {i: (float(i), float(i)) for i in range(100)}
        else:
            self.nodes = nodes

def make_valid_bc_yaml(date="2011-01-28"):
    return {
        "name": "bctides.sed.in.3d",
        "modules": [
            {"name": "GEN", "num_tracers": 1},
            {"name": "SED", "num_tracers": 5},
            {"name": "AGE", "num_tracers": 2}
        ],
        "date": date,
        "mode": "baroclinic",
        "temperature_salinity_computation": True,
        "earth_tides": {
            "cutoff_depth": 40,
            "constituents": [
                {"K1": {
                    "amplitude": 0.0417807,
                    "angular_frequency": 0.00007292,
                    "node_factor": 1.0,
                    "earth_equilibrium_argument": 206.78674
                }},
                {"O1": {
                    "amplitude": 0.0751136,
                    "angular_frequency": 0.00006758,
                    "node_factor": 1.0,
                    "earth_equilibrium_argument": 0.0
                }},
                {"SSA": {
                    "amplitude": 0.0002282,
                    "angular_frequency": 0.00000040,
                    "node_factor": 1.0,
                    "earth_equilibrium_argument": 260.4229976
                }},
            ]
        },
        "boundary_forcing_tides": {
            "constituents": [
                {"Z0": {
                    "angular_frequency": 0.0,
                    "node_factor": 1.0,
                    "equilibrium_argument": 0.0
                }},
                {"O1": {
                    "angular_frequency": 0.00006758,
                    "node_factor": 1.11945,
                    "equilibrium_argument": 7.47302
                }},
                {"K1": {
                    "angular_frequency": 0.00007292,
                    "node_factor": 1.07391,
                    "equilibrium_argument": 206.78674
                }},
                {"M2": {
                    "angular_frequency": 0.00014023,
                    "node_factor": 0.97907,
                    "equilibrium_argument": 217.04138
                }},
            ]
        },
        "open_boundaries": [
            {
                "name": "ocean",
                "flow_direction": "bidirection",
                "variables": [
                    {"elevation": {
                        "source": "tidal",
                        "spatial": "local",
                        "constituents": [
                            {"Z0": [
                                {"amplitude": [1.0]*10},
                                {"phase": [0.0]*10}
                            ]},
                            {"O1": [
                                {"amplitude": [0.2]*10},
                                {"phase": [10.0]*10}
                            ]}
                        ]
                    }},
                    {"velocity": {
                        "source": "tidal",
                        "spatial": "local",
                        "constituents": [
                            {"Z0": [
                                {"u_amplitude": [0.1]*10},
                                {"u_phase": [0.0]*10},
                                {"v_amplitude": [0.1]*10},
                                {"v_phase": [0.0]*10}
                            ]}
                        ]
                    }},
                    {"temperature": {
                        "source": 10.0,
                        "nudge": 1.0
                    }},
                    {"salinity": {
                        "source": 5.0,
                        "nudge": 1.0
                    }},
                    {"tracers": [
                        {"module": "GEN", "source": [1.0], "nudge": 1.0},
                        {"module": "SED", "source": [0.1]*5, "nudge": 1.0},
                        {"module": "AGE", "source": [10.0, 20.0], "nudge": 1.0}
                    ]}
                ]
            }
        ]
    }

def make_dummy_mesh_and_bc_yaml(num_boundaries=1, nodes_per_boundary=10, date="2011-01-28"):
    boundaries = []
    nodes = {}
    for i in range(num_boundaries):
        node_ids = list(range(i * nodes_per_boundary, (i + 1) * nodes_per_boundary))
        boundaries.append(DummyBoundary(node_ids, comment=f'! Open boundary "{i}"'))
        for nid in node_ids:
            nodes[nid] = (float(nid), float(nid))
    mesh = DummyMesh(boundaries, nodes)
    bc_yaml = make_valid_bc_yaml(date)
    # Adjust open_boundaries to match num_boundaries
    bc_yaml["open_boundaries"] = [
        {
            "name": f"{i}",
            "flow_direction": "bidirection",
            "variables": [
                {"elevation": {
                    "source": "tidal",
                    "spatial": "local",
                    "constituents": [
                        {"Z0": [
                            {"amplitude": [1.0]*nodes_per_boundary},
                            {"phase": [0.0]*nodes_per_boundary}
                        ]}
                    ]
                }},
                {"velocity": {
                    "source": "tidal",
                    "spatial": "local",
                    "constituents": [
                        {"Z0": [
                            {"u_amplitude": [0.1]*nodes_per_boundary},
                            {"u_phase": [0.0]*nodes_per_boundary},
                            {"v_amplitude": [0.1]*nodes_per_boundary},
                            {"v_phase": [0.0]*nodes_per_boundary}
                        ]}
                    ]
                }},
                {"temperature": {
                    "source": 10.0,
                    "nudge": 1.0
                }},
                {"salinity": {
                    "source": 5.0,
                    "nudge": 1.0
                }},
                {"tracers": [
                    {"module": "GEN", "source": [1.0], "nudge": 1.0},
                    {"module": "SED", "source": [0.1]*5, "nudge": 1.0},
                    {"module": "AGE", "source": [10.0, 20.0], "nudge": 1.0}
                ]}
            ]
        }
        for i in range(num_boundaries)
    ]
    return mesh, bc_yaml

def test_write_bctides_valid(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1, nodes_per_boundary=10)
    bb = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_test.in"
    bb.write_bctides(str(out_file))
    content = out_file.read_text()
    assert ("2011-01-28" in content) or ("01/28/2011" in content)
    assert "K1" in content
    assert "Z0" in content
    ## module name is not explicitly appears in bctides.in
    #assert "GEN" in content
    #assert "SED" in content
    #assert "AGE" in content
    # The standard header comments mention "earth tidal potential", so only the
    # data fields can be checked for the leaked YAML keyword.
    data = "\n".join(line.split("!")[0] for line in content.splitlines())
    assert "tidal" not in data  # Should not literally write "tidal"
    assert "1 " in content  # number of open boundaries

def test_write_bctides_local_tidal_arrays(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1, nodes_per_boundary=10)
    bb = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_local_tidal.in"
    bb.write_bctides(str(out_file))
    content = out_file.read_text()
    # Check that all amplitudes and phases are present for each node
    assert content.count("1.0   0.0") == 10

def test_write_bctides_missing_fields(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml.pop("earth_tides", None)
    bc_yaml.pop("boundary_forcing_tides", None)
    bb = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_missing_fields.in"
    bb.write_bctides(str(out_file))
    content = out_file.read_text()
    assert "0 40" in content  # Default for missing earth_tides

def test_write_bctides_invalid_elev_source(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml["open_boundaries"][0]["variables"][0]["elevation"]["source"] = "not_supported"
    with pytest.raises(ValueError, match="elevation source .* is not supported"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_elev.in")

def test_write_bctides_invalid_velocity_source(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml["open_boundaries"][0]["variables"][1]["velocity"]["source"] = "not_supported"
    with pytest.raises(ValueError, match="velocity source .* is not supported"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_vel.in")

def test_write_bctides_invalid_tracer_module(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    tracer_mod = "NOT_A_TRACER"
    bc_yaml["open_boundaries"][0]["variables"][4]["tracers"].append(
        {"module": tracer_mod, "source": "th_text", "nudge": 0.5}
    )
    
    with pytest.raises(ValueError) as execinfo: 
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_tracer.in")
        assert "Tracer module '{tracer_mod}' is not defined in the bctides 'modules' section." in str(execinfo.value) 
        assert execinfo.type is ValueError

def test_write_bctides_tracer_constant_list(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bb = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_tracer_list.in"
    bb.write_bctides(str(out_file))
    content = out_file.read_text()
    assert "1.0" in content
    assert "0.1" in content
    assert "10.0" in content



# ---------------------------------------------------------------------------
# Source keywords, validation and writer output for single-boundary cases
# ---------------------------------------------------------------------------

import logging
import textwrap

from schimpy.bctide import (
    CONSTANT_CODE,
    ELEVATION_SOURCES,
    SCALAR_SOURCES,
    VELOCITY_SOURCES,
)
from schimpy.schism_yaml import load as schism_yaml_load

NNODE = 3
TIDE = [{"Z0": [{"amplitude": [1.0] * NNODE}, {"phase": [0.0] * NNODE}]}]
VTIDE = [{"Z0": [{"u_amplitude": [0.1] * NNODE}, {"u_phase": [0.0] * NNODE},
                 {"v_amplitude": [0.2] * NNODE}, {"v_phase": [5.0] * NNODE}]}]
TS = [{"temperature": {"source": "initial_profile", "nudge": 0.5}},
      {"salinity": {"source": "initial_profile", "nudge": 0.5}}]


def single_boundary(variables, flow_direction="bidirection", modules=None,
                    mode="baroclinic"):
    mesh = DummyMesh([DummyBoundary(list(range(NNODE)), comment='! Open boundary "b"')])
    bc_yaml = {
        "name": "bctides.in",
        "date": "2011-01-28",
        "mode": mode,
        "open_boundaries": [
            {"name": "b", "flow_direction": flow_direction, "variables": variables}
        ],
    }
    if mode == "barotropic":
        bc_yaml["temperature_salinity_computation"] = False
    if modules is not None:
        bc_yaml["modules"] = modules
    return mesh, bc_yaml


def write_single(tmp_path, variables, **kwargs):
    mesh, bc_yaml = single_boundary(variables, **kwargs)
    out = tmp_path / "bctides.in"
    boundary(mesh, bc_yaml).write_bctides(str(out))
    return out.read_text().splitlines()


def boundary_lines(lines):
    """Return the type-code line (as ints) and the parameter lines after it."""
    for k, line in enumerate(lines):
        if line.rstrip().endswith("!  b"):
            codes = [int(tok) for tok in line.split("!")[0].split()]
            return codes, lines[k + 1:]
    raise AssertionError("boundary line not found")


def elevation_spec(keyword):
    spec = {"source": keyword}
    if keyword in ("tidal", "tidal_th_nc"):
        spec["constituents"] = TIDE
    return spec


def velocity_spec(keyword):
    spec = {"source": keyword}
    if keyword in ("tidal", "tidal_th_nc", "tidal_th_nc_relaxed"):
        spec["constituents"] = VTIDE
    if keyword in ("th_nc_relaxed", "tidal_th_nc_relaxed"):
        spec.update(inflow_relax=0.8, outflow_relax=0.1)
    if keyword == "flather_1":
        spec.update(eta_m0=[0.1] * NNODE, vnormal_nodes=[[0.0, 0.01]] * NNODE)
    if keyword == "discharge_relation":
        spec["coefficients"] = [1.0, 2.0, 3.0, 4.0]
    return spec


@pytest.mark.parametrize("keyword,code", sorted(ELEVATION_SOURCES.items()))
def test_elevation_keyword_code(tmp_path, keyword, code):
    lines = write_single(tmp_path, [{"elevation": elevation_spec(keyword)}] + TS)
    codes, _ = boundary_lines(lines)
    assert codes[1] == code


@pytest.mark.parametrize("keyword,code", sorted(VELOCITY_SOURCES.items()))
def test_velocity_keyword_code(tmp_path, keyword, code):
    variables = [{"velocity": velocity_spec(keyword)}] + TS
    if keyword in ("th_nc_relaxed", "tidal_th_nc_relaxed"):
        variables.insert(0, {"elevation": {"source": "th_nc"}})
    codes, _ = boundary_lines(write_single(tmp_path, variables))
    assert codes[2] == code


@pytest.mark.parametrize("keyword,code", sorted(SCALAR_SOURCES.items()))
@pytest.mark.parametrize("variable,index", [("temperature", 3), ("salinity", 4)])
def test_scalar_keyword_code(tmp_path, keyword, code, variable, index):
    other = "salinity" if variable == "temperature" else "temperature"
    variables = [{variable: {"source": keyword, "nudge": 0.25}},
                 {other: {"source": "initial_profile", "nudge": 0.5}}]
    codes, params = boundary_lines(write_single(tmp_path, variables))
    assert codes[index] == code


@pytest.mark.parametrize("keyword,code", sorted(SCALAR_SOURCES.items()))
def test_tracer_keyword_code(tmp_path, keyword, code):
    variables = TS + [{"tracers": [{"module": "GEN", "source": keyword, "nudge": 0.3}]}]
    lines = write_single(tmp_path, variables, modules=[{"name": "GEN", "num_tracers": 1}])
    codes, params = boundary_lines(lines)
    assert codes[5] == code
    assert params[-1].strip() == "0.3"


@pytest.mark.parametrize("variable", ["elevation", "velocity", "temperature", "salinity"])
def test_numeric_source_is_constant(tmp_path, variable):
    spec = {"source": 0.5}
    if variable in ("temperature", "salinity"):
        spec["nudge"] = 1.0
        variables = [v for v in TS if variable not in v] + [{variable: spec}]
    else:
        variables = [{variable: spec}] + TS
    codes, params = boundary_lines(write_single(tmp_path, variables))
    index = {"elevation": 1, "velocity": 2, "temperature": 3, "salinity": 4}[variable]
    assert codes[index] == CONSTANT_CODE
    assert "0.5" in [p.strip() for p in params]


def test_constant_elevation_and_velocity_on_separate_lines(tmp_path):
    variables = [{"elevation": {"source": 0.25}}, {"velocity": {"source": -10.0}}] + TS
    codes, params = boundary_lines(write_single(tmp_path, variables))
    assert params[:4] == ["0.25", "-10.0", "0.5", "0.5"]


@pytest.mark.parametrize(
    "variable,keyword",
    [
        ("elevation", "2D elev history"),
        ("elevation", "2D_elev_history"),
        ("elevation", "elev_history"),
        ("elevation", "tidal_and_2D_elev_history"),
        ("elevation", "unspecified"),
        ("elevation", "constant"),
        ("velocity", "flux history"),
        ("velocity", "flux_history"),
        ("velocity", "3D_uv_history"),
        ("velocity", "flather 1"),
        ("temperature", "temperature history"),
        ("temperature", "temperature_history"),
        ("temperature", "initial profile"),
        ("salinity", "salinity_history"),
        ("salinity", "3D_salinity_history"),
    ],
)
def test_legacy_keywords_rejected(tmp_path, variable, keyword):
    spec = {"source": keyword, "nudge": 1.0}
    variables = [v for v in TS if variable not in v] + [{variable: spec}]
    with pytest.raises(ValueError, match=f"{variable} source .* is not supported"):
        write_single(tmp_path, variables)


@pytest.mark.parametrize("variable", ["temperature", "salinity", "tracers"])
def test_outflow_rejects_scalars_in_variables(tmp_path, variable):
    if variable == "tracers":
        spec = [{"module": "GEN", "source": "th_text", "nudge": 1.0}]
    else:
        spec = {"source": "th_text", "nudge": 1.0}
    variables = [{"velocity": {"source": "th_text"}}, {variable: spec}]
    with pytest.raises(ValueError, match="configured as an outflow"):
        write_single(tmp_path, variables, flow_direction="outflow",
                     modules=[{"name": "GEN", "num_tracers": 1}])


def test_outflow_with_velocity_only(tmp_path):
    lines = write_single(tmp_path, [{"velocity": {"source": "th_text"}}],
                         flow_direction="outflow")
    codes, _ = boundary_lines(lines)
    assert codes[1:5] == [0, 1, 0, 0]


def test_misindented_variable_rejected(tmp_path):
    variables = [{"velocity": None, "source": "th_text"}] + TS
    with pytest.raises(ValueError, match="indentation"):
        write_single(tmp_path, variables)


@pytest.mark.parametrize("keyword", ["flather_1", "discharge_relation"])
def test_velocity_requires_no_elevation(tmp_path, keyword):
    variables = [{"elevation": {"source": "th_nc"}},
                 {"velocity": velocity_spec(keyword)}] + TS
    with pytest.raises(ValueError, match="no elevation boundary"):
        write_single(tmp_path, variables)


@pytest.mark.parametrize("keyword", ["th_nc_relaxed", "tidal_th_nc_relaxed"])
def test_relaxed_velocity_requires_elevation(tmp_path, keyword):
    variables = [{"velocity": velocity_spec(keyword)}] + TS
    with pytest.raises(ValueError, match="requires an elevation boundary"):
        write_single(tmp_path, variables)


@pytest.mark.parametrize("key,value", [("inflow_relax", 1.5), ("outflow_relax", -0.1),
                                       ("inflow_relax", None)])
def test_relax_factor_range(tmp_path, key, value):
    vel = velocity_spec("th_nc_relaxed")
    if value is None:
        del vel[key]
    else:
        vel[key] = value
    variables = [{"elevation": {"source": "th_nc"}}, {"velocity": vel}] + TS
    with pytest.raises(ValueError, match=key):
        write_single(tmp_path, variables)


@pytest.mark.parametrize("coefs", [[1.0, 2.0, 3.0], None, [1, 2, 3, "x"]])
def test_discharge_relation_coefficients(tmp_path, coefs):
    vel = velocity_spec("discharge_relation")
    vel["coefficients"] = coefs
    with pytest.raises(ValueError, match="coefficients"):
        write_single(tmp_path, [{"velocity": vel}] + TS)


@pytest.mark.parametrize("nudge", [1.5, -0.2, "x", None])
def test_scalar_nudge_range(tmp_path, nudge):
    spec = {"source": "th_text"}
    if nudge is not None:
        spec["nudge"] = nudge
    variables = [{"temperature": spec}, TS[1]]
    with pytest.raises(ValueError, match="nudge"):
        write_single(tmp_path, variables)


def test_tracer_relax_key_rejected(tmp_path):
    variables = TS + [{"tracers": [{"module": "GEN", "source": "th_text", "relax": 1.0}]}]
    with pytest.raises(ValueError, match="named 'nudge'"):
        write_single(tmp_path, variables, modules=[{"name": "GEN", "num_tracers": 1}])


def test_flather_requires_per_node_values(tmp_path):
    vel = velocity_spec("flather_1")
    vel["eta_m0"] = [0.1]
    with pytest.raises(ValueError, match="eta_m0"):
        write_single(tmp_path, [{"velocity": vel}] + TS)


def test_write_relaxed_velocity(tmp_path):
    variables = [{"elevation": {"source": "th_nc"}},
                 {"velocity": velocity_spec("th_nc_relaxed")}] + TS
    codes, params = boundary_lines(write_single(tmp_path, variables))
    assert codes[1:3] == [4, -4]
    assert params[:3] == ["0.8 0.1", "0.5", "0.5"]


def test_write_tidal_relaxed_velocity(tmp_path):
    variables = [{"elevation": {"source": "th_nc"}},
                 {"velocity": velocity_spec("tidal_th_nc_relaxed")}] + TS
    codes, params = boundary_lines(write_single(tmp_path, variables))
    assert codes[1:3] == [4, -5]
    # SCHISM reads the relax factors before the tidal block
    assert params[0] == "0.8 0.1"
    assert params[1] == "Z0"
    assert params[2:2 + NNODE] == ["0.1 0.0 0.2 5.0"] * NNODE


def test_write_discharge_relation(tmp_path):
    variables = [{"velocity": velocity_spec("discharge_relation")}] + TS
    codes, params = boundary_lines(write_single(tmp_path, variables))
    assert codes[1:3] == [0, -2]
    assert params[0] == "1.0 2.0 3.0 4.0"


def test_write_flather(tmp_path):
    variables = [{"velocity": velocity_spec("flather_1")}] + TS
    codes, params = boundary_lines(write_single(tmp_path, variables))
    assert codes[1:3] == [0, -1]
    expected = (["eta_mean"] + ["0.1"] * NNODE + ["vn_mean"] + ["0.0 0.01"] * NNODE)
    assert params[:len(expected)] == expected


def test_tracer_subset_written_in_module_order(tmp_path):
    modules = [{"name": "SED", "num_tracers": 2}, {"name": "GEN", "num_tracers": 1}]
    variables = TS + [{"tracers": [
        {"module": "SED", "source": [0.1, 0.2], "nudge": 0.4},
    ]}]
    codes, params = boundary_lines(write_single(tmp_path, variables, modules=modules))
    # GEN precedes SED in SCHISM module order
    assert codes[5:] == [0, CONSTANT_CODE]
    assert params[2:4] == ["0.1 0.2", "0.4"]


def test_tracer_constant_length_checked(tmp_path):
    variables = TS + [{"tracers": [{"module": "SED", "source": [0.1], "nudge": 1.0}]}]
    with pytest.raises(ValueError, match="does not match number of tracers"):
        write_single(tmp_path, variables, modules=[{"name": "SED", "num_tracers": 2}])


def test_assignment_logged(tmp_path, caplog):
    variables = [{"elevation": {"source": "th_nc"}},
                 {"velocity": velocity_spec("th_nc_relaxed")}] + TS
    with caplog.at_level(logging.INFO, logger="schimpy.bctide"):
        write_single(tmp_path, variables)
    messages = [r.getMessage() for r in caplog.records]
    assert any(
        "Boundary 'b' (bidirection)" in m
        and "elevation=th_nc(4)" in m
        and "velocity=th_nc_relaxed(-4) inflow_relax=0.8 outflow_relax=0.1" in m
        and "temperature=initial_profile(3) nudge=0.5" in m
        for m in messages
    )


BAYDELTA_STYLE_YAML = textwrap.dedent(
    """
    name: bctides.in.3d
    date: 2011-01-28
    mode: baroclinic
    modules:
      - name: SED
        num_tracers: 3
    open_boundaries:
      - name: ocean
        flow_direction: bidirection
        variables:
          - elevation:
              source: th_nc
          - velocity:
              source: th_nc
          - temperature:
              source: initial_profile
              nudge: 0.01
          - salinity:
              source: initial_profile
              nudge: 0.01
          - tracers:
              - module: SED
                source: [0.007, 0.008, 0.0]
                nudge: 0.01
      - name: sjr
        flow_direction: inflow
        variables:
          - velocity:
              source: th_text
          - temperature:
              source: th_text
              nudge: 1.0
          - salinity:
              source: 0.074318
              nudge: 1.0
          - tracers:
              - module: SED
                source: th_text
                nudge: 1.0
      - name: swp
        flow_direction: outflow
        variables:
          - velocity:
              source: th_text
    """
)


def test_baydelta_style_yaml(tmp_path):
    yaml_file = tmp_path / "bctides.in.3d.yaml"
    yaml_file.write_text(BAYDELTA_STYLE_YAML)
    with open(yaml_file) as stream:
        bc_yaml = schism_yaml_load(stream)
    names = ["ocean", "sjr", "swp"]
    mesh = DummyMesh([
        DummyBoundary(list(range(k * 4, k * 4 + 4)), comment=f'! Open boundary "{n}"')
        for k, n in enumerate(names)
    ])
    out = tmp_path / "bctides.in.3d"
    boundary(mesh, bc_yaml).write_bctides(str(out))
    lines = out.read_text().splitlines()
    headers = [l for l in lines if l.rstrip().endswith(tuple(f"!  {n}" for n in names))]
    assert [h.split("!")[0].split() for h in headers] == [
        ["4", "4", "4", "3", "3", "2"],
        ["4", "0", "1", "1", "2", "1"],
        ["4", "0", "1", "0", "0", "0"],
    ]
