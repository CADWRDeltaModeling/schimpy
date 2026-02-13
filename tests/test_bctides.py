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
                        {"module": "GEN", "source": [1.0], "relax": 1.0},
                        {"module": "SED", "source": [0.1]*5, "relax": 1.0},
                        {"module": "AGE", "source": [10.0, 20.0], "relax": 1.0}
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
                    {"module": "GEN", "source": [1.0], "relax": 1.0},
                    {"module": "SED", "source": [0.1]*5, "relax": 1.0},
                    {"module": "AGE", "source": [10.0, 20.0], "relax": 1.0}
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
    assert "tidal" not in content  # Should not literally write "tidal"
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
    with pytest.raises(ValueError, match="elevation boundary is not supported"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_elev.in")

def test_write_bctides_invalid_velocity_source(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml["open_boundaries"][0]["variables"][1]["velocity"]["source"] = "not_supported"
    with pytest.raises(ValueError, match="velocity boundary is not supported"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_vel.in")

def test_write_bctides_invalid_tracer_module(tmp_path):
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    tracer_mod = "NOT_A_TRACER"
    bc_yaml["open_boundaries"][0]["variables"][4]["tracers"].append(
        {"module": tracer_mod, "source": "constant", "relax": 0.5}
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

# To run these tests, use: pytest test_bctides.py
