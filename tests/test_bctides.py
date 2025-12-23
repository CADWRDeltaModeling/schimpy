# test_bctide_write_bctides.py

import io
import tempfile
import datetime
import pytest
import numpy as np

from schimpy.schism_mesh import BoundaryType

# Import the boundary class directly from the module under test
from schimpy.bctide import boundary

class DummyBoundary:
    """A dummy mesh boundary for testing."""
    def __init__(self, nodes, comment=None, btype=BoundaryType.OPEN):
        self.nodes = nodes
        self.comment = comment
        self.btype = btype

class DummyMesh:
    """A dummy mesh object for testing."""
    def __init__(self, boundaries, nodes=None):
        self.boundaries = boundaries
        # nodes: dict {node_id: (x, y)}
        if nodes is None:
            self.nodes = {i: (float(i), float(i)) for i in range(100)}
        else:
            self.nodes = nodes

def make_valid_bc_yaml(date=None):
    """Return a minimal valid bc_yaml dict for testing."""
    if date is None:
        date = datetime.datetime(2020, 1, 1, 0, 0)
    earth_tidal_constituents = {
            "name": "M2",
            "angular_frequency": 1.9323,
            "node_factor": 1.0,
            "earth_equilibrium_argument": 0.0,
        }
    # If earth_tidals are present, put the number of constituents and cutoff depth
    # immediately after the "date" key (dict insertion order is preserved).
    return {
        "bctides": {
            "date": date,
            "earth_tidals": {
                "tidal_cutoff_depth": 40,
                "tidal_constituents": earth_tidal_constituents,
            },
            "boundary_forcing_tidals": {
            
                "tidal_constituents": [
                    {
                        "name": "M2",
                        "angular_frequency": 1.9323,
                        "node_factor": 1.0,
                        "earth_equilibrium_argument": 0.0,
                    }
                ],
            },
            "open_boundaries": [
                {
                    "name": "ocean",
                    "elevation_boundary": {
                        "source": "constant"
                    },
                    "velocity_boundary": {
                        "source": "constant"
                    },
                    "temperature_boundary": {
                        "source": "constant",
                        "nudge": 0.5
                    },
                    "salinity_boundary": {
                        "source": "constant",
                        "nudge": 0.5
                    }
                }
            ]
        }
    }

def make_dummy_mesh_and_bc_yaml(num_boundaries=1, nodes_per_boundary=3, date=None):
    """Create a dummy mesh and bc_yaml for testing."""
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
    bc_yaml["bctides"]["open_boundaries"] = [
        {
            "name": f"{i}",
            "elevation_boundary": {"source": "constant"},
            "velocity_boundary": {"source": "constant"},
            "temperature_boundary": {"source": "constant", "nudge": 0.5},
            "salinity_boundary": {"source": "constant", "nudge": 0.5}
        }
        for i in range(num_boundaries)
    ]
    return mesh, bc_yaml

def test_write_bctides_valid(tmp_path):
    """Test write_bctides with typical valid input."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=2, nodes_per_boundary=2)
    b = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_test.in"
    b.write_bctides(str(out_file))
    content = out_file.read_text()
    assert "2020-01-01 00:00" in content
    assert "M2" in content
    assert "constant" not in content  # Should not literally write "constant"
    # Should have number of open boundaries
    assert "2 " in content

def test_write_bctides_empty_open_boundaries(tmp_path):
    """Test with zero open boundaries (edge case)."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=0)
    b = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_empty.in"
    b.write_bctides(str(out_file))
    content = out_file.read_text()
    # Should still write date and header info
    assert "2020-01-01 00:00" in content

def test_write_bctides_mismatched_boundaries(tmp_path):
    """Test error when mesh and YAML open boundaries count mismatch."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=2)
    # Remove one open boundary from YAML
    bc_yaml["bctides"]["open_boundaries"] = bc_yaml["bctides"]["open_boundaries"][:1]
    with pytest.raises(ValueError, match="boundary YAML has different number of openbounary"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail.in")

def test_write_bctides_missing_fields(tmp_path):
    """Test missing optional fields (should not raise)."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    # Remove earth_tidals and bounary_forcing_tidals
    bc_yaml["bctides"].pop("earth_tidals", None)
    bc_yaml["bctides"].pop("bounary_forcing_tidals", None)
    b = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_missing_fields.in"
    b.write_bctides(str(out_file))
    content = out_file.read_text()
    assert "0 40" in content  # Default for missing earth_tidals

def test_write_bctides_invalid_tracer(tmp_path):
    """Test error for unsupported tracer module."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    # Add an invalid tracer
    bc_yaml["bctides"]["open_boundaries"][0]["tracers"] = [
        {"tracer": "NOT_A_TRACER", "source": "constant", "nudge": 0.5}
    ]
    with pytest.raises(ValueError, match="is not a supported tracer module"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_tracer.in")

def test_write_bctides_invalid_elev_source(tmp_path):
    """Test error for unsupported elevation boundary source."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml["bctides"]["open_boundaries"][0]["elevation_boundary"]["source"] = "not_supported"
    with pytest.raises(ValueError, match="elevation boundary is not supported"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_elev.in")

def test_write_bctides_invalid_velocity_source(tmp_path):
    """Test error for unsupported velocity boundary source."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml["bctides"]["open_boundaries"][0]["velocity_boundary"]["source"] = "not_supported"
    with pytest.raises(ValueError, match="velocity boundary is not supported"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_vel.in")

def test_write_bctides_invalid_salt_source(tmp_path):
    """Test error for unsupported salinity boundary source."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml["bctides"]["open_boundaries"][0]["salinity_boundary"]["source"] = "not_supported"
    with pytest.raises(ValueError, match="temperature boundary is not supported"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_salt.in")

def test_write_bctides_invalid_tracer_source_type(tmp_path):
    """Test error for unsupported tracer boundary source type."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    # Add a tracer with an unsupported source type (list of strings)
    bc_yaml["bctides"]["open_boundaries"][0]["tracers"] = [
        {"tracer": "GEN", "source": ["bad", "bad"], "nudge": 0.5}
    ]
    with pytest.raises(ValueError, match="boundary is not supported"):
        boundary(mesh, bc_yaml).write_bctides(tmp_path / "fail_tracer_type.in")

def test_write_bctides_tracer_constant_list(tmp_path):
    """Test tracer with constant list as source."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml["bctides"]["open_boundaries"][0]["tracers"] = [
        {"tracer": "GEN", "source": [1.0, 2.0, 3.0], "nudge": 0.5}
    ]
    b = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_tracer_list.in"
    b.write_bctides(str(out_file))
    content = out_file.read_text()
    assert "1.0 2.0 3.0" in content

def test_write_bctides_tracer_constant_scalar(tmp_path):
    """Test tracer with constant scalar as source."""
    mesh, bc_yaml = make_dummy_mesh_and_bc_yaml(num_boundaries=1)
    bc_yaml["bctides"]["open_boundaries"][0]["tracers"] = [
        {"tracer": "GEN", "source": 5.0, "nudge": 0.5}
    ]
    b = boundary(mesh, bc_yaml)
    out_file = tmp_path / "bctides_tracer_scalar.in"
    b.write_bctides(str(out_file))
    content = out_file.read_text()
    assert "5.0" in content

# To run these tests, use: pytest test_bctide_write_bctides.py
