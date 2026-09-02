collect_ignore = ["setup.py"]

import numpy as np
import pytest

from schimpy.schism_mesh import read_mesh
from schimpy.schism_vertical_mesh import SchismLocalVerticalMesh


# Three nodes is the minimum that forms an element, and is enough to cover the
# three regimes: submerged, dry land, and straddling the max(0.1, .) clamp.
TRIANGLE_DP = np.array([2.0, -3.0, 0.05])  # depth positive down, so bed elev = -dp
TRIANGLE_NVRT = 3


@pytest.fixture
def triangle_dp():
    return TRIANGLE_DP


@pytest.fixture
def triangle_gr3(tmp_path):
    path = tmp_path / "hgrid.gr3"
    xy = [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0)]
    lines = ["one triangle", "1 3 ! # of elements and nodes"]
    for i, ((x, y), dp) in enumerate(zip(xy, TRIANGLE_DP), start=1):
        lines.append(f"{i} {x:.8f} {y:.8f} {dp:.8f}")
    lines.append("1 3 1 2 3")
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def triangle_mesh(triangle_gr3):
    mesh = read_mesh(str(triangle_gr3))
    vmesh = SchismLocalVerticalMesh()
    # sigma runs -1 at the bottom to 0 at the surface, so column index 0 is the bed
    vmesh.init(np.tile(np.array([-1.0, -0.5, 0.0]), (3, 1)))
    vmesh.param["nvrt"] = TRIANGLE_NVRT
    mesh._vmesh = vmesh
    return mesh
