import netCDF4
import numpy as np
import pytest
from click.testing import CliRunner

from schimpy.hotstart_fix import fix_hotstart_wet_dry, hotstart_fix_cli
from schimpy.schism_mesh import read_mesh


def _write_hotstart(path, mesh):
    edges = np.asarray(mesh.edges)[:, :2]
    opposite_dry_node = np.flatnonzero(~np.any(edges == 1, axis=1))[0]
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("node", mesh.n_nodes())
        dataset.createDimension("side", mesh.n_edges())
        dataset.createDimension("elem", mesh.n_elems())
        dataset.createVariable("eta2", "f8", ("node",))[:] = [0.0, 2.99, 0.0]
        dataset.createVariable("idry", "i4", ("node",))[:] = [0, 0, 0]
        idry_s = dataset.createVariable("idry_s", "i4", ("side",))
        idry_s[:] = 0
        idry_s[opposite_dry_node] = 1
        dataset.createVariable("idry_e", "i4", ("elem",))[:] = [0]
        dataset.createVariable("tr_nd", "f8", ("node",))[:] = [10.0, 11.0, 12.0]
    return opposite_dry_node


def test_fix_hotstart_preserves_source_fields_and_dry_history(triangle_gr3, tmp_path):
    mesh = read_mesh(str(triangle_gr3))
    source = tmp_path / "source.nc"
    output = tmp_path / "fixed.nc"
    historical_dry_side = _write_hotstart(source, mesh)

    counts = fix_hotstart_wet_dry(source, output, triangle_gr3)

    assert counts["nodes"] == 1
    assert counts["elements"] == 1
    with netCDF4.Dataset(source) as original, netCDF4.Dataset(output) as fixed:
        np.testing.assert_array_equal(original["idry"][:], [0, 0, 0])
        np.testing.assert_array_equal(fixed["idry"][:], [0, 1, 0])
        np.testing.assert_array_equal(fixed["tr_nd"][:], original["tr_nd"][:])
        assert fixed["idry_s"][historical_dry_side] == 1
        assert fixed["idry_e"][0] == 1


def test_fix_hotstart_refuses_in_place_and_existing_output(triangle_gr3, tmp_path):
    source = tmp_path / "source.nc"
    _write_hotstart(source, read_mesh(str(triangle_gr3)))

    with pytest.raises(ValueError, match="must differ"):
        fix_hotstart_wet_dry(source, source, triangle_gr3)
    existing = tmp_path / "existing.nc"
    existing.touch()
    with pytest.raises(FileExistsError, match="already exists"):
        fix_hotstart_wet_dry(source, existing, triangle_gr3)


def test_hotstart_fix_cli(triangle_gr3, tmp_path):
    source = tmp_path / "source.nc"
    output = tmp_path / "fixed.nc"
    _write_hotstart(source, read_mesh(str(triangle_gr3)))

    result = CliRunner().invoke(
        hotstart_fix_cli,
        [str(source), str(output), "--hgrid", str(triangle_gr3)],
    )

    assert result.exit_code == 0, result.output
    assert "Corrected 1 nodes" in result.output
    assert output.exists()