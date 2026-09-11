import shutil
from pathlib import Path

from click.testing import CliRunner

from schimpy.nc_bc_compatible import check_nc_bc_compatible, nc_bc_compatible_cli


def _write_bctides(path, elev_type=4, velocity_type=4):
    path.write_text(f"0 40\n0\n1\n2 {elev_type} {velocity_type} 3 3 ! ocean\n")


def _copy_mesh_inputs(tmp_path):
    data = Path(__file__).parent / "testdata" / "testmesh"
    reference = tmp_path / "reference"
    candidate = tmp_path / "candidate"
    reference.mkdir()
    candidate.mkdir()
    for directory in (reference, candidate):
        shutil.copy(data / "testmesh.gr3", directory / "hgrid.gr3")
        shutil.copy(data / "vgrid.in", directory / "vgrid.in")
        _write_bctides(directory / "bctides.in")
    return reference, candidate


def test_check_nc_bc_compatible_identical_meshes(tmp_path):
    reference, candidate = _copy_mesh_inputs(tmp_path)

    issues = check_nc_bc_compatible(
        "uv3D.th.nc",
        reference_hgrid=reference / "hgrid.gr3",
        candidate_hgrid=candidate / "hgrid.gr3",
        reference_vgrid=reference / "vgrid.in",
        candidate_vgrid=candidate / "vgrid.in",
        reference_bctides=reference / "bctides.in",
        candidate_bctides=candidate / "bctides.in",
        vgrid_version="5.8",
    )

    assert issues == []


def test_check_nc_bc_compatible_reports_level_mismatch(tmp_path):
    reference, candidate = _copy_mesh_inputs(tmp_path)
    candidate_vgrid = candidate / "vgrid.in"
    lines = candidate_vgrid.read_text().splitlines()
    tokens = lines[2].split()
    tokens[1] = str(int(tokens[1]) + 1)
    lines[2] = " ".join(tokens)
    candidate_vgrid.write_text("\n".join(lines) + "\n")

    issues = check_nc_bc_compatible(
        "elev2D.th.nc",
        reference_hgrid=reference / "hgrid.gr3",
        candidate_hgrid=candidate / "hgrid.gr3",
        reference_vgrid=reference / "vgrid.in",
        candidate_vgrid=candidate_vgrid,
        reference_bctides=reference / "bctides.in",
        candidate_bctides=candidate / "bctides.in",
        vgrid_version="5.8",
    )

    assert "levels reference=" in issues[0]


def test_nc_bc_compatible_cli_uses_directory_defaults(tmp_path):
    reference, candidate = _copy_mesh_inputs(tmp_path)

    result = CliRunner().invoke(
        nc_bc_compatible_cli,
        ["--file", "uv3D.th.nc", "--base-dir", str(reference), "--compare-dir", str(candidate), "--vgrid-version", "5.8", "--verbose"],
    )

    assert result.exit_code == 0, result.output
    assert "compatible" in result.output