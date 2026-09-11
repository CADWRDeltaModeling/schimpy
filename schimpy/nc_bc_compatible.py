"""Check whether SCHISM NetCDF boundary forcing can be reused on another mesh."""

from pathlib import Path

import click
import numpy as np

from schimpy.logging_config import configure_logging, resolve_loglevel
from schimpy.schism_mesh import read_mesh
from schimpy.schism_vertical_mesh import read_vmesh


_FILE_SLOTS = {
    "elev2D.th.nc": 1,
    "uv3D.th.nc": 2,
}


def _uncommented_tokens(line):
    return line.split("!", 1)[0].split()


def _boundary_indices(bctides_path, forcing_file):
    """Return zero-based open-boundary indices configured for *forcing_file*."""
    try:
        slot = _FILE_SLOTS[Path(forcing_file).name]
    except KeyError as exc:
        supported = ", ".join(sorted(_FILE_SLOTS))
        raise ValueError(f"Unsupported forcing file {forcing_file!r}; use one of {supported}.") from exc

    with open(bctides_path) as bctides:
        _ = bctides.readline()
        _ = bctides.readline()
        tokens = _uncommented_tokens(bctides.readline())
        if not tokens:
            raise ValueError(f"Could not read the number of open boundaries from {bctides_path}.")
        n_open_boundaries = int(tokens[0])
        boundary_indices = []
        for boundary_index in range(n_open_boundaries):
            tokens = _uncommented_tokens(bctides.readline())
            if len(tokens) < 4:
                raise ValueError(
                    f"Could not read open boundary {boundary_index + 1} from {bctides_path}."
                )
            if int(tokens[slot]) == 4:
                boundary_indices.append(boundary_index)

    return boundary_indices


def check_nc_bc_compatible(
    forcing_file,
    *,
    reference_hgrid,
    candidate_hgrid,
    reference_vgrid,
    candidate_vgrid,
    reference_bctides,
    candidate_bctides,
    tolerance=0.1,
    vgrid_version="5.10",
):
    """Return incompatibilities between the boundary nodes of two SCHISM meshes.

    The file name selects the ``bctides.in`` boundary-condition slot: elevation
    for ``elev2D.th.nc`` and velocity for ``uv3D.th.nc``.  Coordinates are
    compared in open-boundary node order, and vertical level counts are
    calculated as ``nvrt - kbp``.
    """
    if tolerance < 0:
        raise ValueError("tolerance must be non-negative.")

    reference_indices = _boundary_indices(reference_bctides, forcing_file)
    candidate_indices = _boundary_indices(candidate_bctides, forcing_file)
    if reference_indices != candidate_indices:
        return [
            "Open boundaries configured for "
            f"{Path(forcing_file).name} differ: reference={reference_indices}, "
            f"candidate={candidate_indices}."
        ]

    reference_mesh = read_mesh(str(reference_hgrid))
    candidate_mesh = read_mesh(str(candidate_hgrid))
    reference_vmesh = read_vmesh(str(reference_vgrid), vgrid_version=vgrid_version)
    candidate_vmesh = read_vmesh(str(candidate_vgrid), vgrid_version=vgrid_version)
    issues = []

    for boundary_index in reference_indices:
        try:
            reference_nodes = reference_mesh.boundaries[boundary_index].nodes
            candidate_nodes = candidate_mesh.boundaries[boundary_index].nodes
        except IndexError:
            issues.append(f"Open boundary {boundary_index + 1} is missing from one mesh.")
            continue
        if len(reference_nodes) != len(candidate_nodes):
            issues.append(
                f"Open boundary {boundary_index + 1} node count differs: "
                f"reference={len(reference_nodes)}, candidate={len(candidate_nodes)}."
            )
            continue
        for node_ordinal, (reference_node, candidate_node) in enumerate(
            zip(reference_nodes, candidate_nodes), start=1
        ):
            reference_xy = reference_mesh.nodes[reference_node, :2]
            candidate_xy = candidate_mesh.nodes[candidate_node, :2]
            distance = float(np.linalg.norm(reference_xy - candidate_xy))
            reference_levels = reference_vmesh.n_vert_levels() - reference_vmesh.kbps[reference_node]
            candidate_levels = candidate_vmesh.n_vert_levels() - candidate_vmesh.kbps[candidate_node]
            if distance > tolerance or reference_levels != candidate_levels:
                issues.append(
                    f"Open boundary {boundary_index + 1}, node {node_ordinal}: "
                    f"distance={distance:.6g} (tolerance={tolerance:.6g}), "
                    f"levels reference={reference_levels}, candidate={candidate_levels}."
                )
    return issues


@click.command(help="Check whether a NetCDF boundary forcing file is reusable on a candidate SCHISM mesh.")
@click.option("--file", "forcing_file", required=True, type=click.Choice(sorted(_FILE_SLOTS)))
@click.option("--base-dir", "reference_dir", type=click.Path(file_okay=False, path_type=Path))
@click.option("--compare-dir", "candidate_dir", type=click.Path(file_okay=False, path_type=Path))
@click.option("--base-hgrid", "reference_hgrid", type=click.Path(dir_okay=False, path_type=Path))
@click.option("--compare-hgrid", "candidate_hgrid", type=click.Path(dir_okay=False, path_type=Path))
@click.option("--base-vgrid", "reference_vgrid", type=click.Path(dir_okay=False, path_type=Path))
@click.option("--compare-vgrid", "candidate_vgrid", type=click.Path(dir_okay=False, path_type=Path))
@click.option("--base-bctides", "reference_bctides", type=click.Path(dir_okay=False, path_type=Path))
@click.option("--compare-bctides", "candidate_bctides", type=click.Path(dir_okay=False, path_type=Path))
@click.option("--tol", "tolerance", default=0.1, show_default=True, type=float)
@click.option("--vgrid-version", default="5.10", show_default=True, type=click.Choice(["5.8", "5.10"]))
@click.option("--verbose", is_flag=True, help="Show compatible-boundary summary.")
@click.option("--raise", "raise_on_incompatible", is_flag=True, help="Raise an error instead of returning a nonzero exit status.")
@click.help_option("-h", "--help")
def nc_bc_compatible_cli(
    forcing_file,
    reference_dir,
    candidate_dir,
    reference_hgrid,
    candidate_hgrid,
    reference_vgrid,
    candidate_vgrid,
    reference_bctides,
    candidate_bctides,
    tolerance,
    vgrid_version,
    verbose,
    raise_on_incompatible,
):
    """Run the NetCDF boundary compatibility check."""
    level, console = resolve_loglevel(debug=verbose)
    configure_logging(package_name="schimpy", level=level, console=console)

    def resolve(path, directory, filename, label):
        if path is not None:
            return path
        if directory is None:
            raise click.UsageError(f"Specify --{label} or its corresponding directory option.")
        return directory / filename

    reference_hgrid = resolve(reference_hgrid, reference_dir, "hgrid.gr3", "base-hgrid")
    candidate_hgrid = resolve(candidate_hgrid, candidate_dir, "hgrid.gr3", "compare-hgrid")
    reference_vgrid = resolve(reference_vgrid, reference_dir, "vgrid.in", "base-vgrid")
    candidate_vgrid = resolve(candidate_vgrid, candidate_dir, "vgrid.in", "compare-vgrid")
    reference_bctides = resolve(reference_bctides, reference_dir, "bctides.in", "base-bctides")
    candidate_bctides = resolve(candidate_bctides, candidate_dir, "bctides.in", "compare-bctides")
    issues = check_nc_bc_compatible(
        forcing_file,
        reference_hgrid=reference_hgrid,
        candidate_hgrid=candidate_hgrid,
        reference_vgrid=reference_vgrid,
        candidate_vgrid=candidate_vgrid,
        reference_bctides=reference_bctides,
        candidate_bctides=candidate_bctides,
        tolerance=tolerance,
        vgrid_version=vgrid_version,
    )
    if issues:
        message = "\n".join(issues)
        if raise_on_incompatible:
            raise click.ClickException(message)
        click.echo(message, err=True)
        raise click.exceptions.Exit(1)
    if verbose:
        click.echo(f"{forcing_file} is compatible.")
