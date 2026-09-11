"""Targeted repairs for SCHISM-generated hotstart files."""

import shutil
from pathlib import Path

import click
import netCDF4
import numpy as np

from schimpy.schism_mesh import read_mesh


def fix_hotstart_wet_dry(input_hotstart, output_hotstart, hgrid, h0=0.01):
    """Copy a hotstart and force geometrically impossible wet states dry.

    Existing dry node, side, and element states are preserved. Only wet nodes
    with ``dp + eta2 <= h0`` and the sides and elements incident to those nodes
    are changed.

    Parameters
    ----------
    input_hotstart : path-like
        SCHISM hotstart file to repair.
    output_hotstart : path-like
        New hotstart file to create. It must not already exist.
    hgrid : path-like
        Exact horizontal grid used by the run that wrote the hotstart.
    h0 : float, optional
        SCHISM wetting and drying threshold.

    Returns
    -------
    dict
        Counts of corrected nodes, sides, and elements.
    """
    input_path = Path(input_hotstart).resolve()
    output_path = Path(output_hotstart).resolve()
    hgrid_path = Path(hgrid).resolve()

    if input_path == output_path:
        raise ValueError("Input and output hotstart paths must differ")
    if output_path.exists():
        raise FileExistsError(f"Output hotstart already exists: {output_path}")

    mesh = read_mesh(str(hgrid_path))
    edges = np.asarray(mesh.edges)[:, :2]
    elements = np.asarray(mesh.elems)

    with netCDF4.Dataset(input_path) as source:
        expected_dimensions = {
            "node": mesh.n_nodes(),
            "side": mesh.n_edges(),
            "elem": mesh.n_elems(),
        }
        for name, expected in expected_dimensions.items():
            actual = len(source.dimensions[name]) if name in source.dimensions else None
            if actual != expected:
                raise ValueError(
                    f"Hotstart {name} dimension is {actual}; hgrid requires {expected}"
                )
        for name in ("eta2", "idry", "idry_s", "idry_e"):
            if name not in source.variables:
                raise ValueError(f"Hotstart is missing required variable: {name}")

    shutil.copy2(input_path, output_path)
    try:
        with netCDF4.Dataset(output_path, "r+") as output:
            eta = np.asarray(output.variables["eta2"][:]).reshape(-1)
            idry = np.asarray(output.variables["idry"][:]).reshape(-1)
            idry_s = np.asarray(output.variables["idry_s"][:]).reshape(-1)
            idry_e = np.asarray(output.variables["idry_e"][:]).reshape(-1)

            impossible_wet = (idry == 0) & (mesh.nodes[:, 2] + eta <= h0)
            idry[impossible_wet] = 1

            required_dry_s = idry[edges].max(axis=1)
            required_dry_e = np.array(
                [idry[element[element >= 0]].max() for element in elements]
            )
            corrected_sides = (idry_s == 0) & (required_dry_s == 1)
            corrected_elements = (idry_e == 0) & (required_dry_e == 1)
            idry_s[corrected_sides] = 1
            idry_e[corrected_elements] = 1

            output.variables["idry"][:] = idry
            output.variables["idry_s"][:] = idry_s
            output.variables["idry_e"][:] = idry_e
    except Exception:
        output_path.unlink(missing_ok=True)
        raise

    return {
        "nodes": int(impossible_wet.sum()),
        "sides": int(corrected_sides.sum()),
        "elements": int(corrected_elements.sum()),
    }


@click.command()
@click.argument(
    "input_hotstart", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.argument("output_hotstart", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--hgrid",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Exact hgrid.gr3 used by the run that wrote the hotstart.",
)
@click.option("--h0", default=0.01, show_default=True, type=click.FloatRange(min=0))
@click.help_option("-h", "--help")
def hotstart_fix_cli(input_hotstart, output_hotstart, hgrid, h0):
    """Copy a hotstart and repair impossible wet/dry flags."""
    counts = fix_hotstart_wet_dry(input_hotstart, output_hotstart, hgrid, h0)
    click.echo(
        f"Corrected {counts['nodes']} nodes, {counts['sides']} sides, and "
        f"{counts['elements']} elements; wrote {output_hotstart}"
    )


if __name__ == "__main__":
    hotstart_fix_cli()