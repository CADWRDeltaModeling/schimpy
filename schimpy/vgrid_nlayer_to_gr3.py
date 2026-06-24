#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create GR3 files with per-node vertical level and layer counts from a SCHISM vgrid."""

import numpy as np
import click

from schimpy.schism_mesh import read_mesh, write_mesh


def _compute_nlevel_nlayer(mesh):
    if mesh.vmesh is None:
        raise ValueError("Mesh is missing a vertical mesh. Provide both hgrid and vgrid.")

    vmesh = mesh.vmesh
    if not hasattr(vmesh, "kbps") or vmesh.kbps is None:
        raise ValueError("Vertical mesh is missing kbps values required for level count computation.")

    nlevels = np.asarray(vmesh.n_vert_levels() - np.asarray(vmesh.kbps), dtype=float)
    nlayers = np.maximum(nlevels - 1.0, 0.0)
    return nlevels, nlayers


def vgrid_nlayer_to_gr3(
    hgrid,
    vgrid,
    output_nlevel=None,
    output_nlayer=None,
    vgrid_version="5.10",
):
    """Write GR3 files for per-node vertical levels and/or layers from a vgrid."""
    if output_nlevel is None and output_nlayer is None:
        raise ValueError("At least one of --nlevel or --nlayer must be provided.")
    if output_nlevel is not None and output_nlayer is not None and output_nlevel == output_nlayer:
        raise ValueError("--nlevel and --nlayer must write to different output files.")

    mesh = read_mesh(hgrid, vgrid, vgrid_version=vgrid_version)
    nlevels, nlayers = _compute_nlevel_nlayer(mesh)

    if output_nlevel is not None:
        write_mesh(mesh, output_nlevel, node_attr=nlevels)
    if output_nlayer is not None:
        write_mesh(mesh, output_nlayer, node_attr=nlayers)

    return nlevels, nlayers


@click.command()
@click.option(
    "--hgrid",
    default="hgrid.gr3",
    help="Horizontal grid file name.",
)
@click.option(
    "--vgrid",
    default="vgrid.in",
    help="Vertical grid file name.",
)
@click.option(
    "--vgrid_version",
    default="5.10",
    help="SCHISM vgrid format version ('5.8' or '5.10').",
)
@click.option(
    "--nlevel",
    default=None,
    help="Output GR3 file for per-node vertical level counts.",
)
@click.option(
    "--nlayer",
    default=None,
    help="Output GR3 file for per-node per-node layer counts.",
)
def vgrid_nlayer_to_gr3_cli(hgrid, vgrid, vgrid_version, nlevel, nlayer):
    """CLI utility to write nlevel/nlayer GR3 files from a SCHISM vgrid."""
    vgrid_nlayer_to_gr3(
        hgrid=hgrid,
        vgrid=vgrid,
        output_nlevel=nlevel,
        output_nlayer=nlayer,
        vgrid_version=vgrid_version,
    )


if __name__ == "__main__":
    vgrid_nlayer_to_gr3_cli()
