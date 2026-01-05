#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create a mesh file with the number of vertical levels as nodal value"""
from schimpy.schism_mesh import read_mesh, write_mesh
import numpy as np
import click

def create_mesh_n_levels(hgrid, vgrid, output):
    """
    Create a mesh file with the number of vertical levels as nodal values.
    This function reads a horizontal grid (hgrid) and vertical grid (vgrid),
    calculates the number of active vertical levels at each node, and writes
    the result to an output file with the level counts as nodal attributes.
    Parameters
    ----------
    hgrid : str or Path
        Path to the horizontal grid file.
    vgrid : str or Path
        Path to the vertical grid file.
    output : str or Path
        Path to the output mesh file where the mesh with vertical level
        counts will be written.
    Returns
    -------
    None
        The function writes the output to a file and returns nothing.
    Notes
    -----
    The number of vertical levels at each node is calculated as:
        n_levels = total_vertical_levels - kbp
    where kbp is the bottom index for each node.
    Examples
    --------
    >>> create_mesh_levels('hgrid.gr3', 'vgrid.in', 'mesh_levels.gr3')
    """
    
    mesh = read_mesh(hgrid, vgrid)
    print(mesh.vmesh.n_vert_levels())
    n_levels = mesh.vmesh.n_vert_levels() - np.array(mesh.vmesh.kbps)
    write_mesh(mesh, output, node_attr=n_levels)

@click.command()
@click.option(
    "--hgrid",
    default="hgrid.gr3",
    help="Horizontal grid name.",
)
@click.option(
    "--vgrid",
    default="vgrid.in",
    help="Vertical grid file name.",
)
@click.option(
    "--output",
    default="n_levels.gr3",
    help="Output file name.",
)
def create_mesh_n_levels_cli(hgrid, vgrid, output):
    """Command line utility for creting a mesh file with the number of vertical levels"""
    create_mesh_n_levels(hgrid, vgrid, output)

if __name__ == "__main__":
    create_mesh_n_levels_cli()
