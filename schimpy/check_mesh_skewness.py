#!/usr/bin/env python
#! -*- coding: utf-8 -*-
"""Calculate skewness of elements form gr3 file"""

import schimpy.schism_mesh
import schimpy.sms2gr3
import click
import os
import numpy as np

angle_e = np.array((np.pi / 3.0, np.pi * 2.0 / 3.0))


def calculate_skewness(mesh, normalize=True, mask_tri=False):
    """Calculate skewness of elements

    Parameters
    ----------
    mesh: schism_mesh
    normalize: bool, optional
        Normalize skewness if True.
    mask_tri: bool, optional
        If true, mask skewness of triangular elements with zero
    """
    skewness = []
    angle_min = np.pi
    angle_max = -np.pi
    for elem in mesh.elems:
        if mask_tri is True and len(elem) == 3:
            skewness.append(0.0)
            continue
        for i in range(len(elem)):
            node_a = elem[i - 1]
            node_b = elem[i]
            node_c = elem[(i + 1) % len(elem)]
            v_a = mesh.node(node_a)[:2] - mesh.node(node_b)[:2]
            v_b = mesh.node(node_c)[:2] - mesh.node(node_b)[:2]
            angle = np.arccos(
                np.dot(v_a, v_b) / (np.linalg.norm(v_a) * np.linalg.norm(v_b))
            )
            if angle < angle_min:
                angle_min = angle
            if angle > angle_max:
                angle_max = angle
        if angle_min == 0.0:
            raise ValueError("Minimum angle of a polygon is zero")
        if normalize is True:
            e = angle_e[3 - len(elem)]
            skewness.append(np.max((angle_max - e) / (np.pi - e), (e - angle_min) / e))
        else:
            skewness.append(angle_max / angle_min)
    return np.array(skewness)


def write_prop(fpath, data):
    """Write a prop file with a prop vector"""
    output = np.hstack(
        (
            np.arange(data.shape[0], dtype=np.int32).reshape(-1, 1) + 1,
            data.reshape(-1, 1),
        )
    )
    np.savetxt(fpath, output, fmt="%d %g")


def check_skewness(meshfile, propfile, mask_tri, normalize):
    """Calculate skewness of elements from gr3 file."""
    fpath_mesh = meshfile
    fpath_prop = propfile
    if not os.path.exists(fpath_mesh):
        raise ValueError("Mesh file not found")
    if fpath_mesh.endswith(".2dm"):
        fpath_out = "temp.gr3"
        sms2gr3.convert_2dm(fpath_mesh, outfile=fpath_out)
        fpath_mesh = fpath_out
    mesh = schism_mesh.SchismMeshIoFactory().get_reader("gr3").read(fpath_mesh)
    skewness = calculate_skewness(mesh, normalize=normalize, mask_tri=mask_tri)
    write_prop(fpath_prop, skewness)


@click.command()
@click.option(
    "--meshfile",
    default="hgrid.gr3",
    help="Mesh file to analyze. If the file is 2dm format, it will be converted into temp.gr3 and processed.",
)
@click.option(
    "--propfile",
    default="skewness.prop",
    help="Prop output file.",
)
@click.option(
    "--mask_tri",
    is_flag=True,
    help="Mask triangular elements with zero values.",
)
@click.option(
    "--normalize",
    is_flag=True,
    help="Normalize skewness between 0 and 1.",
)
def check_skewness_cli(meshfile, propfile, mask_tri, normalize):

    check_skewness(meshfile, propfile, mask_tri, normalize)


if __name__ == "__main__":
    check_skewness_cli()
