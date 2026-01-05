#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import click
from schimpy.schism_mesh import read_mesh


@click.command()
@click.option('--input_mesh', type=str, default=None, help='Input mesh')
@click.option('--warn', type=float, default=10.0, help='Threshold for warning (areas smaller)')
@click.option('--fail', type=float, default=4.0, help='Threshold for failure (areas smaller)')
def small_areas_cli(input_mesh, warn, fail):
    """Identify small areas"""
    print(input_mesh)
    mesh = read_mesh(input_mesh)
    small_areas(mesh, warn, fail)


def small_areas(
    mesh,
    warn=10.0,
    fail=1.0,
    logger=None,
    write_out_smalls=False,
    prepro_output_dir="./",
):
    if logger is not None:
        logger.info(
            "Checking for elements with small areas. Thresholds: warn={}, fail={}".format(
                warn, fail
            )
        )
    if isinstance(mesh, str):
        mesh = read_mesh(mesh)
    areas = mesh.areas()
    centroids = mesh.centroids()

    sorted_order = np.argsort(areas)
    areas_sorted = areas[sorted_order]
    centroids_sorted = centroids[sorted_order]
    (small,) = np.where(areas_sorted < warn)
    warnings = []
    for s in small:
        x, y = centroids_sorted[s]
        warnings.append(
            "Global element: {:>9d} Area: {:.3f} Centroid: {:.1f},{:.1f}".format(
                sorted_order[s] + 1, areas_sorted[s], x, y
            )
        )

    for item in warnings:
        if logger is None:
            print(item)
        else:
            logger.warning(item)

    if areas_sorted[0] < fail:
        if write_out_smalls:
            import geopandas as gpd
            from shapely.geometry import Point

            sm_centroids = centroids[areas < fail]
            sm_areas = areas[areas < fail]

            gdf = gpd.GeoDataFrame(
                {
                    "FID": range(len(sm_centroids)),
                    "x": sm_centroids[:, 0],
                    "y": sm_centroids[:, 1],
                    "area": sm_areas,
                },
                geometry=[Point(xy) for xy in sm_centroids],
                crs="EPSG:4326",  # Change if you know the mesh CRS
            )
            gdf.to_file(f"{prepro_output_dir}/small_area_centroids.shp")

        raise ValueError(
            "Mesh contains areas smaller than the failure threshold. Consult the log or printout above for areas and warnings"
        )

    return


if __name__ == "__main__":
    small_areas_cli()
