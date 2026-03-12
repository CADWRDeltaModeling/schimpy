from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from schimpy import sav_density_raster


class DummyMesh:
    def __init__(self, nodes):
        self.nodes = np.asarray(nodes, dtype=float)


@pytest.fixture(scope="module")
def test_rasters(tmp_path_factory):
    d = tmp_path_factory.mktemp("sav_raster_data")

    transform = from_origin(0.0, 10.0, 1.0, 1.0)

    type_arr = np.array(
        [
            [0, 0, 1, 1, 1],
            [0, 1, 1, 2, 2],
            [1, 1, 2, 2, 2],
            [1, 2, 2, 3, 3],
            [2, 2, 3, 3, 3],
        ],
        dtype=np.int16,
    )

    dens_arr = np.array(
        [
            [0, 0, 1, 1, 1],
            [0, 1, 2, 2, 2],
            [1, 2, 2, 3, 3],
            [1, 2, 3, 3, 3],
            [1, 2, 3, 3, 3],
        ],
        dtype=np.int16,
    )

    type_file = d / "test_sav_type.tiff"
    dens_file = d / "test_sav_density.tiff"

    with rasterio.open(
        type_file,
        "w",
        driver="GTiff",
        height=type_arr.shape[0],
        width=type_arr.shape[1],
        count=1,
        dtype=type_arr.dtype,
        transform=transform,
        nodata=None,
    ) as dst:
        dst.write(type_arr, 1)

    with rasterio.open(
        dens_file,
        "w",
        driver="GTiff",
        height=dens_arr.shape[0],
        width=dens_arr.shape[1],
        count=1,
        dtype=dens_arr.dtype,
        transform=transform,
        nodata=None,
    ) as dst:
        dst.write(dens_arr, 1)

    return type_file, dens_file


@pytest.fixture
def mapping_csv(tmp_path):
    p = tmp_path / "sav_mapping.csv"
    p.write_text(
        "type,density_bin,density\n"
        "0,0,0\n"
        "0,1,0\n"
        "0,2,0\n"
        "0,3,0\n"
        "1,0,0\n"
        "1,1,5\n"
        "1,2,10\n"
        "1,3,25\n"
        "2,0,0\n"
        "2,1,4\n"
        "2,2,10\n"
        "2,3,15\n"
        "3,0,0\n"
        "3,1,5\n"
        "3,2,10\n"
        "3,3,30\n"
    )
    return p


def test_center_sampling(test_rasters):
    type_tif, dens_tif = test_rasters

    xs = np.array([2.5])
    ys = np.array([7.5])

    t, d = sav_density_raster.sample_center_pixels(xs, ys, type_tif, dens_tif)

    assert t[0] == 2
    assert d[0] == 2


def test_stencil_shape(test_rasters):
    type_tif, dens_tif = test_rasters

    xs = np.array([2.5])
    ys = np.array([7.5])

    type_s, dens_s = sav_density_raster.sample_stencil(
        xs, ys, type_tif, dens_tif, stencil=3
    )

    assert type_s.shape == (1, 9)
    assert dens_s.shape == (1, 9)


def test_stencil_values_match_manual_extract(test_rasters):
    type_tif, dens_tif = test_rasters

    xs = np.array([2.5])
    ys = np.array([7.5])

    type_s, dens_s = sav_density_raster.sample_stencil(
        xs, ys, type_tif, dens_tif, stencil=3
    )

    expected_type = np.array([1, 1, 2, 1, 2, 2, 2, 2, 3], dtype=np.int64)
    expected_dens = np.array([1, 2, 2, 2, 2, 3, 2, 3, 3], dtype=np.int64)

    assert np.array_equal(type_s[0], expected_type)
    assert np.array_equal(dens_s[0], expected_dens)


def test_like_mask_has_at_least_one_like_pixel(test_rasters, mapping_csv):
    type_tif, dens_tif = test_rasters

    xs = np.array([2.5])
    ys = np.array([7.5])
    z = np.array([1.0])

    dens, dbg = sav_density_raster.compute_density(
        xs,
        ys,
        z,
        type_tif,
        dens_tif,
        mapping=mapping_csv,
        stencil=3,
        depth_limit=3.0,
        strict=True,
        return_debug=True,
    )

    center_type = dbg["type_center"]
    center_bin = dbg["bin_center"]
    like = dbg["like_mask"]

    active = (center_type > 0) & (center_bin > 0)
    assert active[0]
    assert like.sum(axis=1)[0] >= 1
    assert np.isfinite(dens[0])


def test_depth_limit_zeroes_density(test_rasters, mapping_csv):
    type_tif, dens_tif = test_rasters

    xs = np.array([2.5])
    ys = np.array([7.5])
    z = np.array([5.0])

    dens = sav_density_raster.compute_density(
        xs,
        ys,
        z,
        type_tif,
        dens_tif,
        mapping=mapping_csv,
        stencil=3,
        depth_limit=3.0,
    )

    assert dens[0] == 0.0


def test_full_pipeline_positive_density(test_rasters, mapping_csv):
    type_tif, dens_tif = test_rasters

    xs = np.array([2.5])
    ys = np.array([7.5])
    z = np.array([1.0])

    dens = sav_density_raster.compute_density(
        xs,
        ys,
        z,
        type_tif,
        dens_tif,
        mapping=mapping_csv,
        stencil=3,
        depth_limit=3.0,
    )

    assert dens.shape == (1,)
    assert dens[0] > 0.0


def test_yaml_friendly_wrapper(test_rasters, mapping_csv):
    type_tif, dens_tif = test_rasters

    nodes = np.array(
        [
            [2.5, 7.5, 1.0],
            [2.5, 7.5, 5.0],
        ],
        dtype=float,
    )
    mesh = DummyMesh(nodes)
    nodes_sel = np.array([0, 1], dtype=int)

    dens = sav_density_raster.sav_density(
        mesh,
        nodes_sel,
        type_tif,
        dens_tif,
        stencil=3,
        depth_limit=3.0,
        mapping=mapping_csv,
        strict=True,
    )

    assert dens.shape == (2,)
    assert dens[0] > 0.0
    assert dens[1] == 0.0

def test_compute_height_sav_fractions():
    z = np.array([1.0, 2.0, 3.0], dtype=float)
    center_type = np.array([1, 1, 1], dtype=int)
    center_bin = np.array([1, 2, 3], dtype=int)

    h = sav_density_raster.compute_height(
        z,
        center_type,
        center_bin,
        sav_fractions=(0.0, 0.3, 0.65, 0.9),
    )

    assert np.allclose(h, [0.3, 1.3, 2.7])


def test_compute_height_floating_and_emergent():
    z = np.array([1.0, 2.0, 3.0], dtype=float)
    center_type = np.array([2, 3, 0], dtype=int)
    center_bin = np.array([1, 2, 3], dtype=int)

    h = sav_density_raster.compute_height(
        z,
        center_type,
        center_bin,
        floating_height=100.0,
        emergent_height=100.0,
    )

    assert np.allclose(h, [100.0, 100.0, 0.0])


def test_sav_height_wrapper(test_rasters):
    type_tif, dens_tif = test_rasters

    nodes = np.array(
        [
            [0.5, 9.5, 2.0],  # type 0, bin 0 -> 0
            [2.5, 9.5, 2.0],  # type 1, bin 1 -> 0.3 * 2.0
            [3.5, 8.5, 2.0],  # type 2, bin 2 -> 100
            [3.5, 6.5, 2.0],  # type 3, bin 3 -> 100
        ],
        dtype=float,
    )
    mesh = DummyMesh(nodes)
    nodes_sel = np.arange(nodes.shape[0], dtype=int)

    h = sav_density_raster.sav_height(
        mesh,
        nodes_sel,
        type_tif,
        dens_tif,
        sav_fractions=(0.0, 0.3, 0.65, 0.9),
        floating_height=100.0,
        emergent_height=100.0,
        strict=True,
    )

    assert np.allclose(h, [0.0, 0.6, 100.0, 100.0])


def test_sav_node_classes_wrapper(test_rasters):
    type_tif, dens_tif = test_rasters

    nodes = np.array(
        [
            [2.5, 7.5, 1.0],
        ],
        dtype=float,
    )
    mesh = DummyMesh(nodes)
    nodes_sel = np.array([0], dtype=int)

    t, b = sav_density_raster.sav_node_classes(
        mesh,
        nodes_sel,
        type_tif,
        dens_tif,
        strict=True,
    )

    assert t.shape == (1,)
    assert b.shape == (1,)
    assert t[0] == 2
    assert b[0] == 2    