#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the LSC2 v2 vgrid generator integration."""

import os
import pytest
import numpy as np


@pytest.fixture
def simple_mesh():
    """Load the simple_triquad test mesh."""
    from schimpy.schism_mesh import read_mesh

    testdir = os.path.join(
        os.path.dirname(__file__),
        "testdata", "prepare_schism", "simple_triquad",
    )
    return read_mesh(os.path.join(testdir, "hgrid.gr3"))


class TestBuildSizefun:
    """Test depth function instantiation from config dicts."""

    def test_bilinear_default(self):
        from schimpy.create_vgrid_lsc2_v2 import _build_sizefun

        section = {}
        sf = _build_sizefun(section)
        assert hasattr(sf, "depth")
        # Default bilinear params
        np.testing.assert_allclose(sf.params, [0.5, 0.005, 0.08, -3.0e-5])

    def test_bilinear_custom_params(self):
        from schimpy.create_vgrid_lsc2_v2 import _build_sizefun

        section = {
            "depth_function": {
                "name": "bilinear",
                "params": {"a": 0.6, "b": 0.01, "c": 0.1, "d": -1e-4},
            }
        }
        sf = _build_sizefun(section)
        np.testing.assert_allclose(sf.params, [0.6, 0.01, 0.1, -1e-4])

    def test_sigma_cap(self):
        from schimpy.create_vgrid_lsc2_v2 import _build_sizefun

        section = {"depth_function": {"name": "sigma_cap", "params": {"dz_max": 2.0}}}
        sf = _build_sizefun(section)
        assert sf.dz_max == 2.0

    def test_unknown_name_raises(self):
        from schimpy.create_vgrid_lsc2_v2 import _build_sizefun

        with pytest.raises(ValueError, match="Unknown depth_function.name"):
            _build_sizefun({"depth_function": {"name": "bogus"}})


class TestBuildPipelineParams:
    """Test pipeline param construction from config dicts."""

    def test_defaults(self):
        from schimpy.create_vgrid_lsc2_v2 import _build_pipeline_params

        pp = _build_pipeline_params({})
        assert pp.L_smooth_passes == 2
        assert pp.hysteresis.delta_up == 0.7
        assert pp.fit.kappa == 0.35
        assert pp.fit.tmin == 0.35

    def test_custom_values(self):
        from schimpy.create_vgrid_lsc2_v2 import _build_pipeline_params

        section = {
            "algorithm": {
                "Lsmooth_passes": 30,
                "Lsmooth_kappa": 0.8,
                "hysteresis": {"delta_up": 1.1, "Jmax": 0},
                "fit": {"kappa": 0.55, "tmin": 0.45, "sm_iters": 2},
            },
            "constraint_taper_rings": 5,
        }
        pp = _build_pipeline_params(section)
        assert pp.L_smooth_passes == 30
        assert pp.L_smooth_kappa == 0.8
        assert pp.hysteresis.delta_up == 1.1
        assert pp.hysteresis.Jmax == 0
        assert pp.fit.kappa == 0.55
        assert pp.fit.tmin == 0.45
        assert pp.fit.sm_iters == 2
        assert pp.constraint_taper_rings == 5


class TestLoadRegionConstraints:
    """Test region constraint loading from polygon dicts."""

    def test_none_returns_none(self):
        from schimpy.create_vgrid_lsc2_v2 import _load_region_constraints

        n_min, n_max = _load_region_constraints(None, None)
        assert n_min is None
        assert n_max is None

    def test_string_raises(self):
        from schimpy.create_vgrid_lsc2_v2 import _load_region_constraints

        with pytest.raises(TypeError, match="must be a dict"):
            _load_region_constraints(None, "some_file.yaml")

    def test_inline_polygons(self, simple_mesh):
        from schimpy.create_vgrid_lsc2_v2 import _load_region_constraints

        # Build a polygon that covers most of the mesh (coords 0-100, 0-100)
        rc = {
            "polygons": [
                {
                    "name": "test_min",
                    "type": "min",
                    "attribute": "5",
                    "vertices": [
                        [0, 0],
                        [100, 0],
                        [100, 100],
                        [0, 100],
                        [0, 0],
                    ],
                }
            ]
        }
        n_min, n_max = _load_region_constraints(simple_mesh, rc)
        # attribute='5' means 5 layers → Nlevels=6
        assert n_min is not None
        assert n_max is None
        assert np.any(n_min == 6)


class TestVgridDispatch:
    """Test the version dispatch logic in create_vgrid."""

    def test_missing_version_raises(self):
        from unittest.mock import MagicMock
        from schimpy.prepare_schism import create_vgrid

        s = MagicMock()
        inputs = {
            "prepro_output_dir": "/tmp/test",
            "vgrid": {
                "vgrid_out": "vgrid.in",
                "eta": 1.5,
                "vgrid_version": "5.10",
            },
        }
        logger = MagicMock()
        with pytest.raises(ValueError, match="vgrid_generator_version"):
            create_vgrid(s, inputs, logger)

    def test_unknown_version_raises(self):
        from unittest.mock import MagicMock
        from schimpy.prepare_schism import create_vgrid

        s = MagicMock()
        inputs = {
            "prepro_output_dir": "/tmp/test",
            "vgrid": {
                "vgrid_generator_version": "v3",
                "vgrid_out": "vgrid.in",
                "eta": 1.5,
                "vgrid_version": "5.10",
            },
        }
        logger = MagicMock()
        with pytest.raises(ValueError, match="Unknown vgrid_generator_version"):
            create_vgrid(s, inputs, logger)

    def test_v1_dispatch(self):
        from unittest.mock import MagicMock, patch
        from schimpy.prepare_schism import create_vgrid

        s = MagicMock()
        inputs = {
            "prepro_output_dir": "/tmp/test",
            "vgrid": {
                "vgrid_generator_version": "v1",
                "vgrid_out": "vgrid.in",
                "minmaxlayerfile": "minmax.shp",
                "eta": 1.5,
                "vgrid_version": "5.10",
            },
        }
        logger = MagicMock()
        with patch("schimpy.prepare_schism.vgrid_gen") as mock_gen:
            create_vgrid(s, inputs, logger)
            mock_gen.assert_called_once()
            # First arg is the mesh
            assert mock_gen.call_args[0][0] is s.mesh

    def test_v2_dispatch(self):
        from unittest.mock import MagicMock, patch
        from schimpy.prepare_schism import create_vgrid

        s = MagicMock()
        inputs = {
            "prepro_output_dir": "/tmp/test",
            "vgrid": {
                "vgrid_generator_version": "v2",
                "vgrid_out": "vgrid.in",
                "eta": 1.5,
                "vgrid_version": "5.10",
                "depth_function": {"name": "bilinear"},
            },
        }
        logger = MagicMock()
        with patch("schimpy.prepare_schism.vgrid_gen_v2") as mock_gen:
            create_vgrid(s, inputs, logger)
            mock_gen.assert_called_once()
            assert mock_gen.call_args[1]["hgrid"] is s.mesh

    def test_v1_warns_on_v2_keys(self):
        from unittest.mock import MagicMock, patch
        from schimpy.prepare_schism import create_vgrid

        s = MagicMock()
        inputs = {
            "prepro_output_dir": "/tmp/test",
            "vgrid": {
                "vgrid_generator_version": "v1",
                "vgrid_out": "vgrid.in",
                "minmaxlayerfile": "minmax.shp",
                "eta": 1.5,
                "vgrid_version": "5.10",
                "depth_function": {"name": "bilinear"},  # v2-only
            },
        }
        logger = MagicMock()
        with patch("schimpy.prepare_schism.vgrid_gen"):
            create_vgrid(s, inputs, logger)
            logger.warning.assert_called()
            warn_msg = logger.warning.call_args[0][0]
            assert "v2-only" in warn_msg
