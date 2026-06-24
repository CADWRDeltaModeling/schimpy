#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regression tests for vgrid_nlayer_to_gr3."""

import os
import tempfile
import unittest

import numpy as np
from click.testing import CliRunner
from schimpy.schism_mesh import read_mesh
from schimpy.vgrid_nlayer_to_gr3 import vgrid_nlayer_to_gr3, vgrid_nlayer_to_gr3_cli


class TestVgridNLayerToGr3(unittest.TestCase):

    def setUp(self):
        self.dir_cur = os.path.dirname(__file__)
        self.testdata_dir = os.path.join(self.dir_cur, "testdata", "testmesh")
        self.hgrid = os.path.join(self.testdata_dir, "testmesh.gr3")
        self.vgrid = os.path.join(self.testdata_dir, "vgrid.in")

    def test_vgrid_nlayer_to_gr3_writes_expected_outputs(self):
        mesh = read_mesh(self.hgrid, self.vgrid, vgrid_version="5.8")
        expected_nlevels = mesh.vmesh.n_vert_levels() - np.asarray(mesh.vmesh.kbps)
        expected_nlayers = np.maximum(expected_nlevels - 1, 0)

        with tempfile.TemporaryDirectory() as tmpdir:
            nlevel_path = os.path.join(tmpdir, "nlevels.gr3")
            nlayer_path = os.path.join(tmpdir, "nlayers.gr3")

            nlevels, nlayers = vgrid_nlayer_to_gr3(
                self.hgrid,
                self.vgrid,
                output_nlevel=nlevel_path,
                output_nlayer=nlayer_path,
                vgrid_version="5.8",
            )

            self.assertTrue(os.path.exists(nlevel_path))
            self.assertTrue(os.path.exists(nlayer_path))
            np.testing.assert_array_equal(nlevels, expected_nlevels)
            np.testing.assert_array_equal(nlayers, expected_nlayers)

            mesh_nlevel = read_mesh(nlevel_path)
            mesh_nlayer = read_mesh(nlayer_path)
            self.assertEqual(mesh_nlevel.n_nodes(), mesh.n_nodes())
            self.assertEqual(mesh_nlayer.n_nodes(), mesh.n_nodes())
            np.testing.assert_array_equal(mesh_nlevel.nodes[:, 2], expected_nlevels)
            np.testing.assert_array_equal(mesh_nlayer.nodes[:, 2], expected_nlayers)

    def test_vgrid_nlayer_to_gr3_cli(self):
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            nlevel_path = os.path.join(tmpdir, "nlevels.gr3")

            result = runner.invoke(
                vgrid_nlayer_to_gr3_cli,
                [
                    "--hgrid",
                    self.hgrid,
                    "--vgrid",
                    self.vgrid,
                    "--vgrid_version",
                    "5.8",
                    "--nlevel",
                    nlevel_path,
                ],
            )
            self.assertEqual(result.exit_code, 0, result.output)
            self.assertTrue(os.path.exists(nlevel_path))

    def test_vgrid_nlayer_to_gr3_requires_output(self):
        with self.assertRaises(ValueError):
            vgrid_nlayer_to_gr3(self.hgrid, self.vgrid)


if __name__ == "__main__":
    unittest.main()
