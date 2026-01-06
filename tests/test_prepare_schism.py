#!/usr/bin/env python
# -*- coding: utf-8 -*-
from schimpy.prepare_schism import prepare_schism_cli
import unittest
import os
from click.testing import CliRunner


class TestPrepareSchism(unittest.TestCase):
    def setUp(self):
        dir_cur = os.path.dirname(__file__)
        dir_test = os.path.join(dir_cur, "testdata/prepare_schism")
        os.chdir(dir_test)
        self.runner = CliRunner()

    def tearDown(self):
        os.chdir(os.path.dirname(__file__))

    def test_prepare_schism(self):
        os.chdir('simple_triquad')
        result = self.runner.invoke(prepare_schism_cli, ['main_good.yaml'])
        self.assertEqual(result.exit_code, 0, f"Command failed with output: {result.output}")


if __name__ == '__main__':
    unittest.main()
