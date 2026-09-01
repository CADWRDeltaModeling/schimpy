#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pytest tests of customized YAML for SCHISM with including and substitution"""
# from __future__ import absolute_import
from schimpy.schism_yaml import load
from click.testing import CliRunner
from schimpy.schism_yaml import schism_yaml_cli
import yaml
import pytest
import os

# Path to test data directory
dir_testdata = os.path.join(os.path.dirname(__file__), "testdata/schism_yaml")


def test_unsubstituted_env():
    """Should raise an exception with unsubstituted environment variables"""
    fname = os.path.join(dir_testdata, "env_wo_sub.yaml")
    with open(fname, "r") as fin:
        with pytest.raises(yaml.composer.ComposerError):
            load(fin)


def test_multiple_level_in_env():
    """Should raise an exception with multiple level variables in env"""
    fname = os.path.join(dir_testdata, "env_w_multiplelevel.yaml")
    with open(fname, "r") as fin:
        with pytest.raises(ValueError):
            load(fin)


def test_multiple_substitution_in_one_item():
    """Test multiple substitution in one item"""
    fname = os.path.join(dir_testdata, "substitute_multiple_times_in_one_item.yaml")
    with open(fname, "r") as f:
        data = load(f)
        assert data["mesh"]["work_dir"] == "dir1/file1"


def test_substitute_wo_including():
    """Test a simple substitution without include"""
    fname = os.path.join(dir_testdata, "substitute_wo_including.yaml")
    with open(fname, "r") as f:
        data = load(f)
        assert data["mesh"]["work_dir"] == "dir1"


def test_substitute_multiple_including():
    """Test substitution in multiple including"""
    fname = os.path.join(dir_testdata, "substitute_w_including.yaml")
    with open(fname, "r") as f:
        data = load(f)
        assert data["mesh"]["work_dir"] == "dir1"

    fname = os.path.join(dir_testdata, "substitute_w_including_typo.yaml")
    with open(fname, "r") as f:
        with pytest.raises(ValueError):
            load(f)


def test_envvar_in_argument():
    """Test environment variables from the command line

    When a config variable is given in the command line, it should be
    substituted into the YAML file.
    """
    fname = os.path.join(dir_testdata, "env_wo_sub.yaml")

    with open(fname, "r") as f:
        params = load(f, envvar={"var5": "/home"})

    assert params["mesh"]["b"] == "/home/baydelta72.gr3"


def test_var_in_include_filename():
    """Test that config variables can be substituted into include filenames"""
    fname = os.path.join(dir_testdata, "substitute_w_var_include.yaml")
    with open(fname, "r") as f:
        data = load(f)
        assert data["mesh"]["work_dir"] == "dir1"


def test_click_cli():
    """Test the click CLI"""
    fname = os.path.join(dir_testdata, "substitute_wo_including.yaml")
    runner = CliRunner()
    result = runner.invoke(schism_yaml_cli, [fname])
    assert result.exit_code == 0
    # Check that the output contains the expected YAML content
    assert "work_dir: dir1" in result.output


def test_click_cli_set():
    """Test the --set option on the CLI injects substitution variables"""
    fname = os.path.join(dir_testdata, "substitute_w_set.yaml")
    runner = CliRunner()
    result = runner.invoke(
        schism_yaml_cli,
        ["--set", "mesh_year=2018", "--set", "mesh_version=118", fname],
    )
    assert result.exit_code == 0
    assert "work_dir: 2018/118" in result.output


def test_include_list_concatenates_lists():
    """A list-style include accumulates list sections such as ``structures``"""
    fname = os.path.join(dir_testdata, "merge_agree.yaml")
    with open(fname, "r") as f:
        data = load(f)
    names = [s["name"] for s in data["hydraulics"]["structures"]]
    assert names == ["weir_north", "breach_west"]


def test_include_list_agreeing_scalar_is_not_accumulated():
    """A scalar repeated with the same value passes through untouched.

    The former behavior applied ``+`` to every repeated key, which silently
    doubled ``nudging`` and concatenated ``outputfile`` into nonsense.
    """
    fname = os.path.join(dir_testdata, "merge_agree.yaml")
    with open(fname, "r") as f:
        data = load(f)
    assert data["hydraulics"]["nudging"] == 0.05
    assert data["hydraulics"]["outputfile"] == "hydraulics.in"


def test_include_list_disagreeing_scalar_warns_and_takes_first():
    """A scalar that disagrees keeps the first file's value and warns"""
    fname = os.path.join(dir_testdata, "merge_disagree.yaml")
    with open(fname, "r") as f:
        with pytest.warns(UserWarning, match="nudging"):
            data = load(f)
    assert data["hydraulics"]["nudging"] == 0.05
