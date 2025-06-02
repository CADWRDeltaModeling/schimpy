#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pytest tests of customized YAML for SCHISM with including and substitution"""
# from __future__ import absolute_import
from schimpy.schism_yaml import load, YamlAction, ArgumentParserYaml
import yaml
import argparse
import sys
import pytest
import copy
import os

# Path to test data directory
dir_testdata = os.path.join(os.path.dirname(__file__), "testdata/schism_yaml")


@pytest.fixture
def argv_backup():
    """Fixture to save and restore sys.argv"""
    original_argv = copy.deepcopy(sys.argv)
    sys.argv = sys.argv[:1]
    yield
    sys.argv = original_argv


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


def test_argparse_action(argv_backup):
    """Test YAML for argparse"""
    fname = os.path.join(dir_testdata, "argparse.yaml")
    # YAML test
    sys.argv.extend(("--yaml", "%s" % fname))
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument("--input1")
    parser.add_argument("--input2")
    parser.add_argument("--input_multi", nargs="+")
    parser.add_argument("--start")
    parser.add_argument("--yaml", action=YamlAction)
    args, _ = parser.parse_known_args()
    assert args.input1 == "1"
    assert args.start == "2009/08/23"
    assert args.input_multi == ["1", "2"]


def test_argparse_action_overwrite1(argv_backup):
    """Test YAML for argparse with overwriting (YAML first)

    When a parameter is specified in both the YAML file and the command line,
    the command line argument should take precedence only when the argument
    comes after the yaml file.
    """
    fname = os.path.join(dir_testdata, "argparse.yaml")
    sys.argv.extend(("--yaml", "%s" % fname))
    sys.argv.extend(("--input1", "2"))
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument("--input1")
    parser.add_argument("--yaml", action=YamlAction)
    args, _ = parser.parse_known_args()
    assert args.input1 == "2"


def test_argparse_action_overwrite2(argv_backup):
    """Test YAML for argparse with overwriting (YAML last)

    When a parameter is specified in both the YAML file and the command line,
    the yaml file argument should take precedence when the argument comes
    before the yaml file.
    """
    fname = os.path.join(dir_testdata, "argparse.yaml")
    sys.argv.extend(("--input1", "2"))
    sys.argv.extend(("--yaml", "%s" % fname))
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument("--input1")
    parser.add_argument("--yaml", action=YamlAction)
    args, _ = parser.parse_known_args()
    assert args.input1 == "1"


def test_argparse_action_whitespace(argv_backup):
    """Test YAML for argparse with whitespace in key"""
    fname = os.path.join(dir_testdata, "argparse_whitespace.yaml")
    sys.argv.extend(("--yaml", "%s" % fname))
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument("--input1")
    parser.add_argument("--yaml", action=YamlAction)
    with pytest.raises(ValueError):
        parser.parse_args()


def test_argparse_action_multiple_level(argv_backup):
    """Test YAML for argparse with multiple level"""
    fname = os.path.join(dir_testdata, "argparse_multiple_level.yaml")
    sys.argv.extend(("--yaml", "%s" % fname))
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument("--input1")
    parser.add_argument("--yaml", action=YamlAction)
    with pytest.raises(ValueError):
        parser.parse_args()


def test_multiple_argparse_action(argv_backup):
    """Test multiple YAML for argparse"""
    fname1 = os.path.join(dir_testdata, "argparse.yaml")
    fname2 = os.path.join(dir_testdata, "argparse2.yaml")
    # YAML test
    sys.argv.extend(("--yaml", fname1, "--yaml", fname2))
    parser = argparse.ArgumentParser()
    parser.add_argument("--input1")
    parser.add_argument("--input2")
    parser.add_argument("--start")
    parser.add_argument("--yaml", action=YamlAction)
    args, _ = parser.parse_known_args()
    assert args.input1 == "4"
    assert args.input2 == "2"
    assert args.start == "2009/08/23"


def test_argparse_yaml(argv_backup):
    """Test ArgumentParserYaml"""
    fname = os.path.join(dir_testdata, "argparse.yaml")
    parser = ArgumentParserYaml(fromfile_prefix_chars="@")
    # YAML test
    sys.argv.append(("@%s" % fname))
    parser.add_argument("--input1")
    parser.add_argument("--input2")
    parser.add_argument("--input_multi", nargs="+")
    parser.add_argument("--start")
    # args = parser.parse_args()  # Left these for reference
    args, _ = parser.parse_known_args()  # Use this to ignore some arguments
    assert args.input1 == "1"
    assert args.start == "2009/08/23"
    assert args.input_multi == ["1", "2"]


def test_argparse_yaml_overwrite1(argv_backup):
    """Test ArgumentParserYaml with overwriting (YAML first)"""
    fname = os.path.join(dir_testdata, "argparse.yaml")
    parser = ArgumentParserYaml(fromfile_prefix_chars="@")
    sys.argv.append(("@%s" % fname))
    sys.argv.extend(("--input1", "2"))
    parser.add_argument("--input1")
    args, _ = parser.parse_known_args()
    assert args.input1 == "2"


def test_argparse_yaml_overwrite2(argv_backup):
    """Test ArgumentParserYaml with overwriting (YAML last)"""
    fname = os.path.join(dir_testdata, "argparse.yaml")
    parser = ArgumentParserYaml(fromfile_prefix_chars="@")
    sys.argv.extend(("--input1", "2"))
    sys.argv.append(("@%s" % fname))
    parser.add_argument("--input1")
    args, _ = parser.parse_known_args()
    assert args.input1 == "1"


def test_argparse_yaml_whitespace(argv_backup):
    """Test ArgumentParserYaml with whitespace in key"""
    fname = os.path.join(dir_testdata, "argparse_whitespace.yaml")
    parser = ArgumentParserYaml(fromfile_prefix_chars="@")
    sys.argv.append(("@%s" % fname))
    parser.add_argument("--input1")
    with pytest.raises(ValueError):
        parser.parse_args()


def test_argparse_yaml_multiple_level(argv_backup):
    """Test ArgumentParserYaml with multiple level"""
    fname = os.path.join(dir_testdata, "argparse_multiple_level.yaml")
    parser = ArgumentParserYaml(fromfile_prefix_chars="@")
    sys.argv.append(("@%s" % fname))
    parser.add_argument("--input1")
    with pytest.raises(ValueError):
        parser.parse_args()


def test_multiple_argparse_yaml(argv_backup):
    """Test multiple YAML for ArgumentParserYaml"""
    fname1 = os.path.join(dir_testdata, "argparse.yaml")
    fname2 = os.path.join(dir_testdata, "argparse2.yaml")
    # YAML test
    sys.argv.extend(("@%s" % fname1, "@%s" % fname2))
    parser = ArgumentParserYaml(fromfile_prefix_chars="@")
    parser.add_argument("--input1")
    parser.add_argument("--input2")
    parser.add_argument("--start")
    args, _ = parser.parse_known_args()
    assert args.input1 == "4"  # Overwrite with the latter value
    assert args.input2 == "2"
    assert args.start == "2009/08/23"


def test_envvar_in_argument(argv_backup):
    """Test environment variables from the command line

    When a config variable is given in the command line, it should be
    substituted into the YAML file.
    """
    fname = os.path.join(dir_testdata, "env_wo_sub.yaml")

    with open(fname, "r") as f:
        params = load(f, envvar={"var5": "/home"})

    assert params["mesh"]["b"] == "/home/baydelta72.gr3"


def test_envvar_in_argparse_argument(argv_backup):
    """Test environment variables from the command line

    When a config variable is given in the command line, it should be
    substituted into the YAML file.
    """
    fname = os.path.join(dir_testdata, "env_wo_sub.yaml")

    sys.argv.append((fname))
    sys.argv.extend(("--var5", "/home"))
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args, extra = parser.parse_known_args()

    def parse_extra_args(extra_args):
        """
        Converts a list like ['--a', '1', '--b=2', '--c', '3']
        into a dictionary: {'a': '1', 'b': '2', 'c': '3'}
        Raises ValueError for unpaired or malformed args.
        """
        result = {}
        i = 0
        while i < len(extra_args):
            arg = extra_args[i]
            if arg.startswith("--"):
                if "=" in arg:
                    key, val = arg[2:].split("=", 1)
                    result[key] = val
                    i += 1
                else:
                    if i + 1 >= len(extra_args) or extra_args[i + 1].startswith("--"):
                        raise ValueError(f"Missing value for argument: {arg}")
                    key = arg[2:]
                    val = extra_args[i + 1]
                    result[key] = val
                    i += 2
            else:
                raise ValueError(f"Unexpected argument format: {arg}")
        return result

    envvar = parse_extra_args(extra)
    with open(args.filename, "r") as f:
        params = load(f, envvar=envvar)

    assert params["mesh"]["b"] == "/home/baydelta72.gr3"
