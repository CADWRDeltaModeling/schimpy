"""Top-level package for schimpy."""

__author__ = """Eli Ateljevich, Kijin Nam"""
__email__ = "Eli.Ateljevich@water.ca.gov; Kijin.Nam@water.ca.gov"

from schimpy._version import get_versions


def load_setup_py_data():
    return get_versions()["version"]
