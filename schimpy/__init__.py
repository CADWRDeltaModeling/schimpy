"""Top-level package for schimpy."""

__author__ = """Eli Ateljevich, Kijin Nam"""
__email__ = "Eli.Ateljevich@water.ca.gov; Kijin.Nam@water.ca.gov"

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("schimpy")
except PackageNotFoundError:
    __version__ = "0.0.0"  # Default version if package metadata is unavailable
