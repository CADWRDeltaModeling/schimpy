"""Top-level package for schimpy."""

__author__ = """Eli Ateljevich, Kijin Nam"""
__email__ = "Eli.Ateljevich@water.ca.gov; Kijin.Nam@water.ca.gov"

try:
    from ._version import __version__
except ImportError:
    __version__ = "0.0.0"  # fallback for weird dev cases
