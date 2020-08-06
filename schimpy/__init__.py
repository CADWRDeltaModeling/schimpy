"""Top-level package for schimpy."""

__author__ = """Eli Ateljevich, Kijin Nam"""
__email__ = 'Eli.Ateljevich@water.ca.gov; Kijin.Nam@water.ca.gov'

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .schism_mesh import *
from .prepare_schism import *
from .split_quad import *
from .metricsplot import *
from .batch_metrics import *
