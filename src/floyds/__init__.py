from floyds.util import *
import floyds.cosmics
from floyds.floydsspecdef import *
from floyds.floydsspecauto import *

__version__ = "unknown"
try:
    from _version import __version__
except ImportError:
    # We're running in a tree that doesn't have a _version.py, so we don't know what our version is.
    pass
