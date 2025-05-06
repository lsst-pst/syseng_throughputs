"""Module for system engineering throughputs curves"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("syseng_throughputs")
except PackageNotFoundError:
    # package is not installed
    pass


from .bandpassUtils import *
from .m5Utils import *
from .sedUtils import *
