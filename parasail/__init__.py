import ctypes
import platform
import os
import sys

import numpy

__version__ = "1.1.17"
__title__ = "parasail"
__description__ = "pairwise sequence alignment library"
__uri__ = "https://github.com/jeffdaily/parasail-python"
__author__ = "Jeff Daily"
__email__ = "jeff.daily@pnnl.gov"
__license__ = "BSD"
__copyright__ = "Copyright (c) 2016 Jeff Daily"

# we need the parasail library loaded so we can query the version

_libname = "libparasail.so"
if platform.system() == 'Darwin':
    _libname = "libparasail.dylib"
elif platform.system() == 'Windows':
    _libname = "parasail.dll"
_libpath = os.path.join(os.path.dirname(__file__), _libname)

_lib = None
if os.path.exists(_libpath):
    _lib = ctypes.CDLL(_libpath)
else:
    _lib = ctypes.CDLL(_libname)

c_int_p = ctypes.POINTER(ctypes.c_int)

_lib.parasail_version.argtypes = [c_int_p, c_int_p, c_int_p]
_lib.parasail_version.restype = None

def version():
    major = ctypes.c_int()
    minor = ctypes.c_int()
    patch = ctypes.c_int()
    _lib.parasail_version(
            ctypes.byref(major),
            ctypes.byref(minor),
            ctypes.byref(patch))
    return major.value, minor.value, patch.value

major,minor,patch = version()

# now that we know the version, import the correct bindings
if major == 1:
    from parasail.bindings_v1 import *
else:
    from parasail.bindings_v2 import *

