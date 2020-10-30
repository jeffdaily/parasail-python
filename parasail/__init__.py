import ctypes
import platform
import os
import re
import sys

import numpy

__version__ = "1.2.1c1"
__title__ = "parasail"
__description__ = "pairwise sequence alignment library"
__uri__ = "https://github.com/jeffdaily/parasail-python"
__author__ = "Jeff Daily"
__email__ = "jeffrey.daily@gmail.com"
__license__ = "BSD"
__copyright__ = "Copyright (c) 2016 Jeff Daily"

# we need the parasail library loaded so we can query the version

_libname = "libparasail.so"
if platform.system() == 'Darwin':
    _libname = "libparasail.dylib"
elif platform.system() == 'Windows':
    _libname = "parasail.dll"
_libpath = os.path.join(os.path.dirname(__file__), _libname)

_verbose = os.environ.get("PARASAIL_VERBOSE", False)

# attempt to load library from package location
_lib = None
try:
    _lib = ctypes.CDLL(_libpath)
except:
    if _verbose: print("failed to open '{}' in package location".format(_libname))

if not _lib:
    try:
        _lib = ctypes.CDLL(_libname)
    except:
        if _verbose: print("failed to open '{}' using ctypes.CDLL defaults".format(_libname))

def _lib_search():
    global _lib
    global _libpath
    _libpath = None
    for env_var in ["PARASAIL_LIBPATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "PATH"]:
        if env_var in os.environ:
            if _verbose: print("searching {} for {}".format(env_var,_libname))
            for path in re.split("[:;]", os.environ[env_var]):
                if _verbose: print("searching {}".format(path))
                for attempt in [os.path.join(root,_libname) for root,dirs,files in os.walk(path) if _libname in files]:
                    try:
                        _lib = ctypes.CDLL(attempt)
                        _libpath = attempt
                        # shortcut the search during bindings import
                        os.environ["PARASAIL_LIBPATH"] = os.path.dirname(attempt)
                        if _verbose: print("found '{}'".format(attempt))
                        return
                    except:
                        if _verbose: print("attempted but failed to load '{}'".format(attempt))
        else:
            if _verbose: print("env var {} not set".format(env_var))

# if library load failed, search for it
if not _lib:
    _lib_search()

# library load still failed, hard error
if not _lib:
    raise Exception("failed to locate and open '{}'".format(_libname))

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

