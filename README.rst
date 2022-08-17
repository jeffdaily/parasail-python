parasail-python
===============

Python Bindings for the Parasail C Library

Travis Build Status:

.. image:: https://travis-ci.org/jeffdaily/parasail-python.svg?branch=master
    :alt: Build Status

PyPI Package:

.. image:: https://badge.fury.io/py/parasail.svg
    :target: https://badge.fury.io/py/parasail

Author: Jeff Daily (jeffrey.daily@gmail.com)

Table of Contents
-----------------

-  `Installation <#installation>`__

   -  `Using pip <#using-pip>`__
   -  `Testing <#tesing>`__
   -  `Building from Source <#building-from-source>`__

-  `Quick Example <#quick-example>`__
-  `Standard Function Naming Convention <#standard-function-naming-convention>`__
-  `Profile Function Naming Convention <#profile-function-naming-convention>`__
-  `Substitution Matrices <#substitution-matrices>`__
-  `SSW Library Emulation <#ssw-library-emulation>`__
-  `Banded Global Alignment <#banded-global-alignment>`__
-  `File Input <#file-input>`__
-  `Tracebacks <#tracebacks>`__
-  `Citing parasail <#citing-parasail>`__
-  `License: Battelle BSD-style <#license-battelle-bsd-style>`__

This package contains Python bindings for
`parasail <https://github.com/jeffdaily/parasail>`__. Parasail is a SIMD
C (C99) library containing implementations of the Smith-Waterman
(local), Needleman-Wunsch (global), and semi-global pairwise sequence
alignment algorithms.

Installation
------------

`back to top <#table-of-contents>`__

Using pip
+++++++++

`back to top <#table-of-contents>`__

The recommended way of installing is to use the latest version available via pip.

::

    pip install parasail
    
Binaries for Windows and OSX should be available via pip.  Using pip on a Linux platform will first download the latest version of the parasail C library sources and then compile them automatically into a shared library.  For an installation from sources, or to learn how the pip installation works on Linux, please read on.

Testing
+++++++

`back to top <#table-of-contents>`__

To run the testsuite use the unittest runner.

::

    python -m unittest discover tests

Building from Source
++++++++++++++++++++

`back to top <#table-of-contents>`__

The parasail python bindings are based on ctypes.  Unfortunately, best practices are not firmly established for providing cross-platform and user-friendly python bindings based on ctypes.  The approach with parasail-python is to install the parasail shared library as "package data" and use a relative path from the parasail/__init__.py in order to locate the shared library.

There are two approaches currently supported.  First, you can compile your own parasail shared library using one of the recommended build processes described in the parasail C library README.md, then copy the parasail.dll (Windows), libparasail.so (Linux), or libparasail.dylib (OSX) shared library to parasail-python/parasail -- the same folder location as parasasail-python/parasail/__init__.py.

The second approach is to let the setup.py script attempt to download and compile the parasail C library for you using the configure script that comes with it.  This happens as a side effect of the bdist_wheel target.

::

    python setup.py bdist_wheel

The bdist_wheel target will first look for the shared library.  If it exists, it will happily install it as package data.  Otherwise, the latest parasail master branch from github will be downloaded, unzipped, configured, made, and the shared library will be copied into the appropriate location for package data installation.

The downloading and building of the parasail C library can be skipped if you set the environment variable PARASAIL_SKIP_BUILD to any value prior to running setup.py or pip install. At runtime during import, the parasail bindings will search for the parasail C library first in the package data location, then in standard system locations, and lastly by searching through the environment variables PARASAIL_LIBPATH, LD_LIBRARY_PATH, DYLD_LIBRARY_PATH, and PATH.. For verbose output during this search, set PARASAIL_VERBOSE=1.

Quick Example
-------------

`back to top <#table-of-contents>`__

The Python interface only includes bindings for the dispatching
functions, not the low-level instruction set-specific function calls.
The Python interface also includes wrappers for the various PAM and
BLOSUM matrices included in the distribution.

Gap open and extension penalties are specified as positive integers.  When any of the algorithms open a gap, only the gap open penalty alone is applied.

.. code:: python

    import parasail
    result = parasail.sw_scan_16("asdf", "asdf", 11, 1, parasail.blosum62)
    result = parasail.sw_stats_striped_8("asdf", "asdf", 11, 1, parasail.pam100)

Be careful using the attributes of the Result object - especially on Result instances constructed on the fly. For example, calling `parasail.sw_trace("asdf", "asdf", 11, 1, parasail.blosum62).cigar.seq` returns a numpy.ndarray that wraps a pointer to memory that is invalid because the Cigar is deallocated before the `seq` statement. You can avoid this problem by assigning Result instances to variables as in the example above.

Standard Function Naming Convention
-----------------------------------

`back to top <#table-of-contents>`__

There are many functions within the parasail library, but most are variations of the familiar main
algorithms.  The following table describes the main algorithms and the shorthand name used for the function.

========================================================================================= =============
Algorithm                                                                                 Function Name
========================================================================================= =============
Smith-Waterman local alignment                                                            sw
Needleman-Wunsch global alignment                                                         nw
Semi-Global, do not penalize gaps at beginning of s1/query                                sg_qb
Semi-Global, do not penalize gaps at end of s1/query                                      sg_qe
Semi-Global, do not penalize gaps at beginning and end of s1/query                        sg_qx
Semi-Global, do not penalize gaps at beginning of s2/database                             sg_db
Semi-Global, do not penalize gaps at end of s2/database                                   sg_de
Semi-Global, do not penalize gaps at beginning and end of s2/database                     sg_dx
Semi-Global, do not penalize gaps at beginning of s1/query and end of s2/database         sg_qb_de
Semi-Global, do not penalize gaps at beginning of s2/database and end of s1/query         sg_qe_db
Semi-Global, do not penalize gaps at beginning of s1/query and beginning of s2/database   sg_qb_db
Semi-Global, do not penalize gaps at end of s2/database and end of s1/query               sg_qe_de
Semi-Global, do not penalize gaps at beginning and end of both sequences                  sg
========================================================================================= =============

A good summary of the various alignment algorithms can be found courtesy of Dr. Dannie Durand's course on
computational genomics `here <http://www.cs.cmu.edu/~durand/03-711/2015/Lectures/PW_sequence_alignment_2015.pdf>`_.
The same document was copied locally to the C library repo in case this link ever breaks (`link <https://github.com/jeffdaily/parasail/blob/master/contrib/PW_sequence_alignment_2015.pdf>`_).

To make it easier to find the function you're looking for, the function names follow a naming convention.  The following will use set notation {} to indicate a selection must be made and brackets [] to indicate an optional part of the name.

- Non-vectorized, reference implementations.

  - Required, select algorithm from table above.
  - Optional return alignment statistics.
  - Optional return DP table or last row/col.
  - Optional use a prefix scan implementation.
  - ``parasail. {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sg_qb_db,sg_qe_de,sw} [_stats] [{_table,_rowcol}] [_scan]``

- Non-vectorized, traceback-capable reference implementations.

  - Required, select algorithm from table above.
  - Optional use a prefix scan implementation.
  - ``parasail. {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sg_qb_db,sg_qe_de,sw} _trace [_scan]``

- Vectorized.

  - Required, select algorithm from table above.
  - Optional return alignment statistics.
  - Optional return DP table or last row/col.
  - Required, select vectorization strategy -- striped is a good place to start, but scan is often faster for global alignment.
  - Required, select solution width. 'sat' will attempt 8-bit solution but if overflow is detected it will then perform the 16-bit operation. Can be faster in some cases, though 16-bit is often sufficient.
  - ``parasail. {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sg_qb_db,sg_qe_de,sw} [_stats] [{_table,_rowcol}] {_striped,_scan,_diag} {_8,_16,_32,_64,_sat}``

- Vectorized, traceback-capable.

  - Required, select algorithm from table above.
  - Required, select vectorization strategy -- striped is a good place to start, but scan is often faster for global alignment.
  - Required, select solution width. 'sat' will attempt 8-bit solution but if overflow is detected it will then perform the 16-bit operation. Can be faster in some cases, though 16-bit is often sufficient.
  - ``parasail. {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sg_qb_db,sg_qe_de,sw} _trace {_striped,_scan,_diag} {_8,_16,_32,_64,_sat}``

Profile Function Naming Convention
----------------------------------

`back to top <#table-of-contents>`__

It has been noted in literature that some performance can be gained by reusing the query sequence when using striped [Farrar, 2007] or scan [Daily, 2015] vector strategies.  There is a special subset of functions that enables this behavior.  For the striped and scan vector implementations *only*, a query profile can be created and reused for subsequent alignments. This can noticeably speed up applications such as database search.

- Profile creation

  - Optional, prepare query profile for a function that returns statistics.  Stats require additional data structures to be allocated.
  - Required, select solution width. 'sat' will allocate profiles for both 8- and 16-bit solutions.
  - ``parasail.profile_create [_stats] {_8,_16,_32,_64,_sat}``

- Profile use

  - Vectorized.

    - Required, select algorithm from table above.
    - Optional return alignment statistics.
    - Optional return DP table or last row/col.
    - Required, select vectorization strategy -- striped is a good place to start, but scan is often faster for global alignment.
    - Required, select solution width. 'sat' will attempt 8-bit solution but if overflow is detected it will then perform the 16-bit operation. Can be faster in some cases, though 16-bit is often sufficient.
    - ``parasail. {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sg_qb_db,sg_qe_de,sw} [_stats] [{_table,_rowcol}] {_striped,_scan} _profile {_8,_16,_32,_64,_sat}``

  - Vectorized, traceback-capable.

    - Required, select algorithm from table above.
    - Required, select vectorization strategy -- striped is a good place to start, but scan is often faster for global alignment.
    - Required, select solution width. 'sat' will attempt 8-bit solution but if overflow is detected it will then perform the 16-bit operation. Can be faster in some cases, though 16-bit is often sufficient.
    - ``parasail. {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sg_qb_db,sg_qe_de,sw} _trace {_striped,_scan} _profile {_8,_16,_32,_64,_sat}``

Please note that the bit size you select for creating the profile *must* match the bit size of the function you call. The example below uses a 16-bit profile and a 16-bit function.

.. code:: python

    profile = parasail.profile_create_16("asdf", parasail.blosum62)
    result1 = parasail.sw_trace_striped_profile_16(profile, "asdf", 10, 1)
    result2 = parasail.nw_scan_profile_16(profile, "asdf", 10, 1)

Substitution Matrices
---------------------

`back to top <#table-of-contents>`__

parasail bundles a number of substitution matrices including PAM and BLOSUM.  To use them, look them up by name (useful for command-line parsing) or use directly. For example

.. code:: python

    print(parasail.blosum62)
    matrix = parasail.Matrix("pam100")

You can also create your own matrices with simple match/mismatch values.
For more complex matrices, you can start by copying a built-in matrix or
start simple and modify values as needed. For example

.. code:: python

    # copy a built-in matrix, then modify like a numpy array
    matrix = parasail.blosum62.copy()
    matrix[2,4] = 200
    matrix[3,:] = 100
    user_matrix = parasail.matrix_create("ACGT", 2, -1)

You can also parse simple matrix files using the function if the file is in the following format::

    #
    # Any line starting with '#' is a comment.
    #
    # Needs a row for the alphabet.  First column is a repeat of the
    # alphabet and assumed to be identical in order to the first alphabet row.
    #
    # Last row and column *must* be a non-alphabet character to represent
    # any input sequence character that is outside of the alphabet.
    #
        A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U   *
    A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2  -4  -5
    T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5  -5
    G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2  -4  -5
    C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2  -4  -5
    S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1  -4  -5
    W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1   1  -5
    R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1  -4  -5
    Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1   1  -5
    K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1   1  -5
    M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1  -4  -5
    B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1  -1  -5
    V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1  -4  -5
    H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1  -1  -5
    D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1  -1  -5
    N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2  -5
    U  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5  -5
    *  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5

.. code:: python

    matrix_from_filename = parasail.Matrix("filename.txt")

SSW Library Emulation
---------------------

`back to top <#table-of-contents>`__

The SSW library (https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) performs Smith-Waterman local alignment using SSE2 instructions and a striped vector.  Its result provides the primary score, a secondary score, beginning and ending locations of the alignment for both the query and reference sequences, as well as a SAM CIGAR.  There are a few parasail functions that emulate this behavior, with the only exception being that parasail does not calculate a secondary score.

.. code:: python

    score_size = 1 # 0, use 8-bit align; 1, use 16-bit; 2, try both
    profile = parasail.ssw_init("asdf", parasail.blosum62, score_size)
    result = parasail.ssw_profile(profile, "asdf", 10, 1)
    print(result.score1)
    print(result.cigar)
    print(result.ref_begin1)
    print(result.ref_end1)
    print(result.read_begin1)
    print(result.read_end1)
    # or skip profile creation
    result = parasail.ssw("asdf", "asdf", 10, 1, parasail.blosum62)

Banded Global Alignment
-----------------------

`back to top <#table-of-contents>`__

There is one version of banded global alignment available.  Though it is not vectorized, it might still be faster than using other parasail global alignment functions, especially for large sequences.  The function signature is similar to the other parasail functions with the only exception being ``k``, the band width.

.. code:: python

    band_size = 3
    result = parasail.nw_banded("asdf", "asdf", 10, 1, band_size, matrix):

File Input
----------

`back to top <#table-of-contents>`__

Parasail can parse FASTA, FASTQ, and gzipped versions of such files if
zlib was found during the C library build. The
function ``parasail.sequences_from_file`` will return a list-like object
containing Sequence instances. A parasail Sequence behaves like an
immutable string but also has extra attributes ``name``, ``comment``,
and ``qual``. These attributes will return an empty string if the input
file did not contain these fields.

Tracebacks
----------

`back to top <#table-of-contents>`__

Parasail supports accessing a SAM CIGAR string from a result.  You must use a traceback-capable alignment function.  Refer to the C interface description above for details on how to use a traceback-capable alignment function.

.. code:: python

    result = parasail.sw_trace("asdf", "asdf", 10, 1, parasail.blosum62)
    cigar = result.cigar
    # cigars have seq, len, beg_query, and beg_ref properties
    # the seq property is encoded
    print(cigar.seq)
    # use decode attribute to return a decoded cigar string
    print(cigar.decode)

Citing parasail
---------------

`back to top <#table-of-contents>`__

If needed, please cite the following paper.

Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global,
and local pairwise sequence alignments. *BMC Bioinformatics*, 17(1),
1-11. doi:10.1186/s12859-016-0930-z

http://dx.doi.org/10.1186/s12859-016-0930-z

License: Battelle BSD-style
---------------------------

`back to top <#table-of-contents>`__

Copyright (c) 2015, Battelle Memorial Institute

1. Battelle Memorial Institute (hereinafter Battelle) hereby grants
   permission to any person or entity lawfully obtaining a copy of this
   software and associated documentation files (hereinafter “the
   Software”) to redistribute and use the Software in source and binary
   forms, with or without modification. Such person or entity may use,
   copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and may permit others to do so, subject to
   the following conditions:

   -  Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimers.

   -  Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

   -  Other than as used herein, neither the name Battelle Memorial
      Institute or Battelle may be used in any form whatsoever without
      the express written consent of Battelle.

   -  Redistributions of the software in any form, and publications
      based on work performed using the software should include the
      following citation as a reference:

   Daily, Jeff. (2016). Parasail: SIMD C library for global,
   semi-global, and local pairwise sequence alignments. *BMC
   Bioinformatics*, 17(1), 1-11. doi:10.1186/s12859-016-0930-z

2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
   OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

