==========
Change Log
==========

All notable changes to this project will be documented in this file.
This project adheres to Semantic Versioning (http://semver.org/).

-----------
Unreleased_
-----------
The Unreleased section will be empty for tagged releases. Unreleased functionality appears in the develop branch.

-------------------
1.3.4_ - 2023-02-16
-------------------
- Version bump to regenerate wheels with v2.6.1 parasail C library

-------------------
1.3.3_ - 2022-08-26
-------------------
- CI improvements
  - build all wheel targets in parallel 
  - run tests as part of build
- fix win_amd64 wheel builds

-------------------
1.3.2_ - 2022-08-25
-------------------
- Version bump to force rebuilding of wheels

-------------------
1.3.1_ - 2022-08-24
-------------------
- Fix macos wheel build
- Add C parasail headers into wheels, add get_include() and get_library()

-----------------
1.3_ - 2022-08-22
-----------------
- Add bindings for sg_qb_db and sg_qe_de

-------------------
1.2.4_ - 2021-02-06
-------------------
- Fix #60 Missing wheels on PyPI

-------------------
1.2.3_ - 2020-12-16
-------------------
- Fix #57 Install error with parasail v1.2.2

-------------------
1.2.2_ - 2020-12-01
-------------------
- Fix #55 Error installing parasail on macos

-------------------
1.2.1_ - 2020-11-03
-------------------
- Fix #54 bug in ssw_init()

-----------------
1.2_ - 2020-02-26
-----------------
- Automated testsuite.
- Fix wheel deployment.
- Update email.

--------------------
1.1.20_ - 2020-02-25
--------------------
- Build wheels with cibuildwheel.
- Fix #11 manylinux1 builds for pypi.

--------------------
1.1.19_ - 2019-11-12
--------------------
- Missed version bump during last release.

--------------------
1.1.18_ - 2019-11-12
--------------------
- Allow install without downloading and building C library, set evn var PARASAIL_SKIP_BUILD=1 to enable.

--------------------
1.1.17_ - 2019-04-29
--------------------
- Add wheel package as dependency.
- Bump version to create new pypi release with latest 2.4.1 C library.

--------------------
1.1.16_ - 2019-02-14
--------------------
- AppVeyor now uses python35 to fix packaging issue. Fixes #33.

--------------------
1.1.15_ - 2019-02-12
--------------------
- Fix travis CI pypi deployment, again.

--------------------
1.1.14_ - 2019-02-11
--------------------
- Fix travis CI pypi deployment.

--------------------
1.1.13_ - 2019-02-11
--------------------
- Bindings up to date with 2.4 C library.
- New semi-global routines.
- Case-sensitive support.
- Users may optionally set string encoding using parasail.set_encoding(name). Default remains 'latin_1'.

--------------------
1.1.12_ - 2018-07-12
--------------------
- Fix autotools build by updating automake version.
- Add functions from latest 2.2 C library.

--------------------
1.1.11_ - 2018-03-06
--------------------
- Fix Cigar lifetime bug. Cigar instance now held by Result to avoid freeing C memory early.
- Fix poor C library performance on OSX by disabling multiarch build.

--------------------
1.1.10_ - 2018-01-16
--------------------
- Bump version to create new pypi release with latest 2.1 C library.

-------------------
1.1.9_ - 2018-01-16
-------------------
- Fix compatibility with 2.1.x C parasail release series
  - configure script, if missing,  must be generated using autoreconf -fi
  - download autotools if missing
- Fix another Python3 issue with Cigar. Fixes #16 again.

-------------------
1.1.8_ - 2018-01-09
-------------------
- Fix relative import statements. Fixes #13.
- Fix Python3 issues with Profile. Fixes #14.
- Fix Python3 issues with Cigar. Fixes #15, #16.

-------------------
1.1.7_ - 2017-10-18
-------------------
- Bump version to create new pypi release with latest 2.0.2 C library.

-------------------
1.1.6_ - 2017-10-03
-------------------
- Fix #10.

-------------------
1.1.5_ - 2017-09-29
-------------------
- Add travis build.
- Add CHANGELOG.rst.
- Add Cigar decode_op() decode_len() @staticmethods
- Add bindings to the ssw emulation functions.
- Add bindings to parasail_sequences_from_file(filename).
- Add tests; use with pytest.
- Fix the integer type for ssw cigar.
- Fix the profile function return values.
- Fix the table and rowcol properties.

-------------------
1.1.4_ - 2017-09-26
-------------------
- Py 2/3 compatible long_description field.

-------------------
1.1.3_ - 2017-09-26
-------------------
- Assign README.rst to long_description field.

-------------------
1.1.2_ - 2017-09-26
-------------------
- Create a new pypi release that didn't fail like the last one did.

-------------------
1.1.1_ - 2017-09-26
-------------------
- Expanded README.rst.
- pypi upload using twine now skips existing files.

-------------------
1.1.0_ - 2017-09-26
-------------------
First tagged release. The 'master' branch always represents the latest stable code. Tagged releases correspond to pypi releases.

.. _Unreleased: https://github.com/jeffdaily/parasail-python/compare/v1.3.4...master
.. _1.3.4:   https://github.com/jeffdaily/parasail-python/compare/v1.3.3...v1.3.4
.. _1.3.3:   https://github.com/jeffdaily/parasail-python/compare/v1.3.2...v1.3.3
.. _1.3.2:   https://github.com/jeffdaily/parasail-python/compare/v1.3.1...v1.3.2
.. _1.3.1:   https://github.com/jeffdaily/parasail-python/compare/v1.3...v1.3.1
.. _1.3:     https://github.com/jeffdaily/parasail-python/compare/v1.2.4...v1.3
.. _1.2.4:   https://github.com/jeffdaily/parasail-python/compare/v1.2.3...v1.2.4
.. _1.2.3:   https://github.com/jeffdaily/parasail-python/compare/v1.2.2...v1.2.3
.. _1.2.2:   https://github.com/jeffdaily/parasail-python/compare/v1.2.1...v1.2.2
.. _1.2.1:   https://github.com/jeffdaily/parasail-python/compare/v1.2...v1.2.1
.. _1.2:     https://github.com/jeffdaily/parasail-python/compare/v1.1.20...v1.2
.. _1.1.20:  https://github.com/jeffdaily/parasail-python/compare/v1.1.19...v1.1.20
.. _1.1.19:  https://github.com/jeffdaily/parasail-python/compare/v1.1.18...v1.1.19
.. _1.1.18:  https://github.com/jeffdaily/parasail-python/compare/v1.1.17...v1.1.18
.. _1.1.17:  https://github.com/jeffdaily/parasail-python/compare/v1.1.16...v1.1.17
.. _1.1.16:  https://github.com/jeffdaily/parasail-python/compare/v1.1.15...v1.1.16
.. _1.1.15:  https://github.com/jeffdaily/parasail-python/compare/v1.1.14...v1.1.15
.. _1.1.14:  https://github.com/jeffdaily/parasail-python/compare/v1.1.13...v1.1.14
.. _1.1.13:  https://github.com/jeffdaily/parasail-python/compare/v1.1.12...v1.1.13
.. _1.1.12:  https://github.com/jeffdaily/parasail-python/compare/v1.1.11...v1.1.12
.. _1.1.11:  https://github.com/jeffdaily/parasail-python/compare/v1.1.10...v1.1.11
.. _1.1.10:  https://github.com/jeffdaily/parasail-python/compare/v1.1.9...v1.1.10
.. _1.1.9:   https://github.com/jeffdaily/parasail-python/compare/v1.1.8...v1.1.9
.. _1.1.8:   https://github.com/jeffdaily/parasail-python/compare/v1.1.7...v1.1.8
.. _1.1.7:   https://github.com/jeffdaily/parasail-python/compare/v1.1.6...v1.1.7
.. _1.1.6:   https://github.com/jeffdaily/parasail-python/compare/v1.1.5...v1.1.6
.. _1.1.5:   https://github.com/jeffdaily/parasail-python/compare/v1.1.4...v1.1.5
.. _1.1.4:   https://github.com/jeffdaily/parasail-python/compare/v1.1.3...v1.1.4
.. _1.1.3:   https://github.com/jeffdaily/parasail-python/compare/v1.1.2...v1.1.3
.. _1.1.2:   https://github.com/jeffdaily/parasail-python/compare/v1.1.1...v1.1.2
.. _1.1.1:   https://github.com/jeffdaily/parasail-python/compare/v1.1.0...v1.1.1
.. _1.1.0:   https://github.com/jeffdaily/parasail-python/releases/tag/v1.1.0

