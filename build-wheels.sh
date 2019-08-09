#!/bin/bash
set -e -x

# use python 3.6, the resultant wheel is a py2-py3 wheell
PYBIN=/opt/python/cp36-cp36m/bin
"${PYBIN}/pip" wheel /io/ -w wheelhouse/

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
done

# Install wheel into all pythons as a check
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install numpy
    "${PYBIN}/pip" install parasail --no-index -f /io/wheelhouse
    seqa="'ACGTTTTTGCA'"
    seqb="'ACGTTTTGCA'"
    "${PYBIN}/python" -c "import parasail; print(parasail.sw_trace_striped_16(${seqa}, ${seqb}, 8, 4, parasail.dnafull).cigar.decode.decode())"
done
