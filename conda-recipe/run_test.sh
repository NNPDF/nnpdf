#!/bin/bash
set -u
set -v
set -e

#Python tests for the installed validphys package
pytest --pyargs --mpl validphys

platformstr=`uname`

# skip n3fit tests on mac
if [[ "$platformstr" != "Darwin" ]]; then
    pytest --pyargs n3fit
else
    echo "Skipping tests on Mac"
fi

mkdir bldtest
cd bldtest
cmake .. -DENABLE_TESTS=ON -DBURN_TAG=OFF
make catch_test -j
./libnnpdf/tests/catch_test
