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

#Check that filter an chi2check run
vp-setupfit --legacy ../nnpdfcpp/config/testconfig.yml
chi2check testconfig NNPDF31_nnlo_as_0118


