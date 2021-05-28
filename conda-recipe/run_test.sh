#!/bin/bash
set -u
set -v
set -e

#Python tests for the installed validphys package
pytest --pyargs --mpl validphys

platformstr=`uname`

mkdir bldtest
cd bldtest
cmake .. -DENABLE_TESTS=ON -DBURN_TAG=OFF
make catch_test -j
./libnnpdf/tests/catch_test

#Check that chi2check run
chi2check testconfig NNPDF31_nnlo_as_0118


#Print linkage data
conda inspect linkages -p $PREFIX $PKG_NAME
