#!/bin/bash
set -u
set -v
set -e

#Download some resources.
#This is here to help isolate the errors and avoid timeouts
vp-get theoryID 162
vp-get pdf NNPDF31_nnlo_as_0118

#Python tests for the installed validphys package
pytest --pyargs --mpl validphys
pytest --pyargs n3fit

mkdir bldtest
cd bldtest
cmake .. -DENABLE_TESTS=ON
make catch_test -j
./libnnpdf/tests/catch_test

#Check that filter an chi2check run
vp-setupfit ../nnpdfcpp/config/testconfig.yml
chi2check testconfig NNPDF31_nnlo_as_0118


#Print linkage data
conda inspect linkages -p $PREFIX $PKG_NAME
