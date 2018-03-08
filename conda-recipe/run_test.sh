#!/bin/bash
set -u
set -v
set -e

#Python tests for the installed validphys package
pytest --pyargs validphys


#Compile in debug mode and run the cpp tests


mkdir bldtest
cd bldtest
cmake .. -DENABLE_TESTS=ON
make catch_test -j
./libnnpdf/tests/catch_test


#Print linkage data
conda inspect linkages -p $PREFIX $PKG_NAME
