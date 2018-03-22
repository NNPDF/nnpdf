#!/bin/bash
set -u
set -v
set -e

#Python tests for the installed validphys package
pytest --pyargs validphys


export LDFLAGS=$(echo $LDFLAGS | sed 's/-Wl,-dead_strip_dylibs//g')

mkdir bldtest
cd bldtest
cmake .. -DENABLE_TESTS=ON
make catch_test -j
./libnnpdf/tests/catch_test


#Print linkage data
conda inspect linkages -p $PREFIX $PKG_NAME
