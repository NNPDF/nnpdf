#!/bin/bash
set -u
set -v
set -e

#Aparently this is undefined in Mac
export CPPFLAGS=${DEBUG_CPPFLAGS:-}
export CXXFLAGS=$DEBUG_CXXFLAGS
export CFLAGS=$DEBUG_CFLAGS

mkdir bldtest
cd bldtest
cmake .. -DCMAKE_BUILD_TYPE=debug -DENABLE_TESTS=ON
make catch_test -j
./libnnpdf/tests/catch_test
