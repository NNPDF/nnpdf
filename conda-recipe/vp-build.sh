#!/bin/bash

mkdir build
cd build
echo "build_version=\"${PKG_VERSION}\"" > ../validphys2/src/validphys/version.py

cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX} -DVP_DEV=OFF -DENABLE_N3=OFF -DBURN_TAG=OFF
make -j${CPU_COUNT}
make install
