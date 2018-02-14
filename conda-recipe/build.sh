#!/bin/bash
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX} -DVP_DEV=0FF
make -j${CPU_COUNT}
make install
