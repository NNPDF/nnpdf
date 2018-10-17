#!/bin/bash
#export LDFLAGS=$(echo $LDFLAGS | sed 's/-Wl,-dead_strip_dylibs//g')

mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX} -DVP_DEV=OFF
make -j${CPU_COUNT}
make install
