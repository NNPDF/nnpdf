#!/bin/bash
# TESTING DELETEME
export LDFLAGS="-Wl,-pie -Wl,-headerpad_max_install_names -Wl,-rpath,${PREFIX}/lib"

mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX} -DVP_DEV=OFF
make -j${CPU_COUNT}
make install
