#!/bin/bash
export LDFLAGS=$(echo $LDFLAGS | sed 's/-Wl,-dead_strip_dylibs//g')
RED='\033[0;31m'
NC='\033[0m' # No Color
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX} -DVP_DEV=OFF
echo  -e "${RED} ${LDFLAGS} ${NC}"
make -j${CPU_COUNT}
make install
