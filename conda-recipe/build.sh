#!/bin/bash

# Install n3fit and validphys
python3 -m pip install . --no-deps --ignore-installed

mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX}
make -j${CPU_COUNT}
make install
