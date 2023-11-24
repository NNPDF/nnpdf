#!/bin/bash

# Install n3fit and validphys
python -m pip install . --no-deps --ignore-installed

mkdir build
cd build
echo "build_version=\"${PKG_VERSION}\"" > ../n3fit/src/n3fit/version.py
echo "build_version=\"${PKG_VERSION}\"" > ../validphys2/src/validphys/version.py

cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX}
make -j${CPU_COUNT}
make install
