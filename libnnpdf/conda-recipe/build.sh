#!/bin/bash
mkdir build
cd build
cmake ../ -DENABLE_PYWRAP=on -DCMAKE_INSTALL_PREFIX=${PREFIX} -DAVX_FOUND=OFF -DAVX2_FOUND=OFF
make
make install

make wrapper
