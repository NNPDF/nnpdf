#!/bin/bash
mkdir build
cd build
cmake ../ -DENABLE_PYWRAP=on -DCMAKE_INSTALL_PREFIX=${PREFIX}
make
make install

make wrapper
