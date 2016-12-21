#!/bin/bash
mkdir build
cd build
cmake ../ -DENABLE_PYWRAP=on
make
make install

make wrapper
