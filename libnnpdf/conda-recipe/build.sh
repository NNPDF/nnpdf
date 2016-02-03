#!/bin/bash
autoreconf -i
./configure --prefix=$PREFIX --enable-pywrap
#This doesn't wor properly with the debian autotools. 
#Just use one core
#make -j${CPU_COUNT}
make
make install

make wrapper
