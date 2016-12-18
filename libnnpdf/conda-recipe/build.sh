#!/bin/bash
./autogen.sh
USE_REAL=${NNPDF_USE_DOUBLE+"--enable-safemode"}
./configure --prefix=$PREFIX --enable-pywrap ${USE_REAL}
#This doesn't wor properly with the debian autotools. 
#Just use one core
#make -j${CPU_COUNT}
make
make install

make wrapper
