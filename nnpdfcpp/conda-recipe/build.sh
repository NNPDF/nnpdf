#!/bin/bash
cmake . -DCOMPILE_validphys=OFF -DCOMPILE_mkthpredictions=ON -DCMAKE_INSTALL_PREFIX=${PREFIX}
make -j${CPU_COUNT}
make install

cd tools/postfit2
python setup.py install --single-version-externally-managed --record record.txt
