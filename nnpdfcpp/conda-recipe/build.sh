#!/bin/bash
cmake . -DCOMPILE_validphys=OFF -DCOMPILE_mkthpredictions=ON
make -j${CPU_COUNT}
make install

cd tools/postfit2
python setup.py install --single-version-externally-managed --record record.txt
