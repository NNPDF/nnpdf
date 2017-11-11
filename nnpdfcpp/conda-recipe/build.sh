#!/bin/bash
cmake . -DCOMPILE_validphys=OFF -DCOMPILE_mkthpredictions=ON
make -j${CPU_COUNT}
mv ./bin/* ${PREFIX}/bin

cd tools/postfit2
python setup.py install --single-version-externally-managed --record record.txt
