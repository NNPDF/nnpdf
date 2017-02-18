#!/bin/bash
cmake . -DCOMPILE_validphys=OFF -DCOMPILE_mkthpredictions=ON
make -j${CPU_COUNT}
mv ./bin/* ${PREFIX}/bin
