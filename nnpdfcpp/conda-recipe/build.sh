#!/bin/bash
cmake . -DCOMPILE_validphys=OFF
make -j${CPU_COUNT}
