#!/bin/bash
cmake .
make -j${CPU_COUNT}
