#!/usr/bin/env bash

#Find conda
source ~/.bashrc

# install cmake
sudo apt get install cmake libgsl-dev sqlite3 libyaml-cpp-dev

# build library and test
cmake .
make
make catch_test
