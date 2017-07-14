#!/usr/bin/env bash

#Find conda
source ~/.bashrc
conda build -q conda-recipe
if [ $? != 0 ]; then
	echo failed to build
	exit 1
fi

# install cmake
conda install cmake git

# build library and test
cmake .
make
make catch_test
