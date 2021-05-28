#!/bin/bash

echo "build_version=\"${PKG_VERSION}\"" > n3fit/src/n3fit/version.py

python -m pip install --no-deps --ignore-installed n3fit