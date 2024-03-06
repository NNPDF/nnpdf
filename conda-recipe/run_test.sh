#!/bin/bash
set -u
set -v
set -e

#Python tests for the installed validphys package
pytest --pyargs --mpl validphys

platformstr=`uname`

# skip n3fit tests on mac
if [[ "$platformstr" != "Darwin" ]]; then
    pytest --pyargs n3fit
else
    echo "Skipping n3fit tests on Mac"
fi
