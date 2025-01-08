#!/bin/bash
set -u
set -v
set -e

# Python tests for the installed validphys package
# Note that the default tolerance in the conda test is higher than the pip test
pytest --pyargs --mpl validphys --mpl-default-tolerance 24

platformstr=`uname`

# skip n3fit tests on mac
if [[ "$platformstr" != "Darwin" ]]; then
    pytest --pyargs n3fit
else
    echo "Skipping n3fit tests on Mac"
fi
