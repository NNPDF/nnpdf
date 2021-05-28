#!/bin/bash
set -u
set -v
set -e

pytest --pyargs n3fit

# check filter works
vp-setupfit --legacy nnpdfcpp/config/testconfig.yml
