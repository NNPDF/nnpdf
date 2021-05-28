#!/bin/bash
set -u
set -v
set -e

pytest --pyargs n3fit

# check filter works
vp-setupfit --legacy nnpdfcpp/config/testconfig.yml
#Check that chi2check run
chi2check testconfig NNPDF31_nnlo_as_0118