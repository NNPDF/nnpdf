#!/usr/bin/python
# nnpydf - nh 08/14

import sys, os, gc, shutil

#nnpydf imports
import thpredictions, commondata, datplot, dat_md
from thpredictions import ThPredictions
from commondata import CommonData
from datplot import *

if (len(sys.argv) < 2):
  print " usage: ", sys.argv[0], " [CommonData file] "
  exit(-1)

#os.mkdir("plots")

for ifile in sys.argv[1:]:
  cData = CommonData(ifile, debug=True)
  dataPlot = dataplot_nokin(cData)
  dataPlot.export("plots")