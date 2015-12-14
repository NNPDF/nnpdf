#!/usr/bin/python
# nnpydf - nh 08/14

import sys, os, gc, shutil

#nnpydf imports
import commondata
from commondata import CommonData

if (len(sys.argv) != 2):
  print " usage: ", sys.argv[0], " [CommonData file]"
  exit(-1)

commonfile = sys.argv[1]
cData = CommonData(commonfile, debug=True)
cData.printData()
