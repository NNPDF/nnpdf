#!/usr/bin/python

import sys
import os
import math
import shutil
from subprocess import PIPE, Popen


if len(sys.argv) < 2:
  print "\nusage: ./fksum [FKTable1] [FKTable2]\n"
  exit()

tolerance = 1E-12

print "*** FK Sum Script"

fkfilename1  = sys.argv[1]
fkfilename2  = sys.argv[1]

# Read FK tables
fktable1 = open(fkfilename1, 'rb')
fktable2 = open(fkfilename2, 'rb')

header = fktable1.readline() + fktable1.readline() + fktable1.readline()
print header

header = fktable2.readline() + fktable2.readline() + fktable2.readline()
print header

# dum
fktable1.readline()
fktable2.readline()

myline = fktable1.readline().split()
ndata = int(myline[0])
nx = int(myline[1])

myline = fktable2.readline().split()
ndata2 = int(myline[0])
nx2 = int(myline[1])

if (ndata != ndata2):
  print "Error! Number of datapoints do not match!"
  exit(-1)

if (nx != nx2):
  print "Error! Number of datapoints do not match!"
  exit(-1)

fktable1.readline()
had1 = fktable1.readline().split()[0]
fktable1.readline()

fktable2.readline()
had2 = fktable2.readline().split()[0]
fktable2.readline()

if (had1 == 0):
  print "Error! Table must be hadronic!"
  exit(-1)

#existing flavourmap/extra stuff
for ifl in range(0,25):
  fktable1.readline()
for ifl in range(0,nx+1):
  fktable1.readline()

#existing flavourmap/extra stuff
for ifl in range(0,25):
  fktable2.readline()
for ifl in range(0,nx+1):
  fktable2.readline()

# Build arrays
  sigma1 = []
  sigma2 = []
  for idat in range(0,ndata):
    sigma1[idat] = []
    sigma2[idat] = []
    for ifl in range(0,14):
      sigma1[idat][ifl] = []
      sigma2[idat][ifl] = []
      for jfl in range(0,14):
        sigma1[idat][ifl][jfl] = []
        sigma2[idat][ifl][jfl] = []
        for ix in range(0,nx):
          sigma1[idat][ifl][jfl][ix] = []
          sigma2[idat][ifl][jfl][ix] = []

fktable1.close()
fktable2.close()

print "Number of active flavours: ", nfl
exit(1)
