#!/usr/bin/python

import sys
import os
import math
import shutil
from subprocess import PIPE, Popen


if len(sys.argv) < 2:
  print "\nusage: ./fktotal [FKTable]\n"
  exit()

# Read FK table
fkfilename  = sys.argv[1]
fktable = open(fkfilename, 'rb')

header = fktable.readline() + fktable.readline() + fktable.readline()
sys.stdout.write(header)

# dum
sys.stdout.write(fktable.readline())

# nx ndata
line = fktable.readline()
myline = line.split()
ndata = int(myline[0])
nx = int(myline[1])
sys.stdout.write(line)

# hadronic
sys.stdout.write(fktable.readline())
line = fktable.readline()
had = line.split()[0]
sys.stdout.write(line)
sys.stdout.write(fktable.readline())

if (had == 0):
  print "Error! Table must be hadronic!"
  exit(-1)

#flavourmap
for ifl in range(0,25):
  sys.stdout.write(fktable.readline())
for ifl in range(0,nx+1):
  sys.stdout.write(fktable.readline())

# Build total xsec array
  sigma_tot = []
  for ifl in range(0,14):
    sigma_tot.append([])
    for jfl in range(0,14):
      sigma_tot[ifl].append([])
      for ix in range(0,nx):
        sigma_tot[ifl][jfl].append([])
        for jx in range(0,nx):
          sigma_tot[ifl][jfl][ix].append(0.0)

for line in iter(lambda: fktable.readline(), ''):
  linesplit = line.split()
  for ifl in range(0,14):
    for jfl in range(0,14):
      sigma_tot[ifl][jfl][int(linesplit[1])][int(linesplit[2])] += float(linesplit[3+ 14*ifl + jfl])


for ix in range(0,nx):
  for jx in range(0,nx):
    nonzero = 0
    fline = ""
    for ifl in range(0,14):
      for jfl in range(0,14):
        fline += "\t "+str(sigma_tot[ifl][jfl][ix][jx])
        if (abs(sigma_tot[ifl][jfl][ix][jx]) > 1E-100):
          nonzero = 1
    if (nonzero == 1):
      print 0,ix,jx,fline

fktable.close()
exit(1)
