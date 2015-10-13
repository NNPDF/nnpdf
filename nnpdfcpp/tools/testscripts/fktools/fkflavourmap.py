#!/usr/bin/python

import sys
import os
import math
import shutil
from subprocess import PIPE, Popen


if len(sys.argv) < 2:
  print "\nusage: ./fkflavourmap [FKTable]\n"
  exit()

tolerance = 1E-12

print "*** FK Flavourmap Finder"

fkfilename  = sys.argv[1]

nonzero = []
avgval = []
for ipdf in range(0,14):
  nonzero.append([])
  avgval.append([])
  for jpdf in range(0,14):
    nonzero[ipdf].append(0)
    avgval[ipdf].append(0)

# Read FK
fktable = open(fkfilename, 'rb')
header = fktable.readline() + fktable.readline() + fktable.readline()
print header

# dum
fktable.readline()

myline = fktable.readline().split()
ndata = int(myline[0])
nx = int(myline[1])

fktable.readline()
had = fktable.readline().split()[0]
fktable.readline()

#existing flavourmap/extra stuff
for ifl in range(0,25):
  fktable.readline()
for ifl in range(0,nx+1):
   fktable.readline()

totentries = 0

while True:
  line=fktable.readline().split()
  if not line: break
  totentries += 1
  for ifl in range(0,14):
    for jfl in range(0,14):
      value = float(line[3+ 14*ifl + jfl])
      avgval[ifl][jfl] = max(abs(value),avgval[ifl][jfl])
      if (abs(value) > tolerance):
        nonzero[ifl][jfl] = 1


print "********* Maximum Value Heirarchy (Log) *********"
print " This table indicates for each flavour channel it's maximum"
print " value with respect to the largest channel (log scale). This"
print " table will show if there are some channels with negligible"
print " contribution. X denotes channels with zero entries."
print " "
maxlist = []
for ifl in range(0,14):
  maxlist.append(max(avgval[ifl]))
maxval = max(maxlist)

for ifl in range(0,14):
  linestring = ""
  for jfl in range(0,14):
    if (avgval[ifl][jfl] > 0):
      calcval = math.log(avgval[ifl][jfl]/maxval)
      linestring = linestring + str(int(calcval))+" "
    else:
      linestring = linestring + " X "
  print linestring

print ""
print " -- Suggested tolerance: 1e", int(math.log(maxval) - 15)
print " -- Current tolerance: ", tolerance
print " -- If tolerances differ, consider editing the tolerance variable to the suggested value"
print " ********* Flavourmap after tolerance *********"
nfl=0
for ifl in range(0,14):
  linestring = ""
  for jfl in range(0,14):
    linestring = linestring + str(nonzero[ifl][jfl])+" "
    if (nonzero[ifl][jfl] == 1):
      nfl = nfl+1
  print linestring


fktable.close()

print "Number of active flavours: ", nfl
exit(1)
