#!/usr/bin/python

import sys
import os
import math
import shutil
from subprocess import PIPE, Popen

import matplotlib
from matplotlib import pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv) < 2:
  print "\nusage: ./postfit [configuration file]\n"
  exit()

config  = sys.argv[1]
mask = []
maskout = []


# building the input directory
fitname = config.replace(".ini","")
dir     = "../../results/" + fitname + "/nnfit/"

# check the total number of replicas
im = 1
while True:
  folder = dir + "replica_" + str(im) + "/"
  if os.path.exists(folder) == True:
    mask.append(folder)
  else:
    break

  im += 1

totalrep = len(mask)
nrep = totalrep
print " - There are",totalrep,"replicas."

# check for outliers
im = 1
while True:
  folder = dir + "replica_" + str(im) + ".postfitveto/"
  
  if os.path.exists(folder) == True:
    maskout.append(folder)
  else:
    if im >= 2000:
      break
  
  im += 1

totalout = len(maskout)
nrepout = totalout
print " - There are",totalout,"outliers."

# include outliers
for i in range(0,totalout):
  mask.append(maskout[i])
nrep+= totalout

# check replica's chi2
print " - Building Scatterplot Array:"

npdf = 7
chi2 = []
alpha =[]
beta = []
arcl = []

for ipdf in range(0,npdf):
  alpha.append([])
  beta.append([])
  arcl.append([])

#for ipdf in range(0,npdf):
#  alpha.append([])
#  beta.append([])
#  arcl.append([])


for irep in range(0,nrep):
  # Loop over replica
  file = mask[irep]
  file += fitname + ".fitinfo"
  
  preproc = mask[irep]
  preproc += fitname + ".preproc"
  
  # Loop over files - read chi^2
  replicafile = open(file, 'rb')
  myline = replicafile.readline().split()
  chi2.append(float(myline[3]))
  
  # read arclengths
  arcline = replicafile.readline().split()
  for ipdf in range(0,npdf):
        arcl[ipdf].append(float(arcline[ipdf]))
  
  replicafile.close()
  
  # Read preproc
  preprocfile = open(preproc, 'rb')
  
  for ipdf in range(0,npdf):
    myline = preprocfile.readline().split()
    alpha[ipdf].append(float(myline[0]))
    beta[ipdf].append(float(myline[1]))
  
  preprocfile.close()

outfilename=fitname+".dat"
outpreproc = open(outfilename, 'w')
outpreproc.write("<chi2> <alpha1> <beta1> <arcl1>......  <alphaN> <betaN> <arclN> \n")

for ipdf in range(0,npdf):
  f, (ax1, ax2) = matplotlib.pyplot.subplots(1, 2, sharey=True)
  f.add_subplot(1,2,1)
  plt.ylabel("Chi^2")
  plt.xlabel("Exponent")
  ax1.scatter(alpha[ipdf], chi2)
  ax1.set_title('Alpha Exponent')
  ax2.scatter(beta[ipdf], chi2)
  ax2.set_title('Beta Exponent')
  pltname="1dfig_"+str(ipdf)+".pdf"
  f.savefig(pltname)

for irep in range(0,nrep):
  line = str(chi2[irep]) + " "
  for ipdf in range(0,npdf):
    line += str(alpha[ipdf][irep]) + " " + str(beta[ipdf][irep]) + " " + str(arcl[ipdf][irep]) + " "
  line += "\n"
  outpreproc.write(line)

outpreproc.close()


exit(1)
