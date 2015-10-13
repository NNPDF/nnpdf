#! /usr/bin/python
import lhapdf, math
import numpy
from matplotlib import pyplot as plt
import scipy.stats

# Which moment
n=3

xmin = 0.1
xmax = 0.99
nstep = 100

xvals = [ xmin + i*((xmax-xmin)/nstep) for i in range(0,nstep+1) ]

#Load replicas
replicas = []
for nrep in range(1,101):
  replicas.append(lhapdf.mkPDF("ClosureTest_MSTW_lvl2_30k_fl", nrep))

#Load replicas
replicas2 = []
for nrep in range(1,101):
  replicas2.append(lhapdf.mkPDF("ClosureTest_MSTW_lvl2", nrep))


#Calculate means
means = []
for x in xvals:
  pdfvals = []
  for rep in replicas:
    pdfvals.append(rep.xfxQ2(21, x, 2))

  means.append(numpy.mean(pdfvals))

#Calculate means2
means2 = []
for x in xvals:
  pdfvals = []
  for rep in replicas2:
    pdfvals.append(rep.xfxQ2(21, x, 2))
  
  means2.append(numpy.mean(pdfvals))


pdfmoment = []
pdfmeans = []
pdferrors = []
pdfup = []
pdfdn = []
for i in xrange(0,len(xvals)):

  pdfvals = []
  for rep in replicas:
    pdfval = rep.xfxQ2(21, xvals[i], 2)
    pdfvals.append(math.pow(pdfval - means[i],n))
  
  pdfmeans.append(numpy.mean(pdfvals))
  pdferrors.append(numpy.std(pdfvals))
  pdfup.append(pdfmeans[-1] + pdferrors[-1])
  pdfdn.append(pdfmeans[-1] - pdferrors[-1])

#line = plt.plot(xvals, pdfmeans, color='g')

line = plt.plot(xvals, pdfup, color='g')
line = plt.plot(xvals, pdfdn, color='g')

print pdfmoment

pdfmeans = []
pdferrors = []
pdfup = []
pdfdn = []
for i in xrange(0,len(xvals)):
  
  pdfvals = []
  for rep in replicas2:
    pdfval = rep.xfxQ2(21, xvals[i], 2)
    pdfvals.append(math.pow(pdfval - means[i],n))
  
  pdfmeans.append(numpy.mean(pdfvals))
  pdferrors.append(numpy.std(pdfvals))
  pdfup.append(pdfmeans[-1] + pdferrors[-1])
  pdfdn.append(pdfmeans[-1] - pdferrors[-1])

#line = plt.plot(xvals, pdfmeans, color='b')

line = plt.plot(xvals, pdfup, color='b')
line = plt.plot(xvals, pdfdn, color='b')

plt.show()
plt.clear()


##########################################

for rep in replicas:
  pdfmom = []
  for i in xrange(0,len(xvals)):
    pdfval = rep.xfxQ2(21, xvals[i], 2)
    pdfmom.append(math.pow(pdfval - means[i],n))

  line = plt.plot(xvals, pdfmom, color='g')

for rep in replicas2:
  pdfmom = []
  for i in xrange(0,len(xvals)):
    pdfval = rep.xfxQ2(21, xvals[i], 2)
    pdfmom.append(math.pow(pdfval - means[i],n))
  
  line = plt.plot(xvals, pdfmom, color='b')


plt.show()

