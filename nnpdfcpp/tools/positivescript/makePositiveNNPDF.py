#!/usr/bin/python

"""
usage: ./MakePositiveNNPDF.py [LHgrid file]

This script loops inside a NNPDF20int LHgrid file, 
replaces negative numbers for replicas [1,100] with 
1e-8 and recomputes replica 0.
"""

import sys

# grid setup
eps    = 1e-08
epsSbar= 1e-06
nrep   = 100
nx     = 100
nq     = 50
nfl    = 13
nlines = 282 # header lines

# Opening files
filename = sys.argv[1]
infile = open(filename,'rb')
oufile = open(filename.replace(".LHgrid","_mc.LHgrid"), 'wb')

# Skipping header
print " - Reading header..."
for l in xrange(0,nlines):
    text = infile.readline()
    text = text.replace( ".LHgrid", "_mc.LHgrid" )
    text = text.replace( "NLO global fit", "NLO global fit for MCs")
    text = text.replace( "mem=0 --> average on replicas.", 
                         "mem=0 --> average on replicas, positive definite")
    text = text.replace( "Evolution: pol. interpolation on the LH grid.", 
                         "mem=1,100 --> positive definite above q = 10 GeV")
    oufile.write(text)
    

# Skipping replica 0
print " - Skipping replica 0..."
for ix in range(0,nx):
    for iq in range(0,nq):
        infile.readline()

# Read replicas 1 to 100 and replace negative values with
array = []
for irep in range(0,nrep):
    negrep = 0
    print " - Reading replica", irep+1,
    array.append([])
    for ix in range(0,nx):
        array[irep].append([])
        for iq in range(0,nq):
            array[irep][ix].append([])
            line = infile.readline().split()
            for ifl in range(0, nfl):
                array[irep][ix][iq].append(float(line[ifl]))		
                if float(line[ifl]) < 0:
                    negrep+=1
                    if ifl == 3:
                        array[irep][ix][iq][ifl] = epsSbar # sbar
                    else:
                        array[irep][ix][iq][ifl] = eps
    print "found", negrep, "negative numbers"
                    
# Compute and print AVG, new replica 0
print " - Computing replica 0"
for ix in range(0,nx):
    for iq in range(0,nq):
        oufile.write(" ")
        for ifl in range(0, nfl):
            sum = 0
            for irep in range(0,nrep):
                sum += array[irep][ix][iq][ifl]
            sum /= nrep
            print >> oufile, "%14.7E" % sum,
        oufile.write("\n")

# Printing replicas 1 to 100
for irep in range(0,nrep):
    print " - Printing replica", irep+1
    for ix in range(0,nx):
        for iq in range(0,nq):
            oufile.write(" ")
            for ifl in range(0, nfl):
                print >> oufile, "%14.7E" % array[irep][ix][iq][ifl],
            oufile.write("\n")

infile.close()
print >> oufile, "  'End:'"
oufile.close()

print " - Finished see", oufile
