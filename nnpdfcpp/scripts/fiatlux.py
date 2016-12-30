#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
""" A postfit script which collects fiatlux replicas in a LHAPDF grid"""

__authors__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@cern.ch'

import os
import shutil
import sys
import argparse
import math
from subprocess import PIPE, Popen


def main(nrep, fit_filename):

    fitname = fit_filename.replace(".yml","")
    dir     = "../results/" + fitname + "/fiatlux/"
    fitname_lux = fitname + '_fiatlux'

    shutil.rmtree(dir + fitname_lux,  ignore_errors=True)
    os.makedirs(dir + fitname_lux)

    # header
    lhapath = Popen(["lhapdf-config","--datadir"],stdout=PIPE).communicate()[0]
    lhapath = lhapath.decode()
    lhapath = lhapath.replace('\n','/')
    with open(lhapath + fitname + '/' + fitname + '.info', 'r') as header:
        oheader = open(dir + '/' + fitname_lux + '/' + fitname_lux + '.info', 'w')
        for line in header.readlines():
            if 'NumMembers' in line:
                oheader.write('NumMembers: %d\n' % nrep)
            elif 'Flavors' in line:
                if not '22' in line:
                    oheader.write(line.replace(']', ', 22]'))
                else:
                    oheader.write(line)
            else:
                oheader.write(line)
        oheader.close()

    ## Preparing replicas
    xpdf = []
    xgrid = []
    qgrid = []
    fgrid = []
    # copy replica files
    for i in range(1,nrep+1):
        replica = dir + 'replica_' + str(i) + ".dat"
        shutil.copyfile(replica, dir + "/" + fitname_lux + "/" + fitname_lux + '_{:04n}.dat'.format(i))

        print("Reading: %s" % replica)
        f = open(replica, 'r')
        xpdf.append([])
        for j in range(0,2): f.readline()

        s = 0
        while True:
            f.readline()
            xs = f.readline()
            qs = f.readline()
            fs = f.readline()

            nx  = len(xs.split())
            nq  = len(qs.split())
            nfl = len(fs.split())

            if nx == 0: break

            xpdf[i-1].append([])

            if i == 1:
                xgrid.append(xs)
                qgrid.append(qs)
                fgrid.append(fs)

            for ix in range(nx):
                xpdf[i-1][s].append([])
                for iq in range(nq):
                    xpdf[i-1][s][ix].append([])
                    line = f.readline().split()
                    for ifl in range(nfl):
                        xpdf[i-1][s][ix][iq].append(float(line[ifl]))
            s+=1
        f.close()

    print("Computing and priting replica 0")
    f = open(dir + "/" + fitname_lux + "/" + fitname_lux + "_0000.dat", 'w')
    f.write("PdfType: central\n")
    f.write("Format: lhagrid1\n---\n")

    for s in range(len(qgrid)):
        f.write(xgrid[s])
        f.write(qgrid[s])
        f.write(fgrid[s])
        for ix in range(len(xgrid[s].split())):
            for iq in range(len(qgrid[s].split())):
                f.write(" ")
                for ifl in range(len(fgrid[s].split())):
                    sum = 0
                    for irep in range(nrep):
                        sum += xpdf[irep][s][ix][iq][ifl]
                    sum /= nrep
                    f.write("%14.7E " % sum)
                f.write("\n")
        f.write("---\n")
    f.close()

    print("\n- Finished see: \n%s" % dir + fitname_lux)

    print("\n- Copying grid to LHAPDF path.")
    src = dir + fitname_lux
    dst = lhapath + fitname_lux
    try:
        shutil.copytree(src,dst)
    except:
        print("Error: this grid already exists, please delete and run the script again")
        exit(-1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('nrep', nargs='?', help="Number of desired replicas", type=int)
    parser.add_argument('fit_filename', nargs='?', help="Fit configuration filename")
    args = parser.parse_args()
    if not all((args.nrep, args.fit_filename)):
        parser.error("Too few arguments: nrep, fit_filename.")
    mainargs = vars(args)
    main(**mainargs)
