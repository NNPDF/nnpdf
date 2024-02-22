#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Test CommonData """

__authors__ = 'Stefano Carrazza, et at.'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@cern.ch'

import argparse
import math

from NNPDF import CommonData


# _____________________________________________________
def main(setname):
    sysfile = "systypes/SYSTYPE_" + setname + "_0.dat"
    setfile = "commondata/DATA_" + setname + ".dat"

    ## load data from commondata and sysfile
    data = CommonData.ReadFile(setfile, sysfile)

    ## read and allocate lists with central value and stat. uncer.
    cvs = [data.GetData(el) for el in range(data.GetNData())]
    stserr = [data.GetStat(el) for el in range(data.GetNData())]
    syscor = []
    sysunc = []
    systot = []
    sysbeam = [0 for el in range(data.GetNData())]
    syslum = [0 for el in range(data.GetNData())]

    syscor2 = [0 for el in range(data.GetNData())]
    sysunc2 = [0 for el in range(data.GetNData())]

    for el in range(data.GetNData()):
        for isys in range(data.GetNSys()):
            if data.GetSys(el, isys).name == "UNCORR":
                sysunc2[el] += pow(data.GetSys(el, isys).add, 2)
            elif "LUMI" in data.GetSys(el, isys).name:
                syslum[el] = data.GetSys(el, isys).add
            elif "BEAM" in data.GetSys(el, isys).name:
                sysbeam[el] = data.GetSys(el, isys).add
            else:
                syscor2[el] += pow(data.GetSys(el, isys).add, 2)

        sysunc.append(math.sqrt(sysunc2[el]))
        syscor.append(math.sqrt(syscor2[el]))

        systot.append(math.sqrt(sysunc2[el] + syscor2[el] + pow(syslum[el], 2)))

    ## print results
    print(" Id      CV     Stat.   CorSys  UncSys  Beam       Lumi        Tot. Unc.")
    print("--------------------------------------------------------------------------")
    for i in range(data.GetNData()):
        #    print ("%3i    %6.2f     %6.2f (%6.2f%%)   %6.2f (%6.2f%%)   %6.2f (%6.2f%%)   %6.2f (%6.2f%%)"
        #            % (i, cvs[i]*100,
        #            stserr[i]*100, stserr[i]/cvs[i]*100,
        #            syscor[i]*100, syscor[i]/cvs[i]*100,
        #            sysunc[i]*100, sysunc[i]/cvs[i]*100,
        #            syslum[i]*100, syslum[i]/cvs[i]*100)
        #            )
        print(
            "%3i   %7.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f(%4.2f%%)  %6.3f(%4.2f%%)"
            % (
                i,
                cvs[i] * 0.001,
                stserr[i] * 0.001,  # stserr[i]/cvs[i]*100,
                syscor[i] * 0.001,  # syscor[i]/cvs[i]*100,
                sysunc[i] * 0.001,  # sysunc[i]/cvs[i]*100,
                sysbeam[i] * 0.001,  # sysbeam[i]/cvs[i]*100,
                syslum[i] * 0.001,
                syslum[i] / cvs[i] * 100,
                systot[i] * 0.001,
                abs(systot[i] / cvs[i]) * 100,
            )
        )


## produce a plot
#  import matplotlib.pyplot as plt
#  plt.figure()
#  plt.errorbar(range(len(cvs)), cvs, yerr=stserr, fmt='o')
#  plt.title('CommonData test with ' + data.GetSetName())
#  plt.ylabel('Prediction')
#  plt.xlabel('Data Point')
#  plt.savefig('%s.pdf' % data.GetSetName())

# _____________________________________________________
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('setname', nargs='?', help='CommonData setname')
    #  parser.add_argument('sysfile',nargs='?',
    #                      help='Sys. type filename')
    args = parser.parse_args()
    if not all((args.setname)):
        parser.error('Too few arguments: setname')
    mainargs = vars(args)
    main(**mainargs)
