#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Test CommonData """

__authors__ = 'Stefano Carrazza, et at.'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__   = 'stefano.carrazza@cern.ch'

import argparse
import math
from NNPDF import CommonData

#_____________________________________________________
def main(filename, sysfile):

  ## load data from commondata and sysfile
  data = CommonData.ReadFile(filename,sysfile)

  ## read and allocate lists with central value and stat. uncer.
  cvs    = [ data.GetData(el) for el in range(data.GetNData()) ]
  stserr = [ data.GetStat(el) for el in range(data.GetNData()) ]
  syscor = []
  sysunc = []
  syslum = [ 0 for el in range(data.GetNData())]

  syscor2 = [ 0 for el in range(data.GetNData()) ]
  sysunc2 = [ 0 for el in range(data.GetNData()) ]

  for el in range(data.GetNData()):
    for isys in range(data.GetNSys()):
        if (data.GetSys(el,isys).name == "UNCORR"):
          sysunc2[el] += (pow(data.GetSys(el,isys).add,2))
        elif (data.GetSys(el,isys).name == "CORR"):
          syscor2[el] += (pow(data.GetSys(el,isys).add,2))
        elif ("LUMI" in data.GetSys(el,isys).name ):
          syslum[el] = (data.GetSys(el,isys).add)
    sysunc.append(math.sqrt(sysunc2[el]))
    syscor.append(math.sqrt(syscor2[el]))

#      syserr.append(data.GetNSys())
#      print (data.GetSys(el,isys) for isys in range(data.GetNSys()))

  ## print the results
  print (" Id   Central         Stat.           Corr.           Uncorr.          Lumi")
  print ("-----------------------------------------------------------------------------------------------")
  for i in range(data.GetNData()):
    print ("%3i    %5.2f     %5.2f (%5.2f%%)   %5.2f (%4.2f%%)   %5.2f (%4.2f%%)   %5.2f (%4.2f%%)"
            % (i, cvs[i]/1e3,
            stserr[i]/1e3, stserr[i]/cvs[i]*100,
            syscor[i]/1e3, syscor[i]/cvs[i]*100,
            sysunc[i]/1e3, sysunc[i]/cvs[i]*100,
            syslum[i]/1e3, syslum[i]/cvs[i]*100)
            )

  ## produce a plot
  import matplotlib.pyplot as plt
  plt.figure()
  plt.errorbar(range(len(cvs)), cvs, yerr=stserr, fmt='o')
  plt.title('CommonData test with ' + data.GetSetName())
  plt.ylabel('Prediction')
  plt.xlabel('Data Point')
  plt.savefig('%s.pdf' % data.GetSetName())

#_____________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('filename',nargs='?',
                      help='CommonData filename')
  parser.add_argument('sysfile',nargs='?',
                      help='Sys. type filename')
  args = parser.parse_args()
  if not all((args.filename, args.sysfile)):
    parser.error('Too few arguments: filename, sysfile')
  mainargs = vars(args)
  main(**mainargs)
