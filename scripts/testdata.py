#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Test CommonData """

__authors__ = 'Stefano Carrazza, et at.'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@cern.ch'

import argparse
from NNPDF import CommonData

#_____________________________________________________
def main(filename, sysfile):

  ## load data from commondata and sysfile
  data = CommonData.ReadFile(filename,sysfile)

  ## read and allocate lists with central value and stat. uncer.
  cvs = [ data.GetData(el) for el in range(data.GetNData()) ]
  sts = [ data.GetStat(el) for el in range(data.GetNData()) ]

  ## print the results
  for i in range(data.GetNData()):
    print ("CV=%5e STAT=%.5e" % (cvs[i],sts[i]))    

  ## produce a plot
  import matplotlib.pyplot as plt
  plt.figure()
  plt.errorbar(range(len(cvs)), cvs, yerr=sts, fmt='o')
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
