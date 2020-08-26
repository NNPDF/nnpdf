"""
CommonData in Python - nh 07/2014
Classes to perform reading/writing of
CommonData files
"""

import os, sys
from math import sqrt

def HERASort( ECoM ):
  """ Sort CoM energy into a valid HERA energy """
  # Correspond to Ep = [460,575,820,920]
  HERACoMEnergies = [225, 252, 301, 319]
  return min(HERACoMEnergies, key=lambda x:abs(x-ECoM))

class CommonData:
  """ Python CommonData class """
  def __init__(self, source, debug = False):
    """ CommonData constructor - takes as argument the filename of a CommonData record """
    # Verify path
    if os.path.exists(source) == False:
      print "CommonData Error: source file \'" + source + "\' not found!"
      sys.exit()

    if (debug == True):
      print "Reading from file: \'"+source+"\'"

    # Open commondata
    datafile = open(source, 'rb')
    metadata = datafile.readline().split()

    self.setname = metadata[0]
    self.nsys = int(metadata[1])
    self.ndata = int(metadata[2])

    self.data = []
    self.kin1 = []
    self.kin2 = []
    self.kin3 = []
    self.stat = []
    self.sys  = []

    self.syserr = []
    self.error  = []

    for line in datafile:

      linesplit = line.split()
      idat = int(linesplit[0])
      self.proc = linesplit[1]

      self.kin1.append(float(linesplit[2]))
      self.kin2.append(float(linesplit[3]))
      self.kin3.append(float(linesplit[4]))

      self.data.append(float(linesplit[5]))
      self.stat.append(float(linesplit[6]))

      self.sys.append([])

      # Replace rapidity with CoM Energy
      if self.proc.startswith("DIS"):
        if self.kin3[-1] != 0:
          self.kin3[-1] = sqrt(self.kin2[-1]/(self.kin1[-1]*self.kin3[-1]))
          if "HERA" in self.setname: # Use HERASORT
            self.kin3[-1] = HERASort(self.kin3[-1])
          elif "NMC" in self.setname: # Agressive round
            self.kin3[-1] = round(self.kin3[-1],0)
          else: # Just round
            self.kin3[-1] = round(self.kin3[-1],0)

      # Need to do proper treatment of correlated errors/systypes etc
      if self.nsys > 0:
        for isys in range(7,self.nsys*2):
          if (isys % 2) != 0:
            self.sys[len(self.sys)-1].append(float(linesplit[isys]))


    # Systematic Error (includes lumi)



    # Total Error
    for idat in range(0,self.ndata):
      error = self.stat[idat]*self.stat[idat]
      for syserror in self.sys[idat]:
        self.syserr.append(sqrt(syserror*syserror))
        error += syserror*syserror
      self.error.append(sqrt(error))

    self.printMetadata()


  def printMetadata(self):
    """ Prints the metadata of a CommonData file """
    print "************COMMONDATA************"
    print "Setname:", self.setname, "PROC:", self.proc
    print "NDAT:", self.ndata,"NSYS:",self.nsys


  def printData(self):
    """ Prints the datapoints present in a CommonData file """
    print "IDAT","KIN1","KIN2","CV", "STAT", "SYS"
    for idat in range(0,self.ndata):
      print idat, self.kin1[idat], self.kin2[idat], self.data[idat], self.stat[idat], self.syserr[idat]
