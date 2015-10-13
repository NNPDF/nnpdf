# Commondata in Python
# nh 07/2014
import os, sys
from math import sqrt


# Sort CoM energy into a valid HERA energy
def HERASort( ECoM ):
  # Correspond to Ep = [460,575,820,920]
  HERACoMEnergies = [225, 252, 301, 319]
  return min(HERACoMEnergies, key=lambda x:abs(x-ECoM))

class CommonData:
  def __init__(self, source):
    
    # Verify path
    if os.path.exists(source) == False:
      print "Error: source file" + source + " not found!"
      sys.exit()

    # Open commondata
    datafile = open(source, 'rb')
    metadata = datafile.readline().split()
    
    self.setname = metadata[0]
    self.nsys = int(metadata[1])
    self.ndata = int(metadata[2])
        
    # Nuclear
    self.nsig = int(metadata[3])
    
    self.data = []
    self.kin1 = []
    self.kin2 = []
    self.kin3 = []
    self.stat = []
    self.sys  = []
    
    self.error = []
    
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
      for isys in range(7,self.nsys*2):
        if (isys % 2) != 0:
          self.sys[len(self.sys)-1].append(float(linesplit[isys]))

    # Total Error
    for idat in range(0,self.ndata):
      error = self.stat[idat]*self.stat[idat]
      for syserror in self.sys[idat]:
        error += syserror*syserror
      self.error.append(sqrt(error))

    self.printMetadata()


  def printMetadata(self):
    print "************COMMONDATA************"
    print "Setname:", self.setname, "PROC:", self.proc
    print "NDAT:", self.ndata,"NSYS:",self.nsys, "NSIG:",self.nsig


  def printData(self):
    print "IDAT","x","Q^2","CV", "VAR"
    for idat in range(0,self.ndata):
      print idat,self.kin1[idat],self.kin2[idat],self.data[idat], self.error[idat]