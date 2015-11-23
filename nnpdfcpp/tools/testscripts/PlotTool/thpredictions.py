# ThPredictions in Python
# nh 07/2014
import os, sys

class ThPredictions:
  def __init__(self,source):

    # Verify path
    if os.path.exists(source) == False:
      print "Error: source file" + source + " not found!"
      sys.exit()

    # Open commondata
    datafile = open(source, 'rb')
    metadata = datafile.readline().split()

    self.setname = metadata[1]
    self.pdfset  = metadata[2]
    self.CV = []
    self.error = []

    for line in datafile:
      linesplit = line.split()
      self.CV.append(float(linesplit[1]))
      self.error.append(float(linesplit[2]))

    datafile.close()
    self.printMetadata()

  def printMetadata(self):
    print "************THPREDICTIONS************"
    print "Setname:", self.setname, "PDFSet:", self.pdfset
    print "NDAT:", len(self.CV)

  def printPredictions(self):
    print "CV", "Error"
    for i in range(0,len(self.CV)-1):
      print self.CV[i], self.error[i]