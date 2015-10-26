#!/usr/bin/python
# nnpydf - nh 08/14

import sys, os, gc, shutil

#nnpydf imports
import thpredictions, commondata, datplot, dat_md
from thpredictions import ThPredictions
from commondata import CommonData
from datplot import *

def genPlots(prefix, dataset, cData, tPred):

  print "Generating plots for: " + dataset
  
  dataPlots = []
  
  # No kinematics, and DIS plots
  dataPlots.append(dataplot_nokin(cData))
  
  # This stuff is a bit unpleasant, but our current kinematics definitions make it neccesary
  # Data plots with all kinematics bins
  if cData.proc.startswith("JET"):
    dataPlots.append(dataplot_dis(cData,2,1,3, logx = True, logy= True, sqrtx = True, xlab="$p_T$", ylab="rapidity", zlab="$E_{CM}$"))
    dataPlots.append(dataplot_rat(cData,2,1,3, logx = True, logy= False, sqrtx = True, xlab="$p_T$", ylab="rapidity", zlab="$E_{CM}$"))
  elif cData.proc.startswith("DIS") and ("H1HERA2HGHY" in cData.setname or "H1HERA2LOWQ2" in cData.setname):
    dataPlots.append(dataplot_dis(cData,1,2,3, logx = True, logy =True, sqrtx = False, xlab="x", ylab="$Q^2$", zlab="$E_{CM}$"))
    dataPlots.append(dataplot_rat(cData,1,2,3, logx = True, logy =False, sqrtx = False, xlab="x", ylab="$Q^2$", zlab="$E_{CM}$"))
  elif cData.proc.startswith("DIS") and "NMC" is not cData.setname:
    dataPlots.append(dataplot_dis(cData,2,1,3, logx = True, logy= True, sqrtx = False, xlab="$Q^2$", ylab="x", zlab="$E_{CM}$"))
    dataPlots.append(dataplot_rat(cData,2,1,3, logx = True, logy= False, sqrtx = False, xlab="$Q^2$", ylab="x", zlab="$E_{CM}$"))
  elif cData.proc.startswith("DYP") and cData.setname.startswith("DYE886R") == False and cData.setname.startswith("ATLASZHIGHMASS49FB") == False:
    dataPlots.append(dataplot_dis(cData,1,2,3, logx = False, logy= True, sqrtx = False, xlab="rapidity", ylab="M^2", zlab="$E_{CM}$"))
    dataPlots.append(dataplot_rat(cData,1,2,3, logx = False, logy= False, sqrtx = False, xlab="rapidity", ylab="M^2", zlab="$E_{CM}$"))
  
  
  # Individual kinematics bins
  if cData.proc.startswith("JET"):
    dataPlots.append(dataplot_kin(cData,2,1, sqrtx = True, logy=True, xlab="$p_T$", ylab="rapidity"))
  elif cData.proc.startswith("DIS") and ("H1HERA2HGHY" in cData.setname or "H1HERA2LOWQ2" in cData.setname):
    dataPlots.append(dataplot_kin(cData,1,2,3, xlab="x", ylab="$Q^2$", zlab="$E_{CM}$"))
  elif cData.proc.startswith("DIS") and "CHARM" not in cData.setname and "NMC" not in cData.setname:
    dataPlots.append(dataplot_kin(cData,2,1,3, xlab="$Q^2$", ylab="x", zlab="$E_{CM}$"))
  elif cData.proc.startswith("EWK") and "WZ" not in cData.setname:
    dataPlots.append(dataplot_kin(cData,1, xlab="(pseudo)rapidity"))
  elif cData.proc.startswith("DYP") and cData.setname.startswith("DYE886R") == False and cData.setname.startswith("ATLASZHIGHMASS49FB") == False:
    dataPlots.append(dataplot_kin(cData,1,2,xlab="rapidity",ylab="$M^2$"))
  
  # Add theory predictions and export
  for plot in dataPlots:
    plot.addTheory(tPred)
    plot.export(prefix)

if (len(sys.argv) != 3):
  print " usage: ", sys.argv[0], " [CommonData file] [THPredictions file]"
  exit(-1)

commonfile = sys.argv[1]
theoryfile = sys.argv[2]

cData = CommonData(commonfile, debug=True)
tPred = ThPredictions(theoryfile)

genPlots("./", cData.setname, cData, tPred)

