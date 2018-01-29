#!/usr/bin/python
# nnpydf - nh 08/14

import sys, os, gc, shutil
import multiprocessing

#nnpydf imports
import thpredictions, commondata, datplot, dat_md
from thpredictions import ThPredictions
from commondata import CommonData
from datplot import *

#############################################################################################################################

def genPlotPage(prefix, dataset, cData, theory):

  print "Generating plots for: " + dataset

  # Make folder
  os.mkdir(prefix)
  
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
  
  # Add theory predictions
  if theory != 0:
    for plot in dataPlots:
      if len(plot.data) == len(theory.CV):
        plot.addTheory(theory)
  
  for plot in dataPlots:
    plot.export(prefix)

#################################################################################################################################

# Make dir for results
plotDir = "./plts/"
if os.path.exists(plotDir):
  shutil.rmtree(plotDir)
os.mkdir(plotDir)

# Fetch datasets
root = "../../../data/commondata/"
datasets = os.listdir(root)


# Load datasets
cData = []
for dataset in datasets:
  if dataset[0:4] == "DATA":
    cData.append(CommonData(root+dataset))

# Make thpredictions if not already present
thDir = "./th/"
thRoot = "../../../data/theory_6/fastkernel/"
if os.path.exists(thDir) != True:
  os.mkdir(thDir)
  for cDat in cData:
    FKPath = thRoot + "FK_"+cDat.setname+".dat"
    if os.path.exists(FKPath) == True:
      print "Making predictions for: ", FKPath
      os.system("FKconvolute NNPDF30_nlo_as_0118 "+FKPath + " > " + thDir + "th_"+cDat.setname+".dat")

    

# Setup thread pool
pool = multiprocessing.Pool()

for cDat in cData:
  prefix = plotDir + cDat.setname + "/"

  # Load ThPredictions (if present)
  thpath = thDir + "th_"+cDat.setname+".dat"
  theory = 0
  if os.path.exists(thpath):
    theory = ThPredictions(thpath)

  pool.apply_async(genPlotPage, [prefix, cDat.setname, cDat, theory])

pool.close()
pool.join()

# Reopen pool
pool = multiprocessing.Pool()

#Export markdown
for cDat in cData:
  prefix = plotDir + cDat.setname + "/"
  dat_md.markdown_data(cDat,prefix)



