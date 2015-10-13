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

def genPlotPage(prefix, dataset, cData):

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
  
  # Find and add theory predictions
  i=0
  for arg in sys.argv:
    if i > 0:
      root = arg + "thpredictions/"
      theory = ThPredictions(root+dataset+".dat")
      
      # Add theory predictions
      for plot in dataPlots:
        if len(plot.data) == len(theory.CV):
          plot.addTheory(theory)
    
    i=i+1
  
  for plot in dataPlots:
    plot.export(prefix)

  # Export markdown


#################################################################################################################################

# Make dir for results
plotDir = "./plts/"
if os.path.exists(plotDir):
  shutil.rmtree(plotDir)
os.mkdir(plotDir)

# Fetch datasets
root = sys.argv[1] + "filter/"
datasets = os.listdir(root)

# Load datasets
cData = []
for dataset in datasets:
  cData.append(CommonData(root+dataset+"/DATA_"+dataset+".dat"))

# Setup thread pool
pool = multiprocessing.Pool()

for cDat in cData:
  prefix = plotDir + cDat.setname + "/"
  pool.apply_async(genPlotPage, [prefix, cDat.setname, cDat])

pool.close()
pool.join()

# Reopen pool
pool = multiprocessing.Pool()

#Export markdown
for cDat in cData:
  print "Exporting Markdown for: " + cDat.setname
  prefix = plotDir + cDat.setname + "/"
  dat_md.markdown_data(cDat,prefix)



