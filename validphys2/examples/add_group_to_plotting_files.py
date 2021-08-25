#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 08:41:37 2021

@author: Rosalyn Pearson

Example script for adding a new group to the PLOTTING files.
In this instance we are adding spiderplot_group which is a 
group for making chi2 spider plots. It is the same as 
nnpdf31_process, but jets, dijets and photon are grouped together.
"""
import glob
from reportengine.compat import yaml
import NNPDF as nnpath

datapath = nnpath.get_data_path()
      
# Selecting the PLOTTING files 
dirstring = f"{datapath}/commondata/PLOTTING_"

for file in glob.glob(f"{dirstring}*"):  
    name = file.replace(dirstring, "")
    name = name.replace(".yaml", "")  # NB some are yaml and some are yml
    name = name.replace(".yml", "")
    with open(file) as openfile:
       fileinfo = yaml.round_trip_load(openfile)
       # The following syntax is specific to spiderplot_group
       if "nnpdf31_process" in fileinfo:
           process = fileinfo["nnpdf31_process"]
       else:
           process = "OTHER"
           
       if process in ["JETS", "DIJET", "PHOTON"]:
           process = "JETS/DIJETS/PHOTON"
       fileinfo.update({"spiderplot_group" : f"{process}"})
       
       
    with open(file, "w") as openfile:
        dumpfile = yaml.round_trip_dump(fileinfo, openfile)
        

    


