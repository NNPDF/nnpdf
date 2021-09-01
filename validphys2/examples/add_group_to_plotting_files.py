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
       # The following syntax is specific to groups for spider plots
       if name in ["NMCPD_dw_ite", "NMC",
                   "SLACP_dwsh", "SLACD_dw_ite"
                   "BCDMSP_dwdh", "BCDMSD_dw_ite",]:
           process = "DIS NC (fixed target)"
       elif name in ["CHORUSNUPb_dw_ite", "CHORUSNBPb_dw_ite",
                     "NTVNUDMNFe_dw_ite", "NTVNBDMNFe_dw_ite"]:
           process = "DIS CC (fixed target)"
       elif name in ["HERACOMBNCEP460", "HERACOMBNCEP575",
                     "HERACOMBNCEP820", "HERACOMBNCEP920",
                     "ERACOMB_SIGMARED_C", "ERACOMB_SIGMARED_B"]:
           process = "DIS NC (colider)"
       elif name in ["HERACOMBCCEM", "HERACOMBCCEP"]:
           process = "DIS CC (collider)"
       elif name in ["DYE886R_dw_ite", "DYE886P",
                     "DYE605_dw_ite", "DYE906R_dw_ite"]:
           process = "DY (fixed target)"
       elif name in ["CDFZRAP_NEW", "D0ZRAP_40", "D0WMASY"]:
           process = "Tevatron W,Z prod. (incl.)"
       elif name in ["ATLASWZRAP36PB", "ATLASZHIGHMASS49FB",
                     "ATLASLOMASSDY11EXT", "ATLASWZRAP11CC",
                     "ATLASWZRAP11CF", "ATLASDY2D8TEV",
                     "ATLAS_DY_2D_8TEV_LOWMASS", "ATLAS_WZ_TOT_13TEV",
                     "CMSWEASY840PB", "CMSWMASY47FB", "CMSDY2D11",
                     "CMSWMU8TEV", "LHCBZ940PB", "LHCBZEE2FB_40",
                     "LHCBWZMU7TEV", "LHCBWZMU8TEV",
                     "LHCB_Z_13TEV_DIMUON", "LHCB_Z_13TEV_DIELECTRON"]:
           process = "LHC W,Z prod. (incl.)"
       elif name in ["ATLAS_WP_JET_8TEV_PT", "ATLAS_WM_JET_8TEV_PT",
                     "ATLASZPT8TEVMDIST", "ATLASZPT8TEVYDIST",
                     "CMSZDIFF12", "ATLAS_WCHARM_WP_DIFF_7TEV",
                     "ATLAS_WCHARM_WM_DIFF_7TEV", "CMSWCHARMTOT",
                     "CMSWCHARMRAT", "CMS_WCHARM_DIFF_UNNORM_13TEV"]:
           process = r"LHC W,Z prod. ($p_T$ and jets)"
       elif name in ["ATLASTTBARTOT7TEV", "ATLASTTBARTOT8TEV",
                     "ATLAS_TTBARTOT_13TEV_FULLLUMI",
                     "ATLAS_TTB_DIFF_8TEV_LJ_TRAPNORM",
                     "ATLAS_TTB_DIFF_8TEV_LJ_TTRAPNORM",
                     "ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORM",
                     "CMSTTBARTOT5TEV", "CMSTTBARTOT7TEV",
                     "CMSTTBARTOT8TEV", "CMSTTBARTOT13TEV",
                     "CMSTOPDIFF8TEVTTRAPNORM",
                     "CMS_TTBAR_2D_DIFF_MTT_TRAP_NORM",
                     "CMS_TTB_DIFF_13TEV_2016_2L_TRAP",
                     "CMS_TTB_DIFF_13TEV_2016_LJ_TRAP"]:
           process = "LHC top-quark pair prod."
       elif name in ["ATLAS_1JET_8TEV_R06_DEC", "ATLAS_2JET_7TEV_R06",
                     "CMS_2JET_7TEV", "CMS_1JET_8TEV"]:
           process = "LHC jet prod."
       elif name in ["ATLASPHT15_SF"]:
           process = r"LHC direct $\gamma$ prod."
       elif name in ["ATLAS_SINGLETOP_TCH_R_7TEV",
                     "ATLAS_SINGLETOP_TCH_R_13TEV",
                     "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORM",
                     "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORM",
                     "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP_NORM",
                     "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP_NORM",
                     "CMS_SINGLETOP_TCH_TOT_7TEV",
                     "CMS_SINGLETOP_TCH_R_8TEV",
                     "CMS_SINGLETOP_TCH_R_13TEV"]:
           process = "LHC single top prod."

       fileinfo.update({"spiderplot_group" : f"{process}"})

       if name in ["NMCPD_dw_ite", "NMCPD_sh_ite", "NMCPD"]:
           process = r"NMC d/p"
       elif name in ["SLACD_dw_ite", "SLACD_sh_ite", "SLACD"]:
           process = r"SLAC d"
       elif name in ["BCDMSD_dw_ite", "BCDMSD_sh_ite", "BCDMSD"]:
           process = r"BCDMS d"
       elif name in ["CHORUSNUPb_dw_ite", "CHORUSNUPb_sh_ite", "CHORUSNUPb"]:
           process = r"CHORUS $\nu$"
       elif name in ["CHORUSNBPb_dw_ite", "CHORUSNBPb_sh_ite", "CHORUSNBPb"]:
           process = r"CHORUS $\bar\nu$"
       elif name in ["NTVNUDNMFe_dw_ite", "NTVNUDNMFe_sh_ite", "NTVNUDNMFe"]:
           process = r"NuTeV $\nu$"
       elif name in ["NTVNBDNMFe_dw_ite", "NTVNBDNMFe_sh_ite", "NTVNBDNMFe"]:
           process = r"NuTeV $\bar\nu$",    
       elif name in ["DYE886R_dw_ite", "DYE886R_sh_ite", "DYE886R"]:
           process = r"E866/NuSea",
       elif name in ["DYE605_dw_ite", "DYE605_sh_ite" , "DYE605"]:
           process = r"E605",
       elif name in ["DYE906R_dw_ite", "DYE906R_sh_ite", "DYE906R"]:
           process = r"E906/SeaQuest"
       else:
           process = "OTHER" 
       
       fileinfo.update({"nuclear_group" : f"{process}"})
       
    with open(file, "w") as openfile:
        dumpfile = yaml.round_trip_dump(fileinfo, openfile)
        

    


