"""
This script reads the sys breakdown from the csv hepdata files HEPData-ins1604271-v1-Table_i.csv
and it prints it in a friendly format for the c++ implementation.
Following the structure of the original hepdata files, the output files are structured as follows
pt	pt_min	pt_max	data	stat+	stat-	sys1+	sys1-	....	sysn+	sysn-
  
"""

import numpy as np

nbins=6
ndata=[34, 33, 32, 30, 24, 18]

for k in range(nbins):
   output=open("bin_" + str(k+1) + ".dat", "w+")
   n = ndata[k]
    
   #Read, pt, ptmin, ptmax, xsec
   pt, ptmin, ptmax, data, statp, statm = np.loadtxt("./HEPData-ins1604271-v1-Table_" + str(k+7) + ".csv", delimiter=",", skiprows=11,  
                                        usecols=(0,1,2,3,4,5), unpack=True)
  
   #Read all the sys
   sys=[]
   for i in range(640):
     a = np.loadtxt("./HEPData-ins1604271-v1-Table_" + str(k+7) + ".csv", delimiter=",", skiprows=11, 
                     usecols=(i+6), unpack=True)
     sys.append(a)
     

   for i in range(n):
      output.write(str(pt[i]) + "\t" + str(ptmin[i]) + "\t" + str(ptmax[i]) + "\t" + str(data[i]) 
                   + "\t" + str(statp[i]) + "\t" + str(statm[i]) + "\t")
      for j in range(640):
        output.write(str(sys[j][i]) + "\t")

      output.write("\n")
   output.close

