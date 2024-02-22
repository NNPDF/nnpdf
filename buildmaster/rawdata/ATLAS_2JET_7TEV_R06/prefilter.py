"""
This script reads the sys breakdown from the csv hepdata files HEPData-ins1268975-v1-Table_i.csv
and it prints it in a friendly format for the c++ implementation.
There are 3 possible scenarios for the treatment of the systematics, 
denoted by "nominal", "stronger" and "weaker".
Each file HEPData-ins1268975-v1-Table_i.csv (https://www.hepdata.net/record/ins1268975) 
is first split by hand into 3 files containing the informations relative to the 3 scenarios, denoted by
HEPData-ins1268975-v1-Table_i.csv (nominal)
HEPData-ins1268975-v1-Table_i_stronger.csv (stronger)
HEPData-ins1268975-v1-Table_i_weaker.csv (weaker)
These 3 files (3 for each bin) are given as input to the script which produces
bin_i_nominal.dat 
bin_i_stronger.dat
bin_i_weaker.dat
where the data are collected keeping the same structure as in the original hepdata files,
but without any string which makes the implementation of the c++ filter problematic. 
Following the structure of the original hepdata files, the output files are structured as follows
pt	ptmin	ptmax	xsec	stat	sys_1+	sys_1-	..	sys_n+	sys_n-
  
where pt, ptmin, ptmax, xsec and stat are the same in all the 3 of them. 
"""

import numpy as np

nbins = 6
ndata = [21, 21, 19, 17, 8, 4]

for k in range(nbins):
    output = open("bin_" + str(k + 1) + "_nominal.dat", "w+")
    n = ndata[k]

    ############## Nominal
    # Read, pt, ptmin, ptmax, xsec
    pt, ptmin, ptmax, data = np.loadtxt(
        "./HEPData-ins1268975-v1-Table_" + str(k + 7) + ".csv",
        delimiter=",",
        skiprows=19,
        usecols=(0, 1, 2, 3),
        unpack=True,
    )
    # Read stat- (stat+ is always equal to (-1)*stat-)
    stat = np.loadtxt(
        "./HEPData-ins1268975-v1-Table_" + str(k + 7) + ".csv",
        delimiter="%,",
        skiprows=19,
        usecols=(1),
        unpack=True,
    )

    # Read all the sys. (skip the last one which is (-1)*quad)
    sys = []
    for i in range(135):
        a = np.loadtxt(
            "./HEPData-ins1268975-v1-Table_" + str(k + 7) + ".csv",
            delimiter="%,",
            skiprows=19,
            usecols=(i + 2),
            unpack=True,
        )
        sys.append(a)

    for i in range(n):
        output.write(
            str(pt[i])
            + "\t"
            + str(ptmin[i])
            + "\t"
            + str(ptmax[i])
            + "\t"
            + str(data[i])
            + "\t"
            + str(-stat[i])
            + "\t"
        )
        for j in range(135):
            output.write(str(sys[j][i]) + "\t")

        # For the last sys (quad) we have sys-=(-1)*sys+
        # Print it to recover the original format of data files
        output.write(str(-sys[134][i]) + "\t")
        output.write("\n")
    output.close

    ############## Stronger
    output = open("bin_" + str(k + 1) + "_stronger.dat", "w+")
    # Read, pt, ptmin, ptmax, xsec
    pt, ptmin, ptmax, data = np.loadtxt(
        "./HEPData-ins1268975-v1-Table_" + str(k + 7) + "_stronger.csv",
        delimiter=",",
        skiprows=19,
        usecols=(0, 1, 2, 3),
        unpack=True,
    )
    # Read stat- (stat+ is always equal to (-1)*stat-)
    stat = np.loadtxt(
        "./HEPData-ins1268975-v1-Table_" + str(k + 7) + "_stronger.csv",
        delimiter="%,",
        skiprows=19,
        usecols=(1),
        unpack=True,
    )

    # Read all the sys. (skip the last one which is (-1)*quad)
    sys.clear()
    sys = []
    for i in range(117):
        a = np.loadtxt(
            "./HEPData-ins1268975-v1-Table_" + str(k + 7) + "_stronger.csv",
            delimiter="%,",
            skiprows=19,
            usecols=(i + 2),
            unpack=True,
        )
        sys.append(a)

    for i in range(n):
        output.write(
            str(pt[i])
            + "\t"
            + str(ptmin[i])
            + "\t"
            + str(ptmax[i])
            + "\t"
            + str(data[i])
            + "\t"
            + str(-stat[i])
            + "\t"
        )
        for j in range(117):
            output.write(str(sys[j][i]) + "\t")

        # For the last sys (quad) we have sys-=(-1)*sys+
        # Print it to recover the original format of data files
        output.write(str(-sys[116][i]) + "\t")
        output.write("\n")
    output.close

    ############## Weaker
    output = open("bin_" + str(k + 1) + "_weaker.dat", "w+")
    # Read, pt, ptmin, ptmax, xsec
    pt, ptmin, ptmax, data = np.loadtxt(
        "./HEPData-ins1268975-v1-Table_" + str(k + 7) + "_weaker.csv",
        delimiter=",",
        skiprows=19,
        usecols=(0, 1, 2, 3),
        unpack=True,
    )
    # Read stat- (stat+ is always equal to (-1)*stat-)
    stat = np.loadtxt(
        "./HEPData-ins1268975-v1-Table_" + str(k + 7) + "_weaker.csv",
        delimiter="%,",
        skiprows=19,
        usecols=(1),
        unpack=True,
    )

    # Read all the sys. (skip the last one which is (-1)*quad)
    sys.clear()
    sys = []
    for i in range(139):
        a = np.loadtxt(
            "./HEPData-ins1268975-v1-Table_" + str(k + 7) + "_weaker.csv",
            delimiter="%,",
            skiprows=19,
            usecols=(i + 2),
            unpack=True,
        )
        sys.append(a)

    for i in range(n):
        output.write(
            str(pt[i])
            + "\t"
            + str(ptmin[i])
            + "\t"
            + str(ptmax[i])
            + "\t"
            + str(data[i])
            + "\t"
            + str(-stat[i])
            + "\t"
        )
        for j in range(139):
            output.write(str(sys[j][i]) + "\t")

        # For the last sys (quad) we have sys-=(-1)*sys+
        # Print it to recover the original format of data files
        output.write(str(-sys[138][i]) + "\t")
        output.write("\n")

    output.close
