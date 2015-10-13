#!/bin/bash

# Script to test the PDF evolution using the fkgenerator project
# by comparing with the hoppet program

# Some cleaning
rm -rf *~ LHApred.dat

# First of all - need to run fkgenerator to generate the FK tables
# for the PDF evolution at the x,Q2 values defined in 
# data/PDFEVOLTEST/DATA_PDFEVOLTEST.dat
# The input file to be used is testevol.ini
cp DATA_PDFEVOLTEST.dat  ../../../data/PDFEVOLTEST/

# Generate the correct FK tables for the PDF evolution test
cp testevol.ini ../../projects/fkgenerator/run/
cd ../../projects/fkgenerator/src/commons
mv LHAgrid.h LHAgrid.h.bkp
cp ../../../../testscripts/testevol/LHAgrid.h .
cd ../../run

# Clean just in case
rm -rf testevol

# Run locally - just a few points
./prepare_run.sh testevol.ini<<EOF
0
EOF
# Use different common for this test
cd ../src/commons
mv LHAgrid.h.bkp LHAgrid.h
cd ../../run

# Copy DIS kinematics into LHA kinematics
cd testevol/results/kinematics
cp DISkinematics.in LHAkinematics.in

# Now run the code to generate the fk tables
cd ../../LHA
sed -ie 's/5000/'5'/' job-1.sh
rm *.she
./run.sh

# Put the sigma tables in the FK format
# Move the PDF evolution FK table in the relevant data folder as well
./FKtables.sh<<EOF
n
EOF
cp ../results/FastestKernel/LHApred.dat ../../../../../testscripts/testevol/.

# Now run hoppet exactly on the same points
# (HOPPET must have been installed as discussed in the README file)
echo "Now running HOPPET"
echo "Remember that hoppet needs to be installed in the same path as the src"
cd ../../../../../external/hoppet/hoppet-1.1.5/example_f77

cp Makefile_test Makefile
make fkgenerator_test

./fkgenerator_test
