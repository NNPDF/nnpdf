#!/bin/bash

# Script to compare the nnpdfcpp code predictions for various 
# DIS observebles with the predictions provided by the x-space
# benchmark code.

# Some cleaning
rm -rf *~ *.res FK_* comparison

# First of all, set the kinematics of the data
# points that should be used for the check
cp DATA_DISOBSTEST.dat  ../../../data/DISOBSTEST/.

# Generate the correct FK tables for the F2c computation
cp testDISobs.ini ../../projects/fkgenerator/run/ 
cd ../../projects/fkgenerator/run/
# Clean just in case
rm -rf testDISobs
# Run locally - just a few points
./prepare_run.sh testDISobs.ini<<EOF
0
EOF

# Now run the code to generate the fk tables
cd testDISobs/DIS
./run.sh
# Put the sigma tables in the FK format
# Move the PDF evolution FK table in the relevant data folder as well
./FKtables.sh<<EOF 
n
EOF

# Move by hand and change name as suitable
cp ../results/FastestKernel/FK_DISOBSTEST.dat ../../../../../testcode/.

# Now go to the testcodes, and get the output of this FKtable
cd ../../../../../testcode/
cp ../testscripts/testDISobs/toyLH_NLO.ini ../config/
#cp ../testscripts/testDISobs/toyLH_LO.ini ../config/

make clean
make

echo "Check that toyLH_NLO.LHgrid is in your LHAPDF path!"

# Compyte the theory predictions using this FK table
# Use the toy LH PDFs as in the original benchmarks
# Get the relevant FK table for check
./fkcheck toyLH_NLO.ini FK_DISOBSTEST.dat
#./fkcheck toyLH_LO.ini FK_DISOBSTEST.dat

# Save the results for later comparison with the LH benchmark tables
mv fkcheck.res ../testscripts/testDISobs/.

# Go to the folder of the x-space benchmark code
cd ../../../external/xspace-bench/
cp ../../trunk/nnpdfcpp/projects/fkgenerator/run/testDISobs/results/kinematics/DISkinematics.in .
make clean
make
./fkcheck
mv bench.dat ../../trunk/nnpdfcpp/testscripts/testDISobs/.

# Compare with the suitable precomputed benchmark tables
cd ../../trunk/nnpdfcpp/testscripts/testDISobs

make clean
make
./comparison
