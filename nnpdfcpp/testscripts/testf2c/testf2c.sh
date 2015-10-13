#!/bin/bash

# Script to compare the nnpdfcpp code predictions for F2c
# with the LH benchmark tables
# We will compare results for FONLL-A, FONLL-B and FONLL-C
# As well as results in the FFN scheme
# Using exactly the same settings as in arXiv:1003.1241
# The tables are also taken from this reference directly

# Some cleaning
rm -rf *~ *.res FK_* comparison

# First of all, set the kinematics of the data
# points that should be used for the check
cp DATA_F2CTEST.dat  ../../../data/F2CTEST/

for hqscheme in "fonlla" "fonllb" "fonllc"
do

# Generate the correct FK tables for the F2c computation
cp testf2c-$hqscheme.ini ../../projects/fkgenerator/run/ 
cd ../../projects/fkgenerator/run/
# Clean just in case
rm -rf testf2c-$hqscheme
# Run locally - just a few points
./prepare_run.sh testf2c-$hqscheme<<EOF
0
EOF

# Now run the code to generate the fk tables
cd testf2c-$hqscheme/DIS
./run.sh
# Put the sigma tables in the FK format
# Move the PDF evolution FK table in the relevant data folder as well
./FKtables.sh<<EOF #> F2cpred.dat
n
EOF
#sort -n F2cpred.dat > F2cpred.tmp
#tail -n 20 F2cpred.tmp > F2cpred.dat 
#rm F2cpred.tmp
#mv F2cpred.dat ../../../../../testscripts/testf2c/.

# Move by hand and change name as suitable
cp ../results/FastestKernel/FK_F2CTEST.dat ../../../../../data/F2CTEST/

# Now go to the testcodes, and get the output of this FKtable
cd ../../../../../testcode/
if [ $hqscheme == "fonllc" ] 
then
    cp ../testscripts/testf2c/toyLH_NNLO.ini ../config/
else
    cp ../testscripts/testf2c/toyLH_NLO.ini ../config/
fi

make clean
make

echo "Check that toyLH .LHgrid are in your LHAPDF path!"

# Compyte the theory predictions using this FK table
# Use the toy LH PDFs as in the original benchmarks
# Get the relevant FK table for check
cp ../../data/F2CTEST/FK_F2CTEST.dat .

if [ $hqscheme == "fonllc" ] 
then
    ./fkcheck toyLH_NNLO.ini FK_F2CTEST.dat
else
    ./fkcheck toyLH_NLO.ini FK_F2CTEST.dat
fi
# Save the results for later comparison with the LH benchmark tables
mv fkcheck.res ../testscripts/testf2c/fkcheck-$hqscheme.res

# Compare with the suitable precomputed benchmark tables
cd ../testscripts/testf2c/


make clean
make
./comparison-$hqscheme

done # End loop over HQ schemes
