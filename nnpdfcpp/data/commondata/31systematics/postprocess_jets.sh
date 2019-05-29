#!/bin/bash
# This script takes the results of process_jets and adds extra C-factor errors
SETS=(ATLAS1JET11 CMSJETS11 ATLAS1JET11_SF CMSJETS11_SF)
DIRECTORY=.

for s in "${SETS[@]}"; do
	FILENAME=$s"_CFAC_ERR.dat"
	echo $FILENAME
	paste ../../NNLOCFAC/CF_QCD_$s.dat | awk 'NR > 9 {print 0, 100*$2;};' > $FILENAME
	echo '' | cat - $FILENAME > temp && mv temp $FILENAME

	paste ../DATA_$s.dat $FILENAME | awk 'NR == 1 {print $1, $2+1, $3}; NR > 1 {print}' > DATA_$s.dat
	rm $FILENAME

	# Generate new systype
	cat ../systypes/SYSTYPE_${s}_DEFAULT.dat | awk 'NR == 1 {print $1+1}; NR > 1 {print}' > SYSTYPE_${s}_DEFAULT.dat
	NSYS=$(wc -l SYSTYPE_${s}_DEFAULT.dat)
	NSYS=$( echo $NSYS | awk '{print $1}')
	echo ${NSYS}"    MULT    SKIP" >> SYSTYPE_${s}_DEFAULT.dat
done

echo "Processing complete"
