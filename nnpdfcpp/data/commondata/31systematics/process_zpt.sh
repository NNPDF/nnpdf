#!/bin/bash
# This script appends an additional 1 percent systematic error to LHC Z pT datasets.
SETS=(ATLASZPT7TEV ATLASZPT8TEVMDIST ATLASZPT8TEVYDIST CMSZDIFF12)

for s in "${SETS[@]}"; do
	echo "Processing "$s

	# Generate new datafile
	cat ../DATA_$s.dat | awk 'NR == 1 {print $1, $2+1, $3}; NR > 1 {print $0, $6*0.01, 1}' > DATA_$s.dat

	# Generate new systype
	cat ../systypes/SYSTYPE_${s}_DEFAULT.dat | awk 'NR == 1 {print $1+1}; NR > 1 {print}' > SYSTYPE_${s}_DEFAULT.dat
	NSYS=$(wc -l SYSTYPE_${s}_DEFAULT.dat)
	NSYS=$( echo $NSYS | awk '{print $1}')
	echo ${NSYS}"    MULT    SKIP" >> SYSTYPE_${s}_DEFAULT.dat
done

echo "Processing complete"
