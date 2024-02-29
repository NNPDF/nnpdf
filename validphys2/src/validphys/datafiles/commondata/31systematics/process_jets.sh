#!/bin/bash
# This script takes scale variation errors for jet data and appends them to datafiles as an additional systematic

SETS=(ATLAS1JET11 ATLASR04JETS2P76TEV ATLASR04JETS36PB CMS1JET276TEV CMSJETS11 CDFR2KT ATLAS1JET11_SF ATLASR04JETS2P76TEV_SF ATLASR04JETS36PB_SF CMS1JET276TEV_SF CMSJETS11_SF CDFR2KT_SF)
DIRECTORY=.

for s in "${SETS[@]}"; do
	# Merge bins
	for i in $DIRECTORY/THSYS-$s*.dat; do
	    # Process $i
	    echo "Processing " $i
	    cat $i | awk ' NR > 9 {print ($2 + sqrt($3*$3))/2.0}' >> $DIRECTORY/THSYS-$s-MERGED.dat
	done
	# Extract data c.v and multiply by % error
	cat ../DATA_$s.dat | awk 'NR > 1 {print $6}' > $DIRECTORY/DAT-$s.dat
	paste $DIRECTORY/DAT-$s.dat $DIRECTORY/THSYS-$s-MERGED.dat | awk '{printf "%4.5e\t%4.5e\n", $1*$2/100.0, $2}' > THSYS-$s.fin
	rm $DIRECTORY/DAT-$s.dat  $DIRECTORY/THSYS-$s-MERGED.dat

	# Add extra line to beginning of final file for formatting and paste against data file
	echo '' | cat - THSYS-$s.fin > temp && mv temp THSYS-$s.fin
	paste ../DATA_$s.dat THSYS-$s.fin | awk 'NR == 1 {print $1, $2+1, $3}; NR > 1 {print}' > DATA_$s.dat
	rm THSYS-$s.fin

	# Generate new systype
	cat ../systypes/SYSTYPE_${s}_DEFAULT.dat | awk 'NR == 1 {print $1+1}; NR > 1 {print}' > SYSTYPE_${s}_DEFAULT.dat
	NSYS=$(wc -l SYSTYPE_${s}_DEFAULT.dat)
	NSYS=$( echo $NSYS | awk '{print $1}')
	echo ${NSYS}"    MULT    SKIP" >> SYSTYPE_${s}_DEFAULT.dat
done

echo "Processing complete"
