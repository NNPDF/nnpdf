#!/bin/bash
# This script takes scale variation errors for jet data and appends them to datafiles as an additional systematic
# Some further actions are requires after running it, specifically the Nsys needs to be incremented by 1 in the commondata
# and the SYSTYPE file must be correspondingly modified

SETS=(ATLAS1JET11 ATLASR04JETS2P76TEV ATLASR04JETS36PB CMS1JET276TEV CMSJETS11 CDFR2KT)
DIRECTORY=.

for s in "${SETS[@]}"; do
	# Merge bins
	rm $DIRECTORY/THSYS-$s-MERGED.dat
	for i in $DIRECTORY/THSYS-$s*.dat; do
	    # Process $i
	    echo "Processing " $i
	    cat $i | awk ' NR > 9 {print $2 + sqrt($3*$3)}' >> $DIRECTORY/THSYS-$s-MERGED.dat
	done
	# Extract data c.v and multiply by % error
	cat ../DATA_$s.dat | awk 'NR > 1 {print $6}' > $DIRECTORY/DAT-$s.dat
	paste $DIRECTORY/DAT-$s.dat $DIRECTORY/THSYS-$s-MERGED.dat | awk '{printf "%4.5e\t%4.5e\n", $1*$2/100.0, $2}' > THSYS-$s.fin

	# Add extra line to beginning of final file for formatting and paste against data file
	echo '+1sys' | cat - THSYS-$s.fin > temp && mv temp THSYS-$s.fin
	paste ../DATA_$s.dat THSYS-$s.fin > DATA_$s.dat
done

echo "Processing complete, please now edit systype and NSys"
