#!/bin/bash

IDX=0
NDAT=0
for var in "$@"
do
    NDAT=$(( $NDAT +$( tail -n 1 $var | awk '{print $1}' ) + 1))
    cDat[$IDX]=$NDAT
    IDX=$(( $IDX+1 ))
done

IDX=-1
for var in "$@"
do
    if [ $IDX -eq -1 ]
    then
	cat $var
    else
	cat $var | awk -v dat="${cDat[$IDX]}" 'NR > 83 {$1=$1+dat;print;}'
    fi
    IDX=$(( $IDX+1 ))
done
