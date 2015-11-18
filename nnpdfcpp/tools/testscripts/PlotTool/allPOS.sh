#!/bin/bash
# Usage : ./appPOS.sh <theoryID>

THID=6
PDFSet=NNPDF30_nlo_as_0118

SETS=(POSDYC POSDYCBD POSDYCBDB POSDYCBS POSDYCgBSB POSDYCD POSDYCDB POSDYCS POSDYCSB POSDYD POSDYS POSDYU POSDYUBD POSDYUBDB POSDYUBS POSDYUBSB POSDYUD POSDYUDB POSDYUS POSDYUSB )
mkdir thpred

for f in "${SETS[@]}"; do
	FKconvolute $PDFSet ../../../data/theory_$THID/fastkernel/FK_$f.dat > ./thpred/pr_$f.dat
	./thcompare.py ../../../data/commondata/DATA_$f.dat ./thpred/pr_$f.dat
done
