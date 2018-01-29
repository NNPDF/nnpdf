#!/bin/bash

export OMP_NUM_THREADS=10
for I in 10
do
   nice -n 19 ./nnfit $I 140113-r1494-001-ag.ini > outfiles/140113-r1494-001-ag_rep_$I.out & 
done
