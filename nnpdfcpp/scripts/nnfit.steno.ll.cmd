#!/bin/bash
#
# @ job_name       = nnfit
# @ job_type       = serial
# @ initialdir     = .
# @ error          = outfiles/nnfit-$(jobid)-$(stepid).err
# @ output         = outfiles/nnfit-$(jobid)-$(stepid).out
# @ environment    = COPY_ALL
# @ class          = tier3
# @ resources = ConsumableCpus(1) ConsumableMemory(2gb)
# @ wall_clock_limit = 47:59:59
# @ queue

echo "Submitting replica" $REP
./nnfit $REP <SET INITIALIZATION FILE>
