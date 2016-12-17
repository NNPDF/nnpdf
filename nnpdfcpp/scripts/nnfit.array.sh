#!/bin/sh
#$ -cwd
#$ -t 1-100
#$ -l h_rt=03:00:00

# Initialise the module environment
. /etc/profile.d/modules.sh

export ROOTSYS=/exports/applications/apps/root/5.22.00/
FILENAME=$SGE_TASK_ID".time"
/usr/bin/time -v -o $FILENAME ./nnfit $SGE_TASK_ID 

