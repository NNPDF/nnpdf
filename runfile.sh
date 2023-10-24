#!/bin/bash

#GLOBAL PATHS
N3FIT='/data/theorie/thasenac/miniconda3/envs/nnpdf-dev/bin/n3fit'

# GLOBAL VARIABLES
NCORES='4'
WALLTIME='4:00:00'
MEMORY='16000mb'
RUNCARD_NAME='NNPDF40_nnlo_as_01200'
REP_MIN=1
REP_MAX=110


function submit_job() {

    # run setup   
    COMMAND=$PWD/'launch_'$RUNCARD_NAME'.sh'
    OUTPUT_PATH=$PWD/'results'/$RUNCARD_NAME
    NREP=$1

    LAUNCH=$N3FIT' '$PWD'/n3fit/runcards/reproduce_nnpdf40/'$RUNCARD_NAME'.yml '$NREP' -o '$OUTPUT_PATH

    [ -e $COMMAND ] && rm $COMMAND
    echo $LAUNCH >> $COMMAND
    chmod +x $COMMAND

    # submission
    qsub -q multicore -W group_list=theorie -l nodes=1:ppn=$NCORES -l pmem=$MEMORY -l walltime=$WALLTIME $COMMAND
    # cleaning
    rm $COMMAND
}

REP_LIST=($(seq $REP_MIN 1 $REP_MAX))
for REP in ${REP_LIST[@]}
    do
    submit_job $REP
    done