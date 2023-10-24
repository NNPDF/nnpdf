#!/bin/bash

################## procduction ##################

#GLOBAL PATHS
EVOLVE_N3FIT='/data/theorie/thasenac/miniconda3/envs/nnpdf-dev/bin/evolven3fit_new'
POSTFIT='/data/theorie/thasenac/miniconda3/envs/nnpdf-dev/bin/postfit'
# FKTABLES_PATH='/data/theorie/thasenac/N3PDF/pineko/data/fktables'
DATA_PATH='/data/theorie/thasenac/miniconda3/envs/nnpdf-dev/share/NNPDF/data'

# GLOBAL VARIABLES
NCORES='32'
WALLTIME='48:00:00'
MEMORY='384gb'
RUNCARD_NAME='eko_400'
EKO_NAME='400_def.tar'
THEORY_ID=400

function submit_job_produce() {

    # run setup   
    COMMAND=$PWD/'launch_'$RUNCARD_NAME'_evolve.sh'

    LAUNCH=$EVOLVE_N3FIT' -n '$NCORES' produce_eko '$THEORY_ID' '$PWD'/ekos/'$EKO_NAME' > '$PWD'/'$RUNCARD_NAME'.log' 

    [ -e $COMMAND ] && rm $COMMAND
    echo $LAUNCH >> $COMMAND
    chmod +x $COMMAND

    # submission
    qsub -q smefit -W group_list=smefit -l nodes=1:ppn=$NCORES -l vmem=$MEMORY -l walltime=$WALLTIME $COMMAND
    # cleaning
    rm $COMMAND
}

#submit_job_produce

#####################################################


# ################## local execution ##################
EKO_NAME='400_def.tar'
NREP=100
NCORES='1'
WALLTIME='8:00:00'
MEMORY='16gb'

function submit_job_evolve() {

    # run setup  
    FIT_NAME=$1
    COMMAND=$PWD/'evolve_'$FIT_NAME'.sh'
    RESULT_PATH=$PWD'/results/'$FIT_NAME
    EKO_PATH=$PWD'/ekos/'$EKO_NAME

    LAUNCH=$EVOLVE_N3FIT' evolve '$RESULT_PATH' -l'$EKO_PATH' -f;'$POSTFIT' '$NREP' '$RESULT_PATH'; cp -r '$RESULT_PATH' /data/theorie/thasenac/miniconda3/envs/nnpdf-dev/share/NNPDF/results/'
    
    [ -e $COMMAND ] && rm $COMMAND
    echo $LAUNCH >> $COMMAND
    chmod +x $COMMAND

    # submission
    qsub -q smefit -W group_list=smefit -l nodes=1:ppn=$NCORES -l vmem=$MEMORY -l walltime=$WALLTIME $COMMAND
    # cleaning
    rm $COMMAND
}

submit_job_evolve NNPDF40_nnlo_as_01200