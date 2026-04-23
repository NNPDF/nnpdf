#!/bin/bash

SETUP_JOB=$(sbatch setup.sh | awk '{print $4}')
FIT_JOB=$(sbatch --dependency=afterok:$SETUP_JOB fit_array.sh | awk '{print $4}')
sbatch --dependency=afterok:$FIT_JOB post.sh
