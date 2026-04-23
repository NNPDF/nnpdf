#!/bin/bash
#SBATCH --job-name=post
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=post_%j.out
#SBATCH --error=post_%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /path/to/runcard

evolven3fit evolve Basic_runcard_bayesian

postfit 100 Basic_runcard_bayesian

cp -r Basic_runcard_bayesian /path/to/environment_nnpdf/share/NNPDF/results

vp-comparefits "Basic_runcard_bayesian" "Basic_runcard_normal" \
    --title "Comparison Report" \
    --author "<NAME>" \
    --keywords "Bayesian_report_40-like" \


