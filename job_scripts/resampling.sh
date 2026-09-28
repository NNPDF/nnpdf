#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=rs1x1000
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=resample_%j.out
#SBATCH --error=resample_%j.err
#SBATCH --export=NONE

set -e

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

# The original fit dir in runcards/examples was deleted; the surviving copy in
# the results dir is used directly (post.sh no longer needs the cp -r step).
FIT=/home/dakshanshchawda/miniconda3/envs/environment_nnpdf/share/NNPDF/results/bay_1x1000

# Redraw 100 pseudo-replicas from the trained posterior (post reset_random fix).
# --start-replica 1 matches the existing 1-indexed replica dirs, so the new
# samples overwrite the exportgrids of replica_1..replica_100.
python /home/dakshanshchawda/nnpdf/n3fit/src/n3fit/bnn_inference.py sample \
    --runcard /home/dakshanshchawda/nnpdf/n3fit/runcards/examples/bay_1x1000.yml \
    --fit-dir "$FIT" \
    --samples 100 \
    --start-replica 1

# Remove the stale degenerate pseudo-replicas from the pre-fix run; postfit
# would otherwise select from them. They are byte-identical copies of one PDF,
# and replica_1..100 keep a copy of the (shared) weights.weights.h5.
for i in $(seq 101 1200); do
    rm -rf "$FIT/nnfit/replica_$i"
done

# Clear old evolution + postfit output, otherwise evolven3fit skips replicas
# that already have a .dat file and postfit refuses to rerun.
rm -rf "$FIT/postfit" "$FIT/evolven3fit"
rm -f "$FIT/evolven3fit.log"
rm -f "$FIT"/nnfit/replica_*/bay_1x1000.dat

echo "Resampling + cleanup done: $(ls -d "$FIT"/nnfit/replica_* | wc -l) replica dirs remain"
