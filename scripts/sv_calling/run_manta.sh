#!/bin/bash
#SBATCH --job-name=manta_sv
#SBATCH --error=/path/to/project/logs/manta.%J.err
#SBATCH --output=/path/to/project/logs/manta.%J.out
#SBATCH --partition=newcompute
#SBATCH --qos=common
#SBATCH --account=common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# Initialize Conda from the installation available on the HPC.
CONDA_BASE="$(conda info --base)" || exit 1
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate manta

cd /path/to/project/analysis/sv_calling
./manta_results/runWorkflow.py -j 8

echo "Manta completed: $(date)"
