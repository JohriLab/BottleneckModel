#!/bin/bash

#SBATCH -n 1
#SBATCH --job-name=sfs_sim
#SBATCH --cpus-per-task=40
#SBATCH --output=logs/out_%A_%a.out
#SBATCH --error=logs/out_%A_%a.err
#SBATCH --mem=40g
#SBATCH -t 01:00:00
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jamescre@live.unc.edu

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules
module load python/3.9.6

# Get batch_id from SLURM_ARRAY_TASK_ID if using job arrays, or from command line argument
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    BATCH_ID=$SLURM_ARRAY_TASK_ID
else
    BATCH_ID=${1:-0}
fi

# Get replicates_per_batch from environment variable, command line argument, or use default
REPLICATES_PER_BATCH=${REPLICATES_PER_BATCH:-${2:-10}}

# Run the simulation with batch parameters
python launch_simulation.py --batch-id $BATCH_ID --replicates-per-batch $REPLICATES_PER_BATCH