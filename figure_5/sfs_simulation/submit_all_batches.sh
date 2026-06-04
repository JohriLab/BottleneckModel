#!/bin/bash

# Automated script to submit all simulation batches
# This script calculates the number of batches needed and submits them all

# Configuration - modify these as needed
REPLICATES_PER_BATCH=1000
TOTAL_REPLICATES=10000

# Calculate number of batches needed
NUM_BATCHES=$(( (TOTAL_REPLICATES + REPLICATES_PER_BATCH - 1) / REPLICATES_PER_BATCH ))

echo "Submitting $NUM_BATCHES jobs..."
echo "Each job will run $REPLICATES_PER_BATCH replicates"
echo "Total replicates: $TOTAL_REPLICATES"

# Create logs directory if it doesn't exist
mkdir -p logs

# Submit jobs using SLURM job array (most efficient)
# This submits all batches as a single job array
# Export REPLICATES_PER_BATCH and pass it to all array tasks
export REPLICATES_PER_BATCH
sbatch --export=REPLICATES_PER_BATCH --array=0-$((NUM_BATCHES-1)) submit_sfs.sh

echo "Submitted job array with $NUM_BATCHES tasks"
echo "Each task will process $REPLICATES_PER_BATCH replicates"
echo "Check status with: squeue -u $USER"
