#!/bin/bash
#SBATCH --job-name=bottleneck_array
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --error=logs/array_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jamescre@live.unc.edu
#SBATCH --time=0:11:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
# Array size will be set dynamically by simple_submit.sh

# Load Python module
module load python/3.9.6

# Change to submission directory
cd $SLURM_SUBMIT_DIR

# Get the command for this array task
# SLURM_ARRAY_TASK_ID is automatically set by SLURM
command=$(sed -n "${SLURM_ARRAY_TASK_ID}p" job_commands.txt)

echo "Array task ${SLURM_ARRAY_TASK_ID}: Executing command: $command"
echo "Start time: $(date)"

# Execute the command
# Execute the command 10 times in a loop
for i in {1..10}; do
    echo "Array task ${SLURM_ARRAY_TASK_ID}: Iteration $i/10"
    eval $command
done

echo "Array task ${SLURM_ARRAY_TASK_ID}: Completed"
echo "End time: $(date)" 