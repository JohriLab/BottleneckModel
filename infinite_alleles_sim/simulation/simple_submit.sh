#!/bin/bash

# Check if job_commands.txt exists
if [ ! -f "job_commands.txt" ]; then
    echo "Error: job_commands.txt not found. Please run generate_job_commands.py first."
    exit 1
fi

# Get total number of commands
total_commands=$(wc -l < job_commands.txt)

echo "Total commands: $total_commands"

# Create logs directory
mkdir -p logs

# Submit array job with the correct size
echo "Submitting array job with ${total_commands} tasks..."
sbatch --array=1-${total_commands}%1000 simple_array_job.sh

echo "Array job submitted!"
echo "Monitor with: squeue -u $USER"
echo "Check logs in: logs/" 