#!/bin/bash

echo "=== Simple Bottleneck Simulation Array Job Workflow ==="
echo "This script will:"
echo "1. Generate parameter combinations"
echo "2. Create job commands file"
echo "3. Submit array job"
echo ""

# Step 1: Generate job commands
echo "Step 1: Generating job commands..."
python3 generate_job_commands.py

if [ $? -ne 0 ]; then
    echo "Error: Failed to generate job commands"
    exit 1
fi

echo ""

# Step 2: Submit array job
echo "Step 2: Submitting array job..."
bash simple_submit.sh

if [ $? -ne 0 ]; then
    echo "Error: Failed to submit array job"
    exit 1
fi

echo ""
echo "=== Workflow completed! ==="
echo "Monitor your jobs with: squeue -u $USER"
echo "Check individual job logs in: logs/"
echo "Results will be saved to: bottleneck_simulation_results.csv" 