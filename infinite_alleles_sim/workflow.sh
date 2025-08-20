#! /bin/bash

rm -r simulation/data

python simulation/generate_job_commands.py

# Run the simulation
./simulation/run_jobs_parallel.sh

analysis_workflow.sh