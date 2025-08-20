#!/bin/bash

# Run commands from job_commands.txt with max 10 concurrent processes

cd simulation

while read line; do
    # Wait if we have 10 or more background jobs
    while [ $(jobs -r | wc -l) -ge 10 ]; do
        sleep 0.1
    done
    
    # Start command in background
    eval "$line" &
done < job_commands.txt

# Wait for all jobs to finish
wait 