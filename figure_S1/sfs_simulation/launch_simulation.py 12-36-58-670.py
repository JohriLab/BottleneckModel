import numpy as np
import pandas as pd
import os
from multiprocessing import Pool
import itertools
import json
import subprocess
import sys
import argparse

# Define parameters
# Hard-coded list of (k, m, d) triples
kmd_pairs = [
    (5000, 0, 5000),
    (5000, 0.0002, 5000),
    (5000, 0.002, 5000),
    (10, 0, 9990)
]
N = 10000000
u = 1e-6
num_samples = 50
num_replicates = 10000  # Total number of replicates
replicates_per_batch = 1000  # Number of replicates per job batch

def run_single_simulation(params):
    """Run a single simulation using subprocess (works well on SLURM)"""
    k, m, d, replicate_idx = params
    
    # Create command to run the simulation
    cmd = [
        sys.executable, "fast_tree_simulation.py",
        "--k", str(k),
        "--m", str(m), 
        "--d", str(d),
        "--N", str(N),
        "--u", str(u),
        "--num-samples", str(num_samples),
        "--folded"
    ]
    
    try:
        # Run the simulation and capture output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return f"SUCCESS: k={k}, m={m}, d={d}, replicate={replicate_idx}"
    except subprocess.CalledProcessError as e:
        return f"ERROR: k={k}, m={m}, d={d}, replicate={replicate_idx} - {e.stderr}"
    except Exception as e:
        return f"ERROR: k={k}, m={m}, d={d}, replicate={replicate_idx} - {str(e)}"

def run_parameter_space(batch_id=None, replicates_per_batch=None):
    """Run simulations across parameter space using multiprocessing
    
    Args:
        batch_id: If provided, only run replicates for this batch (0-indexed)
        replicates_per_batch: Number of replicates per batch (used with batch_id)
    """
    # Determine which replicates to run
    if batch_id is not None and replicates_per_batch is not None:
        start_replicate = batch_id * replicates_per_batch
        end_replicate = min((batch_id + 1) * replicates_per_batch, num_replicates)
        replicate_range = range(start_replicate, end_replicate)
        print(f"Running batch {batch_id}: replicates {start_replicate} to {end_replicate-1}")
    else:
        replicate_range = range(num_replicates)
        print(f"Running all {num_replicates} replicates")
    
    # Create all parameter combinations from (k, m, d) triples
    all_params = []
    for k, m, d in kmd_pairs:
        for replicate_idx in replicate_range:
            all_params.append((k, m, d, replicate_idx))
    
    # Get number of CPUs from SLURM if available, otherwise use system count
    num_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))
    
    print(f"Running {len(all_params)} simulations...")
    print(f"Using {num_cores} CPU cores")
    
    # Use multiprocessing with subprocess calls
    completed = 0
    results = []
    
    with Pool(processes=num_cores) as pool:
        # Use imap to get results as they complete
        for result in pool.imap(run_single_simulation, all_params):
            results.append(result)
            completed += 1
            print(f"Progress: {completed}/{len(all_params)} - {result}")
    
    # Count successes and failures
    successes = sum(1 for r in results if r.startswith("SUCCESS"))
    failures = len(results) - successes
    
    print(f"\nCompleted {completed} simulations.")
    print(f"Successes: {successes}")
    print(f"Failures: {failures}")
    
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run SFS simulations')
    parser.add_argument('--batch-id', type=int, default=None,
                        help='Batch ID (0-indexed) to run specific batch of replicates')
    parser.add_argument('--replicates-per-batch', type=int, default=None,
                        help='Number of replicates per batch (used with --batch-id)')
    
    args = parser.parse_args()
    
    # Use command-line args if provided, otherwise use defaults
    batch_id = args.batch_id
    reps_per_batch = args.replicates_per_batch if args.replicates_per_batch is not None else replicates_per_batch
    
    run_parameter_space(batch_id=batch_id, replicates_per_batch=reps_per_batch)