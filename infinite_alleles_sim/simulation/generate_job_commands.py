#!/usr/bin/env python3
"""
Generate all command combinations for array jobs.
This creates a file with one command per line for SLURM array jobs.
"""

import itertools
from parameter_config import PARAMETERS, REPLICATES_PER_COMBINATION

def generate_commands():
    """Generate all command combinations and write to file."""
    
    # Get parameter names and values
    param_names = list(PARAMETERS.keys())
    param_values = list(PARAMETERS.values())
    
    # Generate all parameter combinations
    param_combinations = list(itertools.product(*param_values))
    
    # Create commands file
    with open('simulation/job_commands.txt', 'w') as f:
        for combo in param_combinations:
            # Create parameter dict for this combination
            params = dict(zip(param_names, combo))
            
            # Generate replicates for this parameter combination
            for replicate in range(1, REPLICATES_PER_COMBINATION + 1):
                # Build command string
                cmd_parts = ['python bottleneck_simulation.py']
                for param_name, param_value in params.items():
                    cmd_parts.append(f'--{param_name} {param_value}')
                cmd_parts.append(f'--replicate {replicate}')
                
                command = ' '.join(cmd_parts)
                f.write(command + '\n')
    
    total_jobs = len(param_combinations) * REPLICATES_PER_COMBINATION
    print(f"Generated {total_jobs} commands in job_commands.txt")
    print(f"Parameter combinations: {len(param_combinations)}")
    print(f"Replicates per combination: {REPLICATES_PER_COMBINATION}")
    
    # Show first few commands as example
    print("\nFirst 5 commands:")
    with open('simulation/job_commands.txt', 'r') as f:
        for i, line in enumerate(f):
            if i < 5:
                print(f"  {i+1}: {line.strip()}")
            else:
                break

if __name__ == "__main__":
    generate_commands() 