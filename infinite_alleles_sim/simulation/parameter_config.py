#!/usr/bin/env python3
"""
Parameter configuration for bottleneck simulation array jobs.
Define your parameter space here.
"""

# Parameter ranges to explore
# PARAMETERS = {
#     'k': [5, 10, 20, 50],  # bottleneck sizes
#     'm': [0.001, 0.01, 0.1, 0.5],  # migration probabilities
# }

PARAMETERS = {
    'k': [10],  # bottleneck sizes
    'm': [0.001, 0.1, 0.4, 0.8],  # migration probabilities
}

# Number of replicates per parameter combination
REPLICATES_PER_COMBINATION = 100

# Total number of jobs will be: len(k) * len(m) * REPLICATES_PER_COMBINATION
# In this case: 4 * 4 * 10000 = 160,000 jobs 