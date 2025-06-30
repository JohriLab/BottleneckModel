import numpy as np
import pandas as pd
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Bottleneck simulation with configurable parameters')
parser.add_argument('--k', type=int, default=5, help='Number of founders (bottleneck size)')
parser.add_argument('--m', type=float, default=0.01, help='Migration probability')
parser.add_argument('--replicate', type=int, default=1, help='Replicate number')
args = parser.parse_args()

# Parameters
d = 100    # number of demes
N = int(1e3)   # carrying capacity per deme
k = args.k    # number of founders (bottleneck size)
m = args.m  # migration probability
u = 1e-7  # mutation rate
num_generations = 20*d*k  # Convert to integer

print(f"Starting replicate {args.replicate} with parameters k = {k}, m = {m}")

Ne = d*N # precalculate Ne for efficiency

# Initialize d x N matrix A
A = np.zeros((d, N), dtype=int)

fixation_tracker = 0

# Pre-allocate arrays to avoid repeated allocation
M = np.zeros((d, k), dtype=int)
D = np.zeros(d, dtype=int)
D_m = np.zeros(d, dtype=int)
C_K = np.zeros((d, k), dtype=int)
K = np.zeros((d, k), dtype=int)
C_A = np.zeros((d, N), dtype=int)
A_new = np.zeros((d, N), dtype=int)
mutation_mask = np.zeros((d, N), dtype=int)

# Pre-allocate index arrays for migration operations
max_migrants = d * k  # Maximum possible number of migrants
i_mig = np.zeros(max_migrants, dtype=int)
j_mig = np.zeros(max_migrants, dtype=int)
i_non = np.zeros(max_migrants, dtype=int)
j_non = np.zeros(max_migrants, dtype=int)

for t in range(num_generations):

    # 1. Generate migration matrix M (d × k)
    # Each element m_ij represents outcome of Bernoulli trial with probability m
    M[:] = np.random.binomial(1, m, size=(d, k))

    # 2. Generate parental deme choice vector D (length d)
    # Each element d_i drawn uniformly from {1, ..., d}
    D[:] = np.random.randint(0, d, size=d)

    # 3. Generate migrant source deme choice vector D_m (length d)
    # For single deme case (d=1), D_m = D since there's only one deme
    if d == 1:
        D_m[:] = D
    else:
        # OPTIMIZED VERSION - the list comprehension is extremely slow
        D_m[:] = np.random.randint(0, d-1, size=d)
        D_m[D_m >= D] += 1  # Shift values to avoid choosing same deme

    # 4. Generate parent choice matrix within demes C^K (d × k)
    # Each element c_ij drawn uniformly from {1, ..., N}
    C_K[:] = np.random.randint(0, N, size=(d, k))

    # 5. Create the founder population matrix K (d × k) - VECTORIZED
    # Use advanced indexing instead of nested loops
    migrant_indices = M == 1
    non_migrant_indices = M == 0
    
    # For non-migrants: use parental deme
    i_non, j_non = np.where(non_migrant_indices)
    K[i_non, j_non] = A[D[i_non], C_K[i_non, j_non]]
    
    # For migrants: use migrant source deme
    i_mig, j_mig = np.where(migrant_indices)
    K[i_mig, j_mig] = A[D_m[i_mig], C_K[i_mig, j_mig]]

    # 6. Generate parent choice matrix C^A for individuals drawn from propagule
    C_A[:] = np.random.randint(0, k, size=(d, N))

    # 7. Create new population matrix A' via Wright-Fisher reproduction - VECTORIZED
    # Use advanced indexing instead of nested loops
    i_indices, j_indices = np.meshgrid(np.arange(d), np.arange(N), indexing='ij')
    A_new[:] = K[i_indices, C_A]

    # 8. Generate mutation mask - each element has probability u of mutating
    mutation_mask[:] = np.random.binomial(1, u, size=(d, N))

    # 9. Apply mutations: where mask is 1, assign new unique alleles     
    num_mutations = np.sum(mutation_mask)
    if num_mutations > 0:
        A_new[mutation_mask == 1] = 1

    # 10. Check if allele is fixed
    if np.sum(A_new) == Ne:
        A_new.fill(0)  # More efficient than creating new array
        print(f"!!!Allele fixed at generation {t+1}!!!")
        fixation_tracker += 1
    
    # 11. Update A for next generation - use copy to avoid creating new array
    A[:] = A_new

# Final calculations
# piT
allele_count = np.sum(A)
p = allele_count / Ne
piT_final = 2 * p * (1 - p)

# piS
row_count = np.sum(A_new, axis=1) / N
piS_deme = 2 * row_count * (1 - row_count)
piS_final = np.mean(piS_deme)

# # FST
# FST_final = 1 - piS_final / piT_final

#Saving results
final_results_df = pd.DataFrame({
    'piS': [piS_final], 
    'piT': [piT_final],
    'allele_count': [allele_count],
    # 'FST': [FST_final],
    'k': [k],
    'm': [m]
})
final_results_df.to_csv(f'bottleneck_simulation_results.csv', mode='a', header=False, index=False)