import argparse
import random
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
from utils_tree_simulation import *



def main():
    parser = argparse.ArgumentParser(description="Tree simulation: take in k, m, d, N as arguments")
    parser.add_argument('--k', type=int, default=10)
    parser.add_argument('--m', type=float, default=0.01)
    parser.add_argument('--d', type=int, default=100)
    parser.add_argument('--N', type=int, default=1000)
    parser.add_argument('--u', type=float, default=1e-6)
    parser.add_argument('--num-samples', type=int, default=50)
    parser.add_argument('--plot', action='store_true', help='Generate tree plot')
    parser.add_argument('--plot-file', default='tree_simulation', help='Output plot filename')
    args = parser.parse_args()

    # Run the simulation
    edges_df = simulate_tree(args.k, args.m, args.N, args.d, args.num_samples)
    
    # Save results
    edges_df.to_csv(args.plot_file+'/tree_edges.csv', index=False)
    print(f"Simulation complete! Tree saved to tree_edges.csv")
    print(f"Total edges: {len(edges_df)}")
    
    # Plot the tree if requested
    sfs = calc_SFS_from_edges(edges_df, args.num_samples)  
    plot_tree_and_sfs(edges_df, sfs, save_file=args.plot_file+'/combined_plot.png')

    E_T_vectorized = calc_T_matrix_vectorized(edges_df, args.num_samples)
    print("E[T] =", E_T_vectorized)
    print("E[pi] =", E_T_vectorized*2*args.u)
    print("E[pi] (SFS) =", calc_pi_from_SFS(sfs, args.num_samples, args.u))

if __name__ == "__main__":
    main()


