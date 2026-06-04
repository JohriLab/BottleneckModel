import argparse
from fast_utils_tree_simulation import *

def main():
    parser = argparse.ArgumentParser(description="Tree simulation: take in k, m, d, N as arguments")
    parser.add_argument('--k', type=int, default=10)
    parser.add_argument('--m', type=float, default=0)
    parser.add_argument('--d', type=int, default=100)
    parser.add_argument('--N', type=int, default=1000)
    parser.add_argument('--u', type=float, default=1e-6)
    parser.add_argument('--num-samples', type=int, default=50)
    parser.add_argument('--plot', action='store_true', help='Generate tree plot')
    parser.add_argument('--plot-file', default='tree_simulation', help='Output plot filename')
    parser.add_argument('--folded', action='store_true', help='Fold the SFS')
    # parser.add_argument('--output-file', default='fast_tree_simulation/simulation_results.h5', help='Output file')
    args = parser.parse_args()

    # Run the simulation
    edges_df = simulate_tree(args.k, args.m, args.N, args.d, args.num_samples)

    # Calculate the SFS
    sfs = calc_SFS_from_edges(edges_df, args.num_samples, args.folded)

    # Calculate π
    pi = calc_pi_from_SFS(sfs, args.num_samples, args.u)

    # Append results to the output file
    save_data_jsonl(args.k, args.m, args.d, args.N, args.u, args.num_samples, pi, sfs)


if __name__ == "__main__":
    main()


