import pandas as pd
import numpy as np
import argparse
import time
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from numba import jit

@jit(nopython=True)
def run_coalescence_simulations(num_simulations, m, inv_N, inv_d, inv_k, u):
    '''run coalescence simulations with jit compilation'''
    times = np.empty(num_simulations, dtype=np.int64)
    mutations = np.empty(num_simulations, dtype=np.int64)
    
    for i in range(num_simulations):
        same_deme = True
        reproductions = 0
        while True:
            if same_deme:
                reproductions += 1
                if np.random.random() < inv_k:
                    break
                else:
                    reproductions += 1
                    same_deme = np.random.random() < m**2 + (1-m)**2
                    if same_deme and (np.random.random() < inv_N):
                        break
            else:
                reproductions += 2
                same_deme = np.random.random() < inv_d
                if same_deme and (np.random.random() < inv_N):
                    break
        
        times[i] = reproductions
        mutations[i] = np.random.binomial(2*reproductions, u)
    
    return times, mutations

def coalescent_simulation(k, m, d=100, N=1000, u=1e-6, num_simulations=1000):
    '''run backward coalescent simulation with binomial mutations'''
    times, mutations = run_coalescence_simulations(num_simulations, m, 1.0/N, 1.0/d, 1.0/k, u)
    
    return {'mean_pi': np.mean(mutations), 'se_pi': np.std(mutations) / np.sqrt(num_simulations), 'mean_T': np.mean(times)}

def calculate_exact_coalescence_time(k, m, d=100, N=1000):
    '''exact coalescence time via first-step analysis'''
    alpha = 2*m*(1-m)
    return (1/k + (1-1/k)*(2 + 2*d*alpha)) / (1 - (1-1/k)*(1-1/N))
    # return (1 + (1-1/k)*(1-1/N)*(2*d*alpha)) / (1 - (1-1/k)*(1-1/N))

def process_single_m(args):
    '''process single migration rate for parallel execution'''
    m, k, d, N, u, num_samples = args
    sim_result = coalescent_simulation(k, m, d, N, u, num_samples)
    exact_T = calculate_exact_coalescence_time(k, m, d, N)
    
    return {
        'k': k, 'm': m, 'd': d, 'N': N, 'u': u,
        'empirical_pi': sim_result['mean_pi'], 'empirical_pi_se': sim_result['se_pi'],
        'empirical_T': sim_result['mean_T'], 'exact_T': exact_T, 'exact_pi': 2*u*exact_T
    }

def run_analysis(k=10, d=100, N=1000, u=1e-6, num_samples=1000,
                output_file='results.csv', plot_file='plot.png', create_plot=True, n_processes=None, num_m_points=19):
    '''run coalescent analysis over migration rate grid'''
    m_values = np.logspace(-6, -0.01, num_m_points)
    n_processes = min(cpu_count(), len(m_values)) if n_processes is None else n_processes
    
    print(f'analysis: k={k}, d={d}, N={N}, u={u}, {len(m_values)} rates, {n_processes} processes')
    
    start_time = time.time()
    with Pool(n_processes) as pool:
        results = pool.map(process_single_m, [(m, k, d, N, u, num_samples) for m in m_values])
    
    print(f'completed in {time.time() - start_time:.1f}s')
    
    # save the stuff
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    print(f'results saved to {output_file}')
    
    if create_plot:
        m_fine = np.logspace(-6, 0, 500)
        theory_pi = [2*u*calculate_exact_coalescence_time(k, m, d, N) for m in m_fine]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(m_fine, theory_pi, 'g-', linewidth=3, label='theory')
        ax.errorbar(df['m'], df['empirical_pi'], yerr=df['empirical_pi_se'], fmt='o', color='blue', capsize=5, label='simulation')
        ax.set_xlabel('$m$', fontsize=16)
        ax.set_ylabel('$\\pi_S$', fontsize=16, rotation=0)
        ax.set_title(f'$k$={k}, $d$={d}, $N$={N}, $u$={u}')
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)
        ax.legend()
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f'plot saved to {plot_file}')
    
    return df

def main():
    parser = argparse.ArgumentParser(description='backward coalescent simulation')
    parser.add_argument('--k', type=int, default=10)
    parser.add_argument('--d', type=int, default=100)
    parser.add_argument('--N', type=int, default=1000)
    parser.add_argument('--u', type=float, default=1e-6)
    parser.add_argument('--num-samples', type=int, default=10000000)
    parser.add_argument('--output-file', default='coalescent_results.csv')
    parser.add_argument('--plot-file', default='coalescent_plot.png')
    parser.add_argument('--no-plot', action='store_true')
    parser.add_argument('--processes', type=int, default=None)
    parser.add_argument('--num-m-points', type=int, default=19)
    
    args = parser.parse_args()
    run_analysis(args.k, args.d, args.N, args.u, args.num_samples,
                args.output_file, args.plot_file, not args.no_plot, args.processes, args.num_m_points)

if __name__ == '__main__':
    main()