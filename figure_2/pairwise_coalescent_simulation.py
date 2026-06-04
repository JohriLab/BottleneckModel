import pandas as pd
import numpy as np
import argparse
import time
import os
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from numba import jit

# Three labeled ticks on log_10(m) (avoids sparse default / culling).
_LOG_M_AXIS_TICKS = (1e-5, 1e-3, 1e-1)


def _apply_log_xaxis_three_ticks(ax, tick_positions=None):
    ax.set_xscale('log')
    ax.set_xticks(tick_positions if tick_positions is not None else _LOG_M_AXIS_TICKS)


@jit(nopython=True)
def run_coalescence_simulations(num_simulations, m, N, d, k, inv_N, inv_d, inv_k, u, same_deme=False):
    '''run coalescence simulations with jit compilation'''
    times = np.empty(num_simulations, dtype=np.int64)
    mutations = np.empty(num_simulations, dtype=np.int64)

    for i in range(num_simulations):
        reproductions = 0
        current_same_deme = same_deme
        if not current_same_deme:
            current_same_deme = np.random.random() < inv_d
            
        while True:
            reproductions += 1
            if current_same_deme:
                if np.random.random() < inv_k:
                    break
                else:
                    current_same_deme = np.random.random() < m**2 + (1-m)**2
                    if current_same_deme and (np.random.random() < inv_N):
                        break
            else:
                current_same_deme = np.random.random() < inv_d
                if current_same_deme and (np.random.random() < inv_N):
                    break
        
        times[i] = reproductions
        mutations[i] = np.random.binomial(2 * reproductions, u)
    
    return times, mutations

def coalescent_simulation(k, m, d=100, N=1000, u=1e-6, num_simulations=1000):
    '''run backward coalescent simulation with binomial mutations'''
    times, mutations = run_coalescence_simulations(num_simulations, m, N, d, k, 1/N, 1/d, 1/k, u, same_deme=False)

    times_s, mutations_s = run_coalescence_simulations(num_simulations, m, N, d, k, 1/N, 1/d, 1/k, u, same_deme=True)
    
    # Calculate FST with safeguards against division by zero
    mean_mutations = np.mean(mutations)
    mean_mutations_s = np.mean(mutations_s)
    fst = 1 - mean_mutations_s / mean_mutations if mean_mutations > 0 else 0
    
    return {'mean_piT': mean_mutations, 'se_piT': np.std(mutations) / np.sqrt(num_simulations), 'mean_T_T': np.mean(times), 'mean_pi_s': mean_mutations_s, 'se_pi_s': np.std(mutations_s) / np.sqrt(num_simulations), 'mean_T_s': np.mean(times_s), 'FST': fst}

def calculate_expected_times(k, m, d, N, u):
    alpha = 2*m*(1-m)
    # beta = (1 - m)**2 + 2*m*(1-m)*(d-2)/(d-1) + m**2*(d-2)/(d-1)
    beta = 1

    Q = np.array([[(1-1/k)*(1-1/N)*(1-alpha), (1-1/k)*(1-1/N) * alpha],
                [(1/d)*(1-1/N), (1-1/d)]])

    I = np.eye(2)
    N_matrix = np.linalg.inv(I - Q)

    # Calculate expected absorption times
    expected_times = N_matrix @ np.ones((2, 1))

    T_S = expected_times[0, 0]
    T_B = expected_times[1, 0]
    T_T = T_S/d + T_B*(1-1/d)

    p = 1/d
    var_T = (2*N_matrix - I) * expected_times - expected_times * expected_times
    var_T_S = var_T[0, 0]
    var_T_B = var_T[1, 0]
    var_T_T = p*var_T_S + (1 - p)*var_T_B + p*(1 - p)*(T_S - T_B)**2

    e_pi_S = 2*u*T_S
    e_pi_T = 2*u*T_T

    var_pi_s = 2*u*(1-u)*T_S + 4*u**2*var_T_S
    var_pi_T = 2*u*(1-u)*T_T + 4*u**2*var_T_T

    FST = 1 - T_S/T_T

    return T_S, T_T, FST, e_pi_S, e_pi_T, var_pi_s, var_pi_T

def process_single_m(args):
    '''process single migration rate for parallel execution'''
    m, k, d, N, u, num_samples = args
    sim_result = coalescent_simulation(k, m, d, N, u, num_samples)
    exact_T_s, exact_T, FST, e_pi_S, e_pi_T, var_pi_s, var_pi_T = calculate_expected_times(k, m, d, N, u)

    return {
        'k': k, 'm': m, 'd': d, 'N': N, 'u': u,
        'empirical_piT': sim_result['mean_piT'], 'empirical_piT_se': sim_result['se_piT'],
        'empirical_pi_s': sim_result['mean_pi_s'], 'empirical_pi_s_se': sim_result['se_pi_s'],
        'empirical_T_T': sim_result['mean_T_T'], 'empirical_T_s': sim_result['mean_T_s'],
        'exact_T': exact_T, 'exact_T_s': exact_T_s, 'FST': FST,
        'exact_pi': e_pi_T, 'exact_pi_s': e_pi_S,
        'exact_pi_s_var': var_pi_s, 'exact_pi_T_var': var_pi_T
    }

def create_plot_with_transformed_x(df, k_values, d, N, u, num_samples, plot_file, 
                                   presentation_mode, log_scale, x_transform_func, 
                                   x_label, x_suffix, plots_dir):
    '''Helper function to create plots with transformed x-axis'''
    m_fine = np.logspace(-6, 0, 500)
    
    # Calculate all theoretical values to determine shared y-axis limits
    all_theory_fst = []
    all_empirical_fst = []
    all_theory_pi_s = []
    all_empirical_pi_s = []
    all_theory_pi_t = []
    all_empirical_pi_t = []
    all_theory_pi_s_var = []
    all_theory_pi_t_var = []
    all_empirical_pi_s_var = []
    all_empirical_pi_t_var = []
    
    for k in k_values:
        df_k = df[df['k'] == k]
        for m in m_fine:
            pi_s_time, pi_t_time, fst, e_pi_S, e_pi_T, var_pi_s, var_pi_T = calculate_expected_times(k, m, d, N, u)
            all_theory_pi_s.append(e_pi_S)
            all_theory_pi_t.append(e_pi_T)
            all_theory_fst.append(fst)
            all_theory_pi_s_var.append(var_pi_s)
            all_theory_pi_t_var.append(var_pi_T)
        all_empirical_pi_s.extend(df_k['empirical_pi_s'].values)
        all_empirical_pi_t.extend(df_k['empirical_piT'].values)
        all_empirical_fst.extend(df_k['FST'].values)
        all_empirical_pi_s_var.extend(df_k['empirical_pi_s_se'].values)
        all_empirical_pi_t_var.extend(df_k['empirical_piT_se'].values)
    
    # Calculate shared y-axis limits
    def calc_ylim(theory_vals, empirical_vals, min_val=0):
        all_vals = theory_vals + empirical_vals
        val_min = min(all_vals)
        val_max = max(all_vals)
        margin = (val_max - val_min) * 0.05
        return (max(min_val, val_min - margin), val_max + margin)
    
    # Convert theoretical variances to standard errors for y-axis calculation
    all_theory_pi_s_se = np.sqrt(np.array(all_theory_pi_s_var)) / np.sqrt(num_samples)
    all_theory_pi_t_se = np.sqrt(np.array(all_theory_pi_t_var)) / np.sqrt(num_samples)
    
    fst_ylim = calc_ylim(all_theory_fst, all_empirical_fst, min_val=0)
    theory_pi_s_upper = (np.array(all_theory_pi_s) + all_theory_pi_s_se).tolist()
    empirical_pi_s_upper = (np.array(all_empirical_pi_s) + np.array(all_empirical_pi_s_var)).tolist()
    pi_s_ylim = (0, max(theory_pi_s_upper + empirical_pi_s_upper) * 1.05)
    theory_pi_t_upper = (np.array(all_theory_pi_t) + all_theory_pi_t_se).tolist()
    empirical_pi_t_upper = (np.array(all_empirical_pi_t) + np.array(all_empirical_pi_t_var)).tolist()
    theory_pi_t_lower = (np.array(all_theory_pi_t) - all_theory_pi_t_se).tolist()
    empirical_pi_t_lower = (np.array(all_empirical_pi_t) - np.array(all_empirical_pi_t_var)).tolist()
    lower_bound_pi_t = min(theory_pi_t_lower + empirical_pi_t_lower)
    upper_bound_pi_t = max(theory_pi_t_upper + empirical_pi_t_upper)
    pi_t_ylim = (lower_bound_pi_t - 0.05 * (upper_bound_pi_t - lower_bound_pi_t),
                 upper_bound_pi_t + 0.05 * (upper_bound_pi_t - lower_bound_pi_t))
    
    # Plot settings
    if presentation_mode:
        font_size = 14
        line_width = 2.5
        markersize = 5
        ribbon_alpha = 0.3
        pi_s_color = "#D55E00"
        pi_t_color = "#009E73"
        fst_color = "#CC79A7"
        emp_color = "#000000"
    else:
        font_size = 16
        line_width = 2.0
        markersize = 4
        ribbon_alpha = 0.2
        pi_s_color = "#D55E00"
        pi_t_color = "#009E73"
        fst_color = "#CC79A7"
        emp_color = "#000000"
    
    plt.rcParams.update({
        'font.size': font_size,
        'axes.titlesize': font_size,
        'axes.labelsize': font_size,
        'xtick.labelsize': font_size,
        'ytick.labelsize': font_size,
        'legend.fontsize': font_size,
        'lines.linewidth': line_width,
        'lines.markersize': markersize,
        'axes.linewidth': 0.8,
        'grid.linewidth': 0.5,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    unit_size = 3.2
    fig, axes = plt.subplots(3, len(k_values), figsize=(unit_size*len(k_values), unit_size*3))
    
    if len(k_values) == 1:
        axes = axes.reshape(-1, 1)

    for r in range(3):
        for c in range(len(k_values)):
            ax = axes[r, c]
            try:
                ax.set_box_aspect(1)
            except Exception:
                ax.set_aspect('equal', adjustable='box')
            ax.tick_params(width=0.8, length=3)
            if r < 2:
                ax.tick_params(labelbottom=False)
            for spine in ax.spines.values():
                spine.set_linewidth(0.8)
    
    for col, k in enumerate(k_values):
        df_k = df[df['k'] == k]
        
        # Calculate theoretical values for current k
        theory_pi_s = []
        theory_pi_t = []
        theory_fst = []
        theory_pi_s_var = []
        theory_pi_t_var = []
        for m in m_fine:
            pi_s_time, pi_t_time, fst, e_pi_S, e_pi_T, var_pi_s, var_pi_T = calculate_expected_times(k, m, d, N, u)
            theory_pi_s.append(e_pi_S)
            theory_pi_t.append(e_pi_T)
            theory_fst.append(fst)
            theory_pi_s_var.append(var_pi_s)
            theory_pi_t_var.append(var_pi_T)
        
        # Transform x-axis values
        x_fine = x_transform_func(m_fine, k, N)
        x_empirical = x_transform_func(df_k['m'].values, k, N)
        
        # Plot 1: piS
        theory_pi_s_se = np.sqrt(np.array(theory_pi_s_var)) / np.sqrt(num_samples)
        empirical_pi_s_se = df_k['empirical_pi_s_se']
        
        axes[0, col].plot(x_fine, theory_pi_s, color=pi_s_color, linewidth=line_width)
        axes[0, col].fill_between(x_fine, 
                                 np.array(theory_pi_s) - theory_pi_s_se,
                                 np.array(theory_pi_s) + theory_pi_s_se,
                                 alpha=ribbon_alpha, color=pi_s_color, linewidth=0, edgecolor='none')
        axes[0, col].errorbar(x_empirical, df_k['empirical_pi_s'], yerr=empirical_pi_s_se, 
                            fmt='o', color=emp_color, capsize=3, elinewidth=1.0, markersize=markersize)
        if col == 0:
            axes[0, col].set_ylabel('$E[\\pi_s]$', fontsize=font_size, rotation=0, labelpad=10)
        axes[0, col].set_title(f'$k={k}$', fontsize=font_size)
        x_ticks = x_transform_func(np.array(_LOG_M_AXIS_TICKS), k, N)
        _apply_log_xaxis_three_ticks(axes[0, col], x_ticks)
        if log_scale:
            axes[0, col].set_yscale('log')
        axes[0, col].set_ylim(pi_s_ylim)
        if col > 0:
            axes[0, col].set_yticklabels([])
        axes[0, col].grid(True, alpha=0.3)
        
        # Plot 2: piT
        theory_pi_t_se = np.sqrt(np.array(theory_pi_t_var)) / np.sqrt(num_samples)
        empirical_pi_t_se = df_k['empirical_piT_se']
        
        axes[1, col].plot(x_fine, theory_pi_t, color=pi_t_color, linewidth=line_width)
        axes[1, col].fill_between(x_fine, 
                                 np.array(theory_pi_t) - theory_pi_t_se,
                                 np.array(theory_pi_t) + theory_pi_t_se,
                                 alpha=ribbon_alpha, color=pi_t_color, linewidth=0, edgecolor='none')
        axes[1, col].errorbar(x_empirical, df_k['empirical_piT'], yerr=empirical_pi_t_se, 
                            fmt='o', color=emp_color, capsize=3, elinewidth=1.0, markersize=markersize)
        axes[1, col].axhline(y=2*d*u, color='grey', linestyle='--', linewidth=1.5, alpha=0.7)
        
        if col == 0:
            axes[1, col].set_ylabel('$E[\\pi_T]$', fontsize=font_size, rotation=0, labelpad=10)
        _apply_log_xaxis_three_ticks(axes[1, col], x_ticks)
        if log_scale:
            axes[1, col].set_yscale('log')
        axes[1, col].set_ylim(pi_t_ylim)
        if col > 0:
            axes[1, col].set_yticklabels([])
        axes[1, col].grid(True, alpha=0.3)
        
        # Plot 3: FST
        axes[2, col].plot(x_fine, theory_fst, color=fst_color, linewidth=line_width)
        axes[2, col].errorbar(x_empirical, df_k['FST'], yerr=0,
                            fmt='o', color=emp_color, capsize=5)
        axes[2, col].set_xlabel(x_label, fontsize=font_size)
        if col == 0:
            axes[2, col].set_ylabel('$F_{ST}$', fontsize=font_size, rotation=0, labelpad=10)
        _apply_log_xaxis_three_ticks(axes[2, col], x_ticks)
        axes[2, col].set_ylim(fst_ylim)
        if col > 0:
            axes[2, col].set_yticklabels([])
        axes[2, col].grid(True, alpha=0.3)
    
    base_name = plot_file.rsplit('.', 1)[0] if '.' in plot_file else plot_file
    extension = plot_file.rsplit('.', 1)[1] if '.' in plot_file else 'png'
    mode_suffix = "_presentation" if presentation_mode else "_paper"
    config_plot_file = os.path.join(plots_dir, f"{base_name}_all_shared{x_suffix}{mode_suffix}.{extension}")
    
    plt.subplots_adjust(hspace=0.0375, wspace=0.0625)
    plt.savefig(config_plot_file, dpi=600, bbox_inches='tight')
    try:
        pdf_outfile = os.path.join(plots_dir, f"{base_name}_all_shared{x_suffix}{mode_suffix}.pdf")
        plt.savefig(pdf_outfile, bbox_inches='tight')
    except Exception:
        pass
    print(f'plot saved to {config_plot_file}')
    plt.close()

def run_analysis(k_values=[10], d=1000, N=1000, u=1e-6, num_samples=1000,
                output_file='results.csv', plot_file='plot.png', create_plot=True, n_processes=None, num_m_points=19, log_scale=False, presentation_mode=False):
    '''run coalescent analysis over migration rate grid for multiple k values'''
    
    # If in presentation mode, read existing data instead of re-running simulation
    if presentation_mode:
        print(f'Reading existing data from {output_file}')
        df = pd.read_csv(output_file)
    else:
        m_values = np.logspace(-6, -0.01, num_m_points)
        n_processes = min(cpu_count(), len(m_values) * len(k_values)) if n_processes is None else n_processes
        
        print(f'analysis: k={k_values}, d={d}, N={N}, u={u}, {len(m_values)} rates, {n_processes} processes')
        
        # Prepare all combinations of k and m values
        all_args = []
        for k in k_values:
            for m in m_values:
                all_args.append((m, k, d, N, u, num_samples))
        
        start_time = time.time()
        with Pool(n_processes) as pool:
            results = pool.map(process_single_m, all_args)
        
        print(f'completed in {time.time() - start_time:.1f}s')
        
        # Convert to DataFrame
        df = pd.DataFrame(results)
        
        # Save to CSV
        df.to_csv(output_file, index=False)
    
    if create_plot:
        # Create plots directory if it doesn't exist
        plots_dir = 'plots'
        os.makedirs(plots_dir, exist_ok=True)
        
        m_fine = np.logspace(-6, 0, 500)
        
        # Calculate all theoretical values to determine shared y-axis limits
        all_theory_fst = []
        all_empirical_fst = []
        all_theory_pi_s = []
        all_empirical_pi_s = []
        all_theory_pi_t = []
        all_empirical_pi_t = []
        all_theory_pi_s_var = []
        all_theory_pi_t_var = []
        all_empirical_pi_s_var = []
        all_empirical_pi_t_var = []
        
        for k in k_values:
            df_k = df[df['k'] == k]
            for m in m_fine:
                pi_s_time, pi_t_time, fst, e_pi_S, e_pi_T, var_pi_s, var_pi_T = calculate_expected_times(k, m, d, N, u)
                all_theory_pi_s.append(e_pi_S)
                all_theory_pi_t.append(e_pi_T)
                all_theory_fst.append(fst)
                all_theory_pi_s_var.append(var_pi_s)
                all_theory_pi_t_var.append(var_pi_T)
            all_empirical_pi_s.extend(df_k['empirical_pi_s'].values)
            all_empirical_pi_t.extend(df_k['empirical_piT'].values)
            all_empirical_fst.extend(df_k['FST'].values)
            # Calculate empirical standard errors (already in the correct units)
            all_empirical_pi_s_var.extend(df_k['empirical_pi_s_se'].values)
            all_empirical_pi_t_var.extend(df_k['empirical_piT_se'].values)
        
        # Calculate shared y-axis limits
        def calc_ylim(theory_vals, empirical_vals, min_val=0):
            all_vals = theory_vals + empirical_vals
            val_min = min(all_vals)
            val_max = max(all_vals)
            margin = (val_max - val_min) * 0.05
            return (max(min_val, val_min - margin), val_max + margin)
        
        # Convert theoretical variances to standard errors for y-axis calculation
        all_theory_pi_s_se = np.sqrt(np.array(all_theory_pi_s_var)) / np.sqrt(num_samples)
        all_theory_pi_t_se = np.sqrt(np.array(all_theory_pi_t_var)) / np.sqrt(num_samples)
        
        fst_ylim = calc_ylim(all_theory_fst, all_empirical_fst, min_val=0)
        # Ensure piS axis includes highest error bar: use mean + SE for both theory and empirical
        theory_pi_s_upper = (np.array(all_theory_pi_s) + all_theory_pi_s_se).tolist()
        empirical_pi_s_upper = (np.array(all_empirical_pi_s) + np.array(all_empirical_pi_s_var)).tolist()
        # Use zero as lower bound
        pi_s_ylim = (0, max(theory_pi_s_upper + empirical_pi_s_upper) * 1.05)
        # Ensure piT axis includes highest error bar: use mean + SE for both theory and empirical
        theory_pi_t_upper = (np.array(all_theory_pi_t) + all_theory_pi_t_se).tolist()
        empirical_pi_t_upper = (np.array(all_empirical_pi_t) + np.array(all_empirical_pi_t_var)).tolist()
        theory_pi_t_lower = (np.array(all_theory_pi_t) - all_theory_pi_t_se).tolist()
        empirical_pi_t_lower = (np.array(all_empirical_pi_t) - np.array(all_empirical_pi_t_var)).tolist()
        lower_bound_pi_t = min(theory_pi_t_lower + empirical_pi_t_lower)
        upper_bound_pi_t = max(theory_pi_t_upper + empirical_pi_t_upper)
        # Add a small margin around bounds
        pi_t_ylim = (lower_bound_pi_t - 0.05 * (upper_bound_pi_t - lower_bound_pi_t),
                     upper_bound_pi_t + 0.05 * (upper_bound_pi_t - lower_bound_pi_t))
        
        # Plot settings - adjust for paper vs presentation
        if presentation_mode:
            font_size = 14
            line_width = 2.5
            markersize = 5
            ribbon_alpha = 0.3  # More opaque ribbons for projector visibility
            # Keep original colors
            pi_s_color = "#D55E00"  # Orange
            pi_t_color = "#009E73"  # Green
            fst_color = "#CC79A7"   # Pink
            emp_color = "#000000"   # Black
        else:
            font_size = 16
            line_width = 2.0
            markersize = 4
            ribbon_alpha = 0.2
            # Original colors
            pi_s_color = "#D55E00"  # Orange
            pi_t_color = "#009E73"  # Green
            fst_color = "#CC79A7"   # Pink
            emp_color = "#000000"   # Black
        
        plt.rcParams.update({
            'font.size': font_size,
            'axes.titlesize': font_size,
            'axes.labelsize': font_size,
            'xtick.labelsize': font_size,
            'ytick.labelsize': font_size,
            'legend.fontsize': font_size,
            'lines.linewidth': line_width,
            'lines.markersize': markersize,
            'axes.linewidth': 0.8,
            'grid.linewidth': 0.5,
            'pdf.fonttype': 42,
            'ps.fonttype': 42,
        })

        # Create the plot with all axes shared
        # Use smaller per-panel size appropriate for papers and make panels square via box aspect
        unit_size = 3.2  # inches per panel side (2x larger panels)
        fig, axes = plt.subplots(3, len(k_values), figsize=(unit_size*len(k_values), unit_size*3))
        
        # Handle single k value case
        if len(k_values) == 1:
            axes = axes.reshape(-1, 1)

        # Make all panels square and apply consistent tick/spine styling
        for r in range(3):
            for c in range(len(k_values)):
                ax = axes[r, c]
                # Square plotting area
                try:
                    ax.set_box_aspect(1)
                except Exception:
                    ax.set_aspect('equal', adjustable='box')
                ax.tick_params(width=0.8, length=3)
                # Hide x tick labels on top two rows
                if r < 2:
                    ax.tick_params(labelbottom=False)
                for spine in ax.spines.values():
                    spine.set_linewidth(0.8)
        
        for col, k in enumerate(k_values):
            # Filter data for current k value
            df_k = df[df['k'] == k]
            
            # Calculate theoretical values for current k
            theory_pi_s = []
            theory_pi_t = []
            theory_fst = []
            theory_pi_s_var = []
            theory_pi_t_var = []
            for m in m_fine:
                pi_s_time, pi_t_time, fst, e_pi_S, e_pi_T, var_pi_s, var_pi_T = calculate_expected_times(k, m, d, N, u)
                theory_pi_s.append(e_pi_S)
                theory_pi_t.append(e_pi_T)
                theory_fst.append(fst)
                theory_pi_s_var.append(var_pi_s)
                theory_pi_t_var.append(var_pi_T)
            
            # Plot 1: piS
            # Convert variances to standard errors for ribbons
            theory_pi_s_se = np.sqrt(np.array(theory_pi_s_var)) / np.sqrt(num_samples)
            empirical_pi_s_se = df_k['empirical_pi_s_se']
            
            # Plot theoretical mean as line
            axes[0, col].plot(m_fine, theory_pi_s, color=pi_s_color, linewidth=line_width)
            
            # Plot theoretical standard error as ribbon
            axes[0, col].fill_between(m_fine, 
                                     np.array(theory_pi_s) - theory_pi_s_se,
                                     np.array(theory_pi_s) + theory_pi_s_se,
                                     alpha=ribbon_alpha, color=pi_s_color, linewidth=0, edgecolor='none')
            
            # Plot empirical mean as points with empirical standard error as error bars
            axes[0, col].errorbar(df_k['m'], df_k['empirical_pi_s'], yerr=empirical_pi_s_se, 
                                fmt='o', color=emp_color, capsize=3, elinewidth=1.0, markersize=markersize)
            if col == 0:  # Only label y-axis for k=1
                axes[0, col].set_ylabel('$E[\\pi_S]$', fontsize=font_size, rotation=0, labelpad=27)
            axes[0, col].set_title(f'$k={k}$', fontsize=font_size)
            _apply_log_xaxis_three_ticks(axes[0, col])
            if log_scale:
                axes[0, col].set_yscale('log')
            axes[0, col].set_ylim(pi_s_ylim)
            if col > 0:  # Hide y-axis tick labels for non-k=1 columns in shared plots
                axes[0, col].set_yticklabels([])
            axes[0, col].grid(True, alpha=0.3)
            
            # Plot 2: piT
            # Convert variances to standard errors for ribbons
            theory_pi_t_se = np.sqrt(np.array(theory_pi_t_var)) / np.sqrt(num_samples)
            empirical_pi_t_se = df_k['empirical_piT_se']
            
            # Plot theoretical mean as line
            axes[1, col].plot(m_fine, theory_pi_t, color=pi_t_color, linewidth=line_width)
            
            # Plot theoretical standard error as ribbon
            axes[1, col].fill_between(m_fine, 
                                     np.array(theory_pi_t) - theory_pi_t_se,
                                     np.array(theory_pi_t) + theory_pi_t_se,
                                     alpha=ribbon_alpha, color=pi_t_color, linewidth=0, edgecolor='none')
            
            # Plot empirical mean as points with empirical standard error as error bars
            axes[1, col].errorbar(df_k['m'], df_k['empirical_piT'], yerr=empirical_pi_t_se, 
                                fmt='o', color=emp_color, capsize=3, elinewidth=1.0, markersize=markersize)
            
            # # Add grey dashed line at y = 2*d*u
            # axes[1, col].axhline(y=2*d*u, color='grey', linestyle='--', linewidth=1.5, alpha=0.7)
            
            if col == 0:  # Only label y-axis for k=1
                axes[1, col].set_ylabel('$E[\\pi_T]$', fontsize=font_size, rotation=0, labelpad=27)
            _apply_log_xaxis_three_ticks(axes[1, col])
            if log_scale:
                axes[1, col].set_yscale('log')
            axes[1, col].set_ylim(pi_t_ylim)
            if col > 0:  # Hide y-axis tick labels for non-k=1 columns in shared plots
                axes[1, col].set_yticklabels([])
            axes[1, col].grid(True, alpha=0.3)
            
            # Plot 3: FST
            axes[2, col].plot(m_fine, theory_fst, color=fst_color, linewidth=line_width)
            axes[2, col].errorbar(df_k['m'], df_k['FST'], yerr=0,
                                fmt='o', color=emp_color, capsize=5)
            axes[2, col].set_xlabel('$m$', fontsize=font_size)
            if col == 0:  # Only label y-axis for k=1
                axes[2, col].set_ylabel('$F_{ST}$', fontsize=font_size, rotation=0, labelpad=20)
            _apply_log_xaxis_three_ticks(axes[2, col])
            axes[2, col].set_ylim(fst_ylim)
            if col > 0:  # Hide y-axis tick labels for non-k=1 columns in shared plots
                axes[2, col].set_yticklabels([])
            axes[2, col].grid(True, alpha=0.3)
        
        # Generate filename for the all-shared configuration
        base_name = plot_file.rsplit('.', 1)[0] if '.' in plot_file else plot_file
        extension = plot_file.rsplit('.', 1)[1] if '.' in plot_file else 'png'
        mode_suffix = "_presentation" if presentation_mode else "_paper"
        config_plot_file = os.path.join(plots_dir, f"{base_name}_all_shared{mode_suffix}.{extension}")
        
        # Reduce whitespace between subplots (25% of prior spacing)
        plt.subplots_adjust(hspace=0.0375, wspace=0.0625)
        plt.savefig(config_plot_file, dpi=600, bbox_inches='tight')
        # Also save as PDF alongside the raster image
        try:
            pdf_outfile = os.path.join(plots_dir, f"{base_name}_all_shared{mode_suffix}.pdf")
            plt.savefig(pdf_outfile, bbox_inches='tight')
        except Exception:
            pass
        print(f'plot saved to {config_plot_file}')
        plt.close()  # Close figure to free memory
    
    return df

def main():
    parser = argparse.ArgumentParser(description='backward coalescent simulation')
    parser.add_argument('--k', type=int, nargs='+', default=[1, 2, 5, 100], 
                       help='k values to analyze (can specify multiple, e.g., --k 5 10 20)')
    parser.add_argument('--d', type=int, default=100)
    parser.add_argument('--N', type=int, default=1e10)
    parser.add_argument('--u', type=float, default=1e-6)
    parser.add_argument('--num-samples', type=int, default=200000)
    parser.add_argument('--output-file', default='coalescent_results.csv')
    parser.add_argument('--plot-file', default='coalescent_plot.png')
    parser.add_argument('--no-plot', action='store_true')
    parser.add_argument('--processes', type=int, default=None)
    parser.add_argument('--num-m-points', type=int, default=20)
    parser.add_argument('--log-scale', action='store_true', help='Use log scale for piS and piT y-axes')
    
    args = parser.parse_args()
    
    should_plot = not args.no_plot
    
    # Run analysis once - this will do the simulation and create the paper plot
    df = run_analysis(args.k, args.d, args.N, args.u, args.num_samples,
                     args.output_file, args.plot_file, should_plot, 
                     args.processes, args.num_m_points, args.log_scale, 
                     presentation_mode=False)
    
    # If plotting is enabled, also create presentation version
    if should_plot:
        run_analysis(args.k, args.d, args.N, args.u, args.num_samples,
                    args.output_file, args.plot_file, True, 
                    args.processes, args.num_m_points, args.log_scale, 
                    presentation_mode=True)

if __name__ == '__main__':
    main()