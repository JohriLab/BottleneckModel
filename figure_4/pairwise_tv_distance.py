import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ============================================================================
# PARAMETERS - Edit these values to change plot settings
# ============================================================================
# Define sets of (d, k) values to compute total variation distance for
d_values = [2, 5, 10, 50, 100, 500, 1000, 5000, 10000]  # d values to test
k_values = [2, 5, 10, 50, 100, 500, 1000, 5000, 10000]      # k values to test
T_max = 1000000  # Maximum time T to compute distributions for
m = 0.0002  # migration rate parameter for full_model_distribution
# d_values = [1000, 5000, 9500, 9900]  # d values to test
# k_values = [100, 500, 5000, 9000]      # k values to test
# T_max = 1000000  # Maximum time T to compute distributions for
# m = 0.0002  # migration rate parameter for full_model_distribution
N = 1e10  # deme size for full_model_distribution
# ============================================================================

def full_model_distribution(m, N, d, k, T_max):
    """
    p_times[t - 1, i] is P(coalesce at time t | start in state i), t = 1, ..., T_max.
    """

    a = 1 - 1/k
    b = 1 - 1/N
    c = 1/d
    alpha = 2*m*(1-m)

    Q = np.array([[a*b*(1-alpha), a*b*alpha], [c*b, 1-c]])

    n = Q.shape[0]
    powers = np.empty((T_max, 2, 2), dtype=Q.dtype)
    powers[0] = np.eye(2)
    for t in range(1, T_max):
        powers[t] = powers[t - 1] @ Q

    I = np.eye(n)
    ones = np.ones(n, dtype=Q.dtype)
    iq_ones = (I-Q) @ ones
    # p_times = np.einsum("tij,j->ti", powers, iq_ones)
    # p_times = powers @ iq_ones

    mu = np.array([[c, 1-c]])
    p_times = mu @ powers @iq_ones

    q_T = powers[-1] @ Q
    tail_mass = mu @q_T @ ones

    return p_times, tail_mass[0]

def compute_expectation(m, N, d, k):
    """
    returns the expectation of T_T for fixed n, N, d, k
    """

    a = 1 - 1/k
    b = 1 - 1/N
    c = 1/d
    alpha = 2*m*(1-m)

    Q = np.array([[a*b*(1-alpha), a*b*alpha], [c*b, 1-c]])

    n = Q.shape[0]
    I = np.eye(n)
    F = np.linalg.inv(I-Q)

    mu = np.array([[c, 1-c]])
    ones = np.ones(n, dtype=Q.dtype)
    expectation = mu @ F @ ones

    return expectation[0]

def wf_distribution(Ne, T_max):
    """
    Computes a geometric distribution with p = Ne (i.e. the distribution of time to coalescence under the WF model).
    """
    t = np.arange(0, T_max)

    p = 1/Ne
    p_times = p * (1-p)**(t)
    return p_times

def compute_total_variation_distance(d, k, T_max, m, N):
    """
    Compute total variation distance between model and WF distributions.
    
    Total variation distance: TV(P, Q) = 0.5 * sum(|P(i) - Q(i)|)
    
    Parameters:
    -----------
    d : int
        Number of demes
    k : int
        k value
    T_max : int
        Maximum time to compute distributions for
    m : float
        Migration parameter passed to full_model_distribution
    N : float
        Deme size passed to full_model_distribution
    
    Returns:
    --------
    tv_distance : float
        Total variation distance between the two distributions
    """
    model_probs, _tail = full_model_distribution(m, N, d, k, T_max)
    model_probs = np.asarray(model_probs, dtype=float).reshape(-1)

    Ne = compute_expectation(m, N, d, k)
    WF_probs = np.asarray(wf_distribution(Ne, T_max), dtype=float).reshape(-1)

    if model_probs.shape[0] != WF_probs.shape[0]:
        raise ValueError(
            f"PMF length mismatch: model {model_probs.shape[0]} vs WF {WF_probs.shape[0]}"
        )

    # Normalize to get probability mass functions on t = 1, ..., T_max (same as figure_3)
    model_probs = model_probs / np.sum(model_probs)
    WF_probs = WF_probs / np.sum(WF_probs)
    
    # Compute total variation distance
    tv_distance = 0.5 * np.sum(np.abs(model_probs - WF_probs))
    
    return tv_distance

def plot_tv_distance_heatmap(pivot_table, output_file='tv_distance_heatmap.png'):
    """
    Create a heatmap visualization of the total variation distance pivot table.
    
    Parameters:
    -----------
    pivot_table : pandas.DataFrame
        Pivot table with d as index, k as columns, and TV Distance as values
    output_file : str
        Output filename for the heatmap
    """
    plt.rcParams.update({
        'font.size': 9,
        'axes.titlesize': 9,
        'axes.labelsize': 9,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'lines.linewidth': 2.0,
        'axes.linewidth': 0.8,
        'grid.linewidth': 0.5,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })
    
    # Convert pivot table to numpy array for plotting
    data = pivot_table.values
    
    # Create figure with appropriate size for word document
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Create heatmap using imshow
    # Put the first row (smallest d) at the bottom so the largest d is at the top.
    im = ax.imshow(
        data,
        cmap='YlOrRd',
        aspect='auto',
        interpolation='nearest',
        origin='lower',
    )
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(pivot_table.columns)))
    ax.set_yticks(np.arange(len(pivot_table.index)))
    ax.set_xticklabels(pivot_table.columns)
    ax.set_yticklabels(pivot_table.index)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Total Variation Distance', rotation=270, labelpad=15)
    
    # Add text annotations for each cell
    for i in range(len(pivot_table.index)):
        for j in range(len(pivot_table.columns)):
            value = data[i, j]
            # Use white text for dark backgrounds, black for light backgrounds
            text_color = 'white' if value > np.nanmax(data) * 0.5 else 'black'
            ax.text(j, i, f'{value:.4f}', 
                   ha='center', va='center', 
                   color=text_color, fontsize=8)
    
    # Set labels (italicize d and k)
    ax.set_xlabel(r'$k$')
    ax.set_ylabel(r'$d$', rotation=0, labelpad=10)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    try:
        pdf_file = output_file.rsplit('.', 1)[0] + '.pdf'
        plt.savefig(pdf_file, bbox_inches='tight')
        print(f'Heatmap saved to {pdf_file}')
    except Exception:
        pass
    print(f'Heatmap saved to {output_file}')
    plt.close()

def create_tv_distance_table(d_values, k_values, T_max, m, N):
    """
    Create a table of total variation distances for all (d, k) pairs.
    
    Parameters:
    -----------
    d_values : list
        List of d values to test
    k_values : list
        List of k values to test
    T_max : int
        Maximum time to compute distributions for
    m : float
        Migration parameter for full_model_distribution
    N : float
        Deme size for full_model_distribution
    
    Returns:
    --------
    df : pandas.DataFrame
        DataFrame with d, k, and TV distance columns
    """
    results = []
    
    print(f'Computing total variation distances for {len(d_values)} d values and {len(k_values)} k values...')
    print(f'Total combinations: {len(d_values) * len(k_values)}')
    
    for d in d_values:
        for k in k_values:
            print(f'  Computing TV distance for d={d}, k={k}...')
            tv_dist = compute_total_variation_distance(d, k, T_max, m, N)
            results.append({'d': d, 'k': k, 'TV Distance': tv_dist})
    
    df = pd.DataFrame(results)
    return df

def main():
    print('=' * 60)
    print('Total Variation Distance Calculation')
    print('=' * 60)
    print(f'd values: {d_values}')
    print(f'k values: {k_values}')
    print(f'T_max: {T_max}')
    print()
    
    # Compute total variation distances for all (d, k) pairs
    df = create_tv_distance_table(d_values, k_values, T_max, m, N)
    
    # Create a pivot table for better visualization
    pivot_table = df.pivot(index='d', columns='k', values='TV Distance')
    
    print('\n' + '=' * 60)
    print('Total Variation Distance Table')
    print('=' * 60)
    print('\nFull table:')
    print(df.to_string(index=False))
    
    print('\n\nPivot table (d rows, k columns):')
    print(pivot_table.to_string())

    # Create and save heatmap visualization
    print('\nCreating heatmap visualization...')
    plot_tv_distance_heatmap(pivot_table, 'tv_distance_heatmap.png')

if __name__ == '__main__':
    main()