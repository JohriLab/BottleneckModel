import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

u = 1e-6
d = 100
N = 1000
k_values = [1, 2, 5, 10]

def calculate_expected_times(k, m, d, N):
    alpha = 2*m*(1-m)
    beta = (1 - m)**2 + 2*m*(1-m)*(d-2)/(d-1) + m**2*(d-2)/(d-1)

    Q = np.array([[(1-1/k)*(1-1/N)*(1-alpha), (1-1/k)*(1-1/N) * alpha],
                [1-beta*(1-1/d), beta*(1-1/d)]])

    I = np.eye(2)
    N_matrix = np.linalg.inv(I - Q)

    # Calculate expected absorption times
    expected_times = N_matrix @ np.ones((2, 1))

    return expected_times

# Define a range of m values on a log scale
m_values = np.logspace(-6, 0, 1000)

# Calculate theoretical expectations for each k value
all_results = []
for k in k_values:
    results = []
    for m_val in m_values:
        times = calculate_expected_times(k, m_val, d, N)
        results.append({
            'k': k,
            'm': m_val,
            'expected_time_within': times[0, 0],
            'expected_time_between': times[1, 0],
            'expected_time_total': times[0, 0]/d + times[1, 0]*(1-1/d),
            'expected_piS': 2*u*times[0, 0],
            'expected_piB': 2*u*times[1, 0],
            'expected_piT': 2*u*times[0, 0]/d + 2*u*times[1, 0]*(1-1/d),
            'FST': 1 - times[0, 0]/(times[0, 0]/d + times[1, 0]*(1-1/d))
        })
    all_results.extend(results)

df = pd.DataFrame(all_results)

# Load empirical data
empirical_data = pd.read_csv('pi_data_plotting_backwards.csv')

# Create faceted plot with 4 rows and 4 columns (one for each k value)
fig, axes = plt.subplots(4, 4, figsize=(20, 16))

# Define plot types and their properties
plot_types = [
    {'name': 'FST', 'column': 'FST', 'color': 'red'},
    {'name': 'πS (within-deme)', 'column': 'expected_piS', 'color': 'blue'},
    {'name': 'πB (between-deme)', 'column': 'expected_piB', 'color': 'orange'},
    {'name': 'πT (total)', 'column': 'expected_piT', 'color': 'green'}
]

# Create plots for each k value and each plot type
for col_idx, k_val in enumerate(k_values):
    # Filter data for current k value
    k_theoretical = df[df['k'] == k_val]
    k_empirical = empirical_data[empirical_data['k'] == k_val]
    
    for row_idx, plot_type in enumerate(plot_types):
        ax = axes[row_idx, col_idx]
        
        # Plot theoretical curve
        ax.plot(k_theoretical['m'], k_theoretical[plot_type['column']], 
                label=f'{plot_type["name"]} (theoretical)', 
                color=plot_type['color'], linewidth=2)
        
        # For πT plot, add empirical data
        if plot_type['name'] == 'πT (total)':
            ax.errorbar(k_empirical['m'], k_empirical['pi'], 
                       yerr=k_empirical['se_pi'],
                       label=f'{plot_type["name"]} (empirical)', 
                       color='black', marker='o', markersize=4, 
                       capsize=3, linestyle='None')
        
        # Formatting
        ax.set_xscale('log')
        ax.set_ylabel(plot_type['name'])
        ax.grid(True, alpha=0.3)
        
        # Add k value as title for top row
        if row_idx == 0:
            ax.set_title(f'k = {k_val}', fontsize=12, fontweight='bold')
        
        # Add x-axis label only for bottom row
        if row_idx == 3:
            ax.set_xlabel('Migration rate (m)')

plt.tight_layout()
plt.savefig('expectations_plot.png', dpi=300, bbox_inches='tight')
plt.close()

