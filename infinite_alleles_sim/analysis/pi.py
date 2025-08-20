import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Read in pi and sfs data
pi_df = pd.read_csv('analysis/pi.csv')
sfs_df = pd.read_csv('analysis/sfs.csv')

# If sfs_df has no column names, add them
if sfs_df.shape[1] == 4:
    sfs_df.columns = ['count', 'k', 'm', 'replicate']

sfs = (sfs_df['count'].value_counts().sort_index() / len(sfs_df)).to_frame('frequency')
    

# If pi_df has no column names, add them
if pi_df.shape[1] == 5:
    pi_df.columns = ['piS', 'piT', 'k', 'm', 'replicate']

output_df = pd.DataFrame()

# Calculate and print mean piT, piS, FST for each k and m combination
grouped = pi_df.groupby(['k', 'm'])
for (k_val, m_val), group in grouped:
    piT_mean = group['piT'].mean()
    piS_mean = group['piS'].mean()
    FST_mean = 1 - piS_mean / piT_mean
    print(f"k = {k_val}, m = {m_val}")
    print(f"  piT_mean: {piT_mean}")
    print(f"  piS_mean: {piS_mean}")
    print(f"  FST_mean: {FST_mean}")
    piT_se = group['piT'].std(ddof=1) / (len(group) ** 0.5)
    piS_se = group['piS'].std(ddof=1) / (len(group) ** 0.5)
    # Calculate FST standard error only for entries where piT is not 0
    valid = group['piT'] != 0
    if valid.sum() > 1:
        FST_vals = 1 - group.loc[valid, 'piS'] / group.loc[valid, 'piT']
        FST_se = FST_vals.std(ddof=1) / (len(FST_vals) ** 0.5)
    else:
        FST_se = np.nan
    num_replicates = len(group)
    new_row = pd.DataFrame({
        'm': [m_val],
        'k': [k_val],
        'piT': [piT_mean],
        'piT_se': [piT_se],
        'piS': [piS_mean],
        'piS_se': [piS_se],
        'FST': [FST_mean],
        'FST_se': [FST_se],
        'num_replicates': [num_replicates]
    })
    output_df = pd.concat([output_df, new_row], ignore_index=True)

output_df.to_csv('analysis/pi_data_plotting.csv', index=False)

# # Get unique k and m values for faceting
# k_values = sorted(sfs_df['k'].unique())
# m_values = sorted(sfs_df['m'].unique())
# n_k = len(k_values)
# n_m = len(m_values)

# print(f"Found {n_k} k values: {k_values}")
# print(f"Found {n_m} m values: {m_values}")

# # Create faceted SFS plots with k on x-axis and m on y-axis
# fig, axes = plt.subplots(n_m, n_k, figsize=(4*n_k, 4*n_m), sharey=True, sharex=True)

# # Handle case where there's only one k or m value
# if n_k == 1:
#     axes = axes.reshape(-1, 1)
# if n_m == 1:
#     axes = axes.reshape(1, -1)

# # Unbinned SFS plots faceted by k and m
# for i, m_val in enumerate(m_values):
#     for j, k_val in enumerate(k_values):
#         ax = axes[i, j]
        
#         # Filter sfs_df for this k and m combination and exclude the max count
#         subset = sfs_df[(sfs_df['k'] == k_val) & (sfs_df['m'] == m_val) & (sfs_df['count'] < sfs_df['count'].max())]
        
#         if not subset.empty:
#             # Unbinned SFS: plot the raw count frequencies
#             sfs_unbinned = (subset['count'].value_counts().sort_index()).to_frame('frequency')
#             sfs_unbinned.index.name = 'count'
            
#             ax.bar(sfs_unbinned.index.astype(str), sfs_unbinned['frequency'], color='skyblue', edgecolor='black')
#             ax.set_title(f"k={k_val}, m={m_val}")
            
#             # Only add x-label for bottom row
#             if i == n_m - 1:
#                 ax.set_xlabel('Count')
#             # Only add y-label for leftmost column
#             if j == 0:
#                 ax.set_ylabel('Frequency')
                
#             ax.tick_params(axis='x', rotation=45)
#         else:
#             ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
#             ax.set_title(f"k={k_val}, m={m_val}")

# plt.tight_layout()
# plt.savefig('analysis/sfs_faceted_by_k_and_m.png', dpi=300, bbox_inches='tight')

# # Create binned SFS plots with the same faceting
# fig_binned, axes_binned = plt.subplots(n_m, n_k, figsize=(4*n_k, 4*n_m), sharey=True, sharex=True)

# # Handle case where there's only one k or m value
# if n_k == 1:
#     axes_binned = axes_binned.reshape(-1, 1)
# if n_m == 1:
#     axes_binned = axes_binned.reshape(1, -1)

# for i, m_val in enumerate(m_values):
#     for j, k_val in enumerate(k_values):
#         ax_binned = axes_binned[i, j]
        
#         # Filter sfs_df for this k and m combination and exclude the max count
#         subset = sfs_df[(sfs_df['k'] == k_val) & (sfs_df['m'] == m_val) & (sfs_df['count'] < sfs_df['count'].max())]
        
#         if not subset.empty:
#             # Bin the counts into 20 bins
#             min_count = subset['count'].min()
#             max_count = subset['count'].max()
#             # If all counts are the same, just make one bin
#             if min_count == max_count:
#                 bins = [min_count, max_count + 1]
#                 labels = [str(min_count)]
#             else:
#                 bins = np.linspace(min_count, max_count, 21)
#                 labels = [f"{int(bins[j])}-{int(bins[j+1])-1}" for j in range(len(bins)-1)]
#             # Bin the data
#             binned = pd.cut(subset['count'], bins=bins, labels=labels, include_lowest=True, right=False)
#             sfs_binned = binned.value_counts().sort_index().to_frame('frequency')
#             sfs_binned.index.name = 'count_bin'
            
#             ax_binned.bar(sfs_binned.index.astype(str), sfs_binned['frequency'], color='salmon', edgecolor='black')
#             ax_binned.set_title(f"k={k_val}, m={m_val}")
            
#             # Only add x-label for bottom row
#             if i == n_m - 1:
#                 ax_binned.set_xlabel('Count Bin')
#             # Only add y-label for leftmost column
#             if j == 0:
#                 ax_binned.set_ylabel('Frequency')
                
#             ax_binned.tick_params(axis='x', rotation=45)
#         else:
#             ax_binned.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax_binned.transAxes)
#             ax_binned.set_title(f"k={k_val}, m={m_val}")

# plt.tight_layout()
# plt.savefig('analysis/sfs_binned_faceted_by_k_and_m.png', dpi=300, bbox_inches='tight')

# # Create log-scaled binned SFS plots with the same faceting
# fig_binned_log, axes_binned_log = plt.subplots(n_m, n_k, figsize=(4*n_k, 4*n_m), sharey=True, sharex=True)

# # Handle case where there's only one k or m value
# if n_k == 1:
#     axes_binned_log = axes_binned_log.reshape(-1, 1)
# if n_m == 1:
#     axes_binned_log = axes_binned_log.reshape(1, -1)

# for i, m_val in enumerate(m_values):
#     for j, k_val in enumerate(k_values):
#         ax_binned_log = axes_binned_log[i, j]
        
#         # Filter sfs_df for this k and m combination and exclude the max count
#         subset = sfs_df[(sfs_df['k'] == k_val) & (sfs_df['m'] == m_val) & (sfs_df['count'] < sfs_df['count'].max())]
        
#         if not subset.empty:
#             # Bin the counts into 20 bins
#             min_count = subset['count'].min()
#             max_count = subset['count'].max()
#             # If all counts are the same, just make one bin
#             if min_count == max_count:
#                 bins = [min_count, max_count + 1]
#                 labels = [str(min_count)]
#             else:
#                 bins = np.linspace(min_count, max_count, 21)
#                 labels = [f"{int(bins[j])}-{int(bins[j+1])-1}" for j in range(len(bins)-1)]
#             # Bin the data
#             binned = pd.cut(subset['count'], bins=bins, labels=labels, include_lowest=True, right=False)
#             sfs_binned = binned.value_counts().sort_index().to_frame('frequency')
#             sfs_binned.index.name = 'count_bin'
            
#             # Add small constant to avoid log(0) issues
#             frequencies = sfs_binned['frequency'] + 1e-10
            
#             ax_binned_log.bar(sfs_binned.index.astype(str), frequencies, color='lightcoral', edgecolor='black')
#             ax_binned_log.set_yscale('log')
#             ax_binned_log.set_title(f"k={k_val}, m={m_val}")
            
#             # Only add x-label for bottom row
#             if i == n_m - 1:
#                 ax_binned_log.set_xlabel('Count Bin')
#             # Only add y-label for leftmost column
#             if j == 0:
#                 ax_binned_log.set_ylabel('Frequency (log scale)')
                
#             ax_binned_log.tick_params(axis='x', rotation=45)
#         else:
#             ax_binned_log.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax_binned_log.transAxes)
#             ax_binned_log.set_title(f"k={k_val}, m={m_val}")

# plt.tight_layout()
# plt.savefig('analysis/sfs_binned_faceted_by_k_and_m_log_scale.png', dpi=300, bbox_inches='tight')

# import os

# # Ensure the 'sfs' directory exists
# os.makedirs('analysis/sfs', exist_ok=True)

# # Save individual SFS files for each k and m combination
# for k_val in k_values:
#     for m_val in m_values:
#         # Subset for this k and m combination and exclude the max count
#         subset = sfs_df[(sfs_df['k'] == k_val) & (sfs_df['m'] == m_val) & (sfs_df['count'] < sfs_df['count'].max())]
        
#         if not subset.empty:
#             # Unbinned SFS
#             sfs_unbinned = (subset['count'].value_counts().sort_index()).to_frame('occurrences')
#             sfs_unbinned.index.name = 'count'

#             # Binned SFS
#             min_count = subset['count'].min()
#             max_count = subset['count'].max()
#             if min_count == max_count:
#                 bins = [min_count, max_count + 1]
#                 labels = [str(min_count)]
#             else:
#                 bins = np.linspace(min_count, max_count, 21)
#                 labels = [f"{int(bins[j])}-{int(bins[j+1])-1}" for j in range(len(bins)-1)]
#             binned = pd.cut(subset['count'], bins=bins, labels=labels, include_lowest=True, right=False)
#             sfs_binned = (binned.value_counts().sort_index()).to_frame('occurrences')
#             sfs_binned.index.name = 'count_bin'
#         else:
#             sfs_unbinned = pd.DataFrame({'occurrences': []})
#             sfs_binned = pd.DataFrame({'occurrences': []})

#         # Write unbinned SFS to CSV file in 'sfs' folder, filename includes k and m values
#         k_str = f"{k_val:.6g}".replace('.', 'p')
#         m_str = f"{m_val:.6g}".replace('.', 'p')
#         out_path_unbinned = f"sfs/sfs_k_{k_str}_m_{m_str}.csv"
#         sfs_unbinned.reset_index().to_csv(out_path_unbinned, index=False)

#         # Write binned SFS to CSV file in 'sfs' folder, filename includes k and m values
#         out_path_binned = f"sfs/sfs_binned_k_{k_str}_m_{m_str}.csv"
#         sfs_binned.reset_index().to_csv(out_path_binned, index=False)

# print("Analysis complete! Check the generated plots and CSV files.")


