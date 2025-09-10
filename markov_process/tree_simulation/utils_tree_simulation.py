import argparse
import random
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
from itertools import combinations

def simulate_tree(k, m, N, d, num_samples):
    '''simulate a coalescent tree with migration, deme drift, and coalescence'''
    '''
    Args:
        k: number of lineages
        m: migration rate
        N: population size
        d: number of demes
        num_samples: number of samples

    Returns:
        edges_df: DataFrame with edge list
    '''

    p_coal = 1 - (1-1/k)*(1-1/N)
    
    sample_id = 0
    samples = []
    for i in range(num_samples):
        samples.append({
            'id': sample_id,
            'deme': random.randint(1, d),
            'time': 0,
            'leaves_below': 1
        })
        sample_id += 1

    edges = []
    current_time = 0

    while len(samples) > 1:
        current_time += 1

        ### Migration: ###
        for sample in samples:
            if np.random.random() < m:
                current_deme = sample['deme']
                new_deme = random.randint(1, d-1)
                if new_deme >= current_deme:
                    new_deme += 1
                sample['deme'] = new_deme

        ### Deme Drift: ###
        unique_demes = list(set(sample['deme'] for sample in samples))
        deme_parent_mapping = {}
        for deme in unique_demes:
            parent_deme = random.randint(1, d)
            deme_parent_mapping[deme] = parent_deme
        
        for sample in samples:
            sample['deme'] = deme_parent_mapping[sample['deme']]

        ### Coalescence: ###
        deme_groups = defaultdict(list)
        for sample in samples:
            deme_groups[sample['deme']].append(sample)
            
        new_samples = []
        for deme, deme_samples in deme_groups.items():
            if len(deme_samples) > 1:
                coalescing_pairs = []
                
                for pair in combinations(deme_samples, 2):
                    if np.random.random() < p_coal:
                        coalescing_pairs.append(pair)
                
                if coalescing_pairs:
                    # Use the simple NetworkX approach just for component finding
                    G = nx.Graph()
                    for pair in coalescing_pairs:
                        G.add_edge(pair[0]['id'], pair[1]['id'])
                    
                    components = list(nx.connected_components(G))
                    
                    # Convert back to sample objects
                    sample_components = []
                    for component in components:
                        sample_component = []
                        for node_id in component:
                            for sample in deme_samples:
                                if sample['id'] == node_id:
                                    sample_component.append(sample)
                                    break
                        sample_components.append(sample_component)
                    
                    # Create coalescent nodes for each component
                    coalesced_lineages = set()
                    for component in sample_components:
                        total_leaves = sum(lineage['leaves_below'] for lineage in component)
                        
                        # Create coalescent node
                        coal_node = {
                            'id': sample_id,
                            'deme': deme,
                            'time': current_time,
                            'leaves_below': total_leaves
                        }
                        sample_id += 1
                        new_samples.append(coal_node)
                        
                        # Add edges from coalescent node to children
                        for child in component:
                            branch_length = current_time - child['time']
                            edges.append({
                                'parent': coal_node['id'],
                                'child': child['id'],
                                'branch_length': branch_length,
                                'coalescence_time': current_time,
                                'deme': deme,
                                'leaves_below': child['leaves_below']
                            })
                            coalesced_lineages.add(child['id'])
                    
                    # Add remaining lineages that didn't coalesce
                    remaining = [s for s in deme_samples if s['id'] not in coalesced_lineages]
                    new_samples.extend(remaining)
                else:
                    new_samples.extend(deme_samples)
            else:
                new_samples.extend(deme_samples)
        
        samples = new_samples

    return pd.DataFrame(edges)

def plot_tree_tskit_style(edges_df, title="Coalescent Tree", figsize=(12, 8), save_file=None, ax=None):
    """
    Plot tree in tskit-style rectangular format with clean lines.
    Args:
        edges_df: DataFrame with edge list
        title: Title of the plot
        figsize: Size of the figure
        save_file: Optional filename to save the plot
        ax: Optional Axes object to plot on

    Returns:
        fig: Figure object
        ax: Axes object
    """
    # Create a directed graph
    G = nx.DiGraph()
    
    # Add edges to the graph
    for _, row in edges_df.iterrows():
        G.add_edge(row['parent'], row['child'], 
                  branch_length=row['branch_length'],
                  coalescence_time=row['coalescence_time'])
    
    # Find the root
    root = None
    for node in G.nodes():
        if G.in_degree(node) == 0:
            root = node
            break
    
    if root is None:
        print("Error: No root found in the tree")
        return
    
    # Get all leaf nodes (samples) and sort them
    leaves = [node for node in G.nodes() if G.out_degree(node) == 0]
    leaves.sort()
    
    # Calculate node positions
    positions = {}
    
    # Position leaves at the top (time 0) - FLIPPED
    for i, leaf in enumerate(leaves):
        positions[leaf] = (i, 0)
    
    # Calculate internal node positions
    def get_node_time(node):
        """Get the coalescence time for a node."""
        for _, row in edges_df.iterrows():
            if row['parent'] == node:
                return row['coalescence_time']
        return 0
    
    def calculate_internal_positions(node):
        """Calculate positions for internal nodes."""
        if node in positions:
            return positions[node]
        
        children = list(G.neighbors(node))
        if not children:
            return positions[node]
        
        # Recursively calculate positions for children
        child_positions = [calculate_internal_positions(child) for child in children]
        
        # Position internal node at the midpoint of its children
        x = sum(pos[0] for pos in child_positions) / len(child_positions)
        y = get_node_time(node)
        
        positions[node] = (x, y)
        return positions[node]
    
    # Calculate all positions
    calculate_internal_positions(root)
    
    # Create the plot
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    # Draw the tree with rectangular style
    for edge in G.edges():
        parent, child = edge
        parent_pos = positions[parent]
        child_pos = positions[child]
        
        # Draw vertical line from child to parent's time level
        ax.plot([child_pos[0], child_pos[0]], [child_pos[1], parent_pos[1]], 
                'k-', linewidth=2, solid_capstyle='round')
        
        # Draw horizontal line from child to parent
        ax.plot([child_pos[0], parent_pos[0]], [parent_pos[1], parent_pos[1]], 
                'k-', linewidth=2, solid_capstyle='round')
    
    # Draw nodes
    for node in G.nodes():
        pos = positions[node]
        if G.out_degree(node) > 0:  # Internal node (coalescent event)
            # Count how many lineages coalesced at this node
            num_lineages = G.out_degree(node)
            
            # Color red if more than 2 lineages coalesced
            if num_lineages > 2:
                ax.plot(pos[0], pos[1], 'ro', markersize=8, markeredgecolor='black', markeredgewidth=1)
            else:
                ax.plot(pos[0], pos[1], 'ko', markersize=8, markeredgecolor='black', markeredgewidth=1)
            
            ax.text(pos[0], pos[1] - 0.3, f'{node}', fontsize=11, ha='center', va='top', fontweight='bold')
        else:  # Leaf node (sample)
            ax.plot(pos[0], pos[1], 'ks', markersize=8, markeredgecolor='black', markeredgewidth=1)
            ax.text(pos[0], pos[1] + 0.3, f'{node}', fontsize=11, ha='center', va='bottom')
    
    # Customize the plot
    ax.set_xlabel('Sample Index', fontsize=12)
    ax.set_ylabel('Time (generations ago)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Set axis limits
    all_x = [pos[0] for pos in positions.values()]
    all_y = [pos[1] for pos in positions.values()]
    ax.set_xlim(min(all_x) - 0.5, max(all_x) + 0.5)
    ax.set_ylim(min(all_y) - 0.5, max(all_y) + 1)
    
    # DON'T invert y-axis - keep it normal so time goes from present (top) to past (bottom)
    # ax.invert_yaxis()  # Commented out to flip the tree
    
    # Remove top and right spines for cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if ax is None:
        plt.tight_layout()
        
        if save_file:
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Tree plot saved to {save_file}")
            plt.close()
        else:
            return fig, ax
    else:
        return fig, ax

def plot_sfs(sfs, title="Site Frequency Spectrum", figsize=(10, 6), save_file=None):
    """
    Plot the Site Frequency Spectrum as a bar plot.
    
    Args:
        sfs: Array of SFS values
        title: Title for the plot
        figsize: Figure size tuple
        save_file: Optional filename to save the plot
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create bar plot
    x_positions = np.arange(len(sfs))
    ax.bar(x_positions, sfs)
    
    # Customize
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Proportion')
    ax.set_title(title)
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{i+1}' for i in range(len(sfs))])
    
    plt.tight_layout()
    
    if save_file:
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        print(f"SFS plot saved to {save_file}")
        plt.close()
    else:
        plt.show()

def plot_tree_and_sfs(edges_df, sfs, title="Coalescent Tree and SFS", figsize=(12, 10), save_file=None):
    """
    Plot the tree on top and SFS below in a single figure.
    
    Args:
        edges_df: DataFrame with edge list
        sfs: Array of SFS values
        title: Title for the combined plot
        figsize: Figure size tuple
        save_file: Optional filename to save the plot
    """
    # Create subplots: tree on top, SFS below
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, height_ratios=[2, 1])
    
    # Plot tree on top subplot
    plot_tree_tskit_style(edges_df, title="", figsize=figsize, save_file=None, ax=ax1)
    
    # Plot SFS on bottom subplot
    x_positions = np.arange(len(sfs))
    ax2.bar(x_positions, sfs)
    ax2.set_xlabel('Frequency')
    ax2.set_ylabel('Proportion')
    ax2.set_title('Site Frequency Spectrum')
    ax2.set_xticks(x_positions)
    ax2.set_xticklabels([f'{i+1}' for i in range(len(sfs))])
    
    # Set overall title
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    
    if save_file:
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        print(f"Combined plot saved to {save_file}")
        plt.close()
    else:
        plt.show()

def calc_T_matrix_vectorized(edges_df, num_samples):
    """
    Fully vectorized version using numpy operations.
    Fastest for medium-sized trees.
    """
    # Build node mapping
    all_nodes = set()
    for _, row in edges_df.iterrows():
        all_nodes.add(row['parent'])
        all_nodes.add(row['child'])
    
    node_to_idx = {node: idx for idx, node in enumerate(sorted(all_nodes))}
    n_nodes = len(all_nodes)
    
    # Initialize distance matrix
    dist_matrix = np.full((n_nodes, n_nodes), np.inf)
    np.fill_diagonal(dist_matrix, 0)
    
    # Fill adjacency matrix
    for _, row in edges_df.iterrows():
        parent_idx = node_to_idx[row['parent']]
        child_idx = node_to_idx[row['child']]
        branch_length = row['branch_length']
        
        dist_matrix[parent_idx, child_idx] = branch_length
        dist_matrix[child_idx, parent_idx] = branch_length
    
    # Floyd-Warshall (vectorized)
    for k in range(n_nodes):
        dist_matrix = np.minimum(dist_matrix, dist_matrix[:, k:k+1] + dist_matrix[k:k+1, :])
    
    # Find leaf nodes
    parents = set(row['parent'] for _, row in edges_df.iterrows())
    leaf_nodes = [node for node in all_nodes if node not in parents]
    leaf_indices = np.array([node_to_idx[node] for node in leaf_nodes])
    
    # Vectorized pairwise distance calculation
    leaf_distances = dist_matrix[np.ix_(leaf_indices, leaf_indices)]
    
    # Get upper triangle (excluding diagonal)
    upper_triangle = np.triu(leaf_distances, k=1)
    
    # Calculate mean
    non_zero_distances = upper_triangle[upper_triangle != 0]
    return np.mean(non_zero_distances) if len(non_zero_distances) > 0 else 0

def calc_pi_from_SFS(sfs, num_samples, u):

    """
    Calculate π (nucleotide diversity) from the Site Frequency Spectrum.
    
    π = Σ(i=1 to n-1) SFS[i] * i * (n-i) / (n choose 2)
    """
    pi = 0
    for i in range(len(sfs)):
        # SFS[i] represents frequency i+1 (since SFS is 0-indexed)
        frequency = i + 1
        pi += sfs[i] * frequency * (num_samples - frequency)
    
    # Normalize by total number of pairs
    total_pairs = num_samples * (num_samples - 1) / 2
    pi = 2 * u * pi / total_pairs
    
    return pi

def calc_SFS_from_edges(edges_df, num_samples):
    unfolded_length = num_samples - 1
    folded_length = int(np.floor(num_samples / 2))

    unfolded_sfs = np.zeros(unfolded_length)
    for bin in range(len(unfolded_sfs)):
        # Count edges where leaves_below == bin + 1
        unfolded_sfs[bin] = edges_df.loc[edges_df['leaves_below'] == bin + 1, 'branch_length'].sum()
    
    folded_sfs = np.zeros(folded_length)
    for bin in range(len(folded_sfs)):
        if bin == folded_length - 1 and num_samples % 2 == 0:
            # For the middle bin when n is even
            folded_sfs[bin] = unfolded_sfs[bin]
        else:
            folded_sfs[bin] = unfolded_sfs[bin] + unfolded_sfs[num_samples - bin - 2]

    return folded_sfs