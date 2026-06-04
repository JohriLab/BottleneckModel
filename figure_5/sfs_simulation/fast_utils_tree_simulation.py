import random
import numpy as np
import pandas as pd
from collections import defaultdict
import networkx as nx
from itertools import combinations

def simulate_tree(k, m, N, d, num_samples):
    '''simulate a coalescent tree with migration, deme drift, and coalescence'''
    '''
    Args:
        k: bottleneck size
        m: migration rate
        N: population size
        d: number of demes
        num_samples: number of samples

    Returns:
        edges_df: DataFrame with edge list
    '''
    
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

        ### Coalescence (1/k): ###
        deme_groups = defaultdict(list)
        for sample in samples:
            deme_groups[sample['deme']].append(sample)
            
        new_samples = []
        for deme, deme_samples in deme_groups.items():
            if len(deme_samples) > 1:
                coalescing_pairs = []
                
                for pair in combinations(deme_samples, 2):
                    if np.random.random() < 1/k:
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

        ### Migration: ###
        migration_mask = np.random.random(len(samples)) < m
        for i, sample in enumerate(samples):
            if migration_mask[i]:
                current_deme = sample['deme']
                # new_deme = random.randint(1, d-1)
                # if new_deme >= current_deme:
                #     new_deme += 1
                new_deme = random.randint(1, d)
                sample['deme'] = new_deme

        ### Deme Drift: ###
        unique_demes = list(set(sample['deme'] for sample in samples))
        deme_parent_mapping = {}
        for deme in unique_demes:
            parent_deme = random.randint(1, d)
            deme_parent_mapping[deme] = parent_deme
        
        for sample in samples:
            sample['deme'] = deme_parent_mapping[sample['deme']]

        ### Coalescence (1/N): ###
        deme_groups = defaultdict(list)
        for sample in samples:
            deme_groups[sample['deme']].append(sample)
            
        new_samples = []
        for deme, deme_samples in deme_groups.items():
            if len(deme_samples) > 1:
                coalescing_pairs = []
                
                for pair in combinations(deme_samples, 2):
                    if np.random.random() < 1/N:
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

def calc_SFS_from_edges(edges_df, num_samples, folded=True):
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

    if folded:
        return folded_sfs
    else:
        return unfolded_sfs
    return unfolded_sfs

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
    pi = u * pi / total_pairs
    
    return pi

def save_data_jsonl(k, m, d, N, u, num_samples, pi, sfs, output_file='simulation_results.jsonl'):
    """Save data in JSON Lines format - human readable
    
    Uses file locking to prevent corruption when multiple processes write concurrently.
    """
    import json
    import fcntl  # For file locking on Unix systems
    
    data = {
        'k': k, 'm': m, 'd': d, 'N': N, 'u': u,
        'num_samples': num_samples, 'pi': pi,
        'sfs': sfs.tolist()
    }
    
    # Use file locking to prevent concurrent write corruption
    with open(output_file, 'a') as f:
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)  # Exclusive lock
            f.write(json.dumps(data) + '\n')
            f.flush()  # Ensure data is written before releasing lock
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)  # Release lock

