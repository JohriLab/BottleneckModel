import numpy as np
import pandas as pd
from fast_utils_tree_simulation import *


def load_jsonl_data(output_file='simulation_results.jsonl'):
    """Load JSONL data"""
    import pandas as pd
    import json
    
    data = []
    with open(output_file, 'r') as f:
        for line in f:
            data.append(json.loads(line))
    
    df = pd.DataFrame(data)
    df['sfs'] = df['sfs'].apply(lambda x: np.array(x))
    return df