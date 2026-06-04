import numpy as np
import pandas as pd


def load_jsonl_data(output_file='fast_tree_simulation/simulation_results.jsonl'):
    """Load JSONL data"""
    import pandas as pd
    import json
    
    data = []
    skipped_lines = []
    
    with open(output_file, 'rb') as f:  # Read as binary to detect null bytes
        for line_num, line_bytes in enumerate(f, 1):
            # Check for null bytes
            if b'\x00' in line_bytes:
                null_count = line_bytes.count(b'\x00')
                # Try to strip null bytes and see if there's valid JSON after
                line_cleaned = line_bytes.replace(b'\x00', b'').decode('utf-8', errors='ignore').strip()
                if not line_cleaned:
                    skipped_lines.append((line_num, f"Line contains only {null_count} null bytes"))
                    continue
                else:
                    # Try to parse the cleaned line
                    try:
                        parsed = json.loads(line_cleaned)
                        data.append(parsed)
                        print(f"Warning: Line {line_num} had {null_count} null bytes, but cleaned version parsed successfully")
                        continue
                    except json.JSONDecodeError:
                        skipped_lines.append((line_num, f"Line contains {null_count} null bytes and cannot be parsed even after cleaning"))
                        continue
            
            # Decode and process normal lines
            try:
                line = line_bytes.decode('utf-8')
            except UnicodeDecodeError:
                skipped_lines.append((line_num, "Line contains invalid UTF-8 encoding"))
                continue
            
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            
            try:
                parsed = json.loads(line)
                data.append(parsed)
            except json.JSONDecodeError as e:
                skipped_lines.append((line_num, f"JSON decode error: {e}"))
                continue
    
    if skipped_lines:
        print(f"\nSkipped {len(skipped_lines)} problematic lines:")
        for line_num, reason in skipped_lines[:10]:  # Show first 10
            print(f"  Line {line_num}: {reason}")
        if len(skipped_lines) > 10:
            print(f"  ... and {len(skipped_lines) - 10} more")
    
    if not data:
        raise ValueError(f"No valid JSON data found in {output_file}")
    
    print(f"Successfully loaded {len(data)} valid JSON records")
    df = pd.DataFrame(data)
    df['sfs'] = df['sfs'].apply(lambda x: np.array(x))
    return df