#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python "${SCRIPT_DIR}/simple_analyze_tree_sim.py" \
  --grid --output "${SCRIPT_DIR}/plots/sfs_grid_7panel.pdf" --input-file "${SCRIPT_DIR}/../sfs_simulation/simulation_results.jsonl" \
  --panel 5000 5000 0 0.1841 \
  --panel 5000 5000 0.00002 0.1575 \
  --panel 5000 5000 0.0001 0.0936 \
  --panel 5000 5000 0.0002 0.0565 \
  --panel 5000 5000 0.0003 0.0376 \
  --panel 5000 5000 0.0005 0.0202 \
  --panel 5000 5000 0.001 0.0081 \
  --panel 5000 5000 0.002 0.0035