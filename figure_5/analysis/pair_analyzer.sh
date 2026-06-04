python analysis/simple_analyze_tree_sim.py --grid --yscale linear --output plots/sfs_grid.pdf --input-file "sfs_simulation/simulation_results.jsonl" \
  --panel 9900 100 0.0002 0.0097 \
  --panel 9500 500 0.0002 0.0330 \
  --panel 5000 5000 0.0002 0.0574\
  --panel 1000 9000 0.0002 0.0460

python analysis/simple_analyze_tree_sim.py --grid --yscale log --output plots/sfs_grid_logy.pdf --input-file "sfs_simulation/simulation_results.jsonl" \
  --panel 9900 100 0.0002 0.0097 \
  --panel 9500 500 0.0002 0.0330 \
  --panel 5000 5000 0.0002 0.0574\
  --panel 1000 9000 0.0002 0.0460
