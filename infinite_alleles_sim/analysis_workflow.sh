#! /bin/bash

# Combine the data
# ./analysis/data_combiner.sh

python analysis/pi.py

Rscript analysis/pi_plot.R

Rscript analysis/sfs_sampled.R