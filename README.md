# BottleneckModel
Contains the code used to run simulations and create figures for ...

Figures 2-4 can be run directly in the folder. Figure 2 runs a new set of simulations according to the parameters given. Figures 3 and 4 are perform calculations.

Figures 5 and S1 can be generated from the simulation data in these folders. This analysis can be performed by running ./analysis/pair_analyzer.sh in folder figure_5 or figure_S1. Running new replicates/parameters should be done by setting parameters in launch_simulation.py and then running ./submit_sfs.sh on a longleaf cluster. Modification to the SBATCH scripts will be needed.