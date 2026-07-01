#!/bin/bash

# ==============================================================================
# SLURM script for ML elliptic blade test using global refinement
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                        ### adjust ###
#SBATCH --nodes=1                                  ### adjust ###
#SBATCH --ntasks=64                                ### adjust ###

#SBATCH --cpus-per-task=1
#SBATCH --mem=500000
#SBATCH --time=00-07:55:00
#SBATCH --job-name=blade_global-ph
#SBATCH --output=blade_global-ph_%j.out

# for performance measurements:
##SBATCH --cpu-freq=highm1
##SBATCH --exclusive

# ==============================================================================
# General

export EXEC="srun"

# also consider
# export EXEC="srun --distribution=block:block"
# export EXEC="srun --distribution=block:cyclic" # default

# max number of root partition
export NPR=$SLURM_NTASKS
# note: a lower value of NPR reduces communication on coarser levels, 
# but increases the workload per partition

# ==============================================================================
# Execution

date

for np in 32 64 # 128 256 512 1024 2048 ...
do
    NPROC=$np ./blade_global-ph.sh
done

# ------------------------------------------------------------------------------

date

