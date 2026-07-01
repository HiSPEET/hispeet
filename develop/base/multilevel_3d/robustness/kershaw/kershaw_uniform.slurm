#!/bin/bash

# ==============================================================================
# SLURM script for ML elliptic kershaw test
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                        ### adjust ###
#SBATCH --nodes=1                                  ### adjust ###
#SBATCH --ntasks=12                                ### adjust ###

#SBATCH --cpus-per-task=1
#SBATCH --mem=500000
#SBATCH --time=00-00:55:00
#SBATCH --job-name=kershaw_uniform
#SBATCH --output=kershaw_uniform_%j.out

# for performance measurements:
##SBATCH --cpu-freq=highm1
##SBATCH --exclusive

# ==============================================================================
# General

PROGRAM="../../../ml_elliptic_test_static"
CASE="kershaw_uniform"

# ==============================================================================
# Execution

date

srun -n $SLURM_NTASKS $PROGRAM $CASE 2>&1 | tee ${CASE}-np_${SLURM_NTASKS}.log

# ------------------------------------------------------------------------------

date

