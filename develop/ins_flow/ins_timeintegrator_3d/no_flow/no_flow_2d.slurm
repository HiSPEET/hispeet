#!/bin/bash

# ==============================================================================
# SLURM script for 2D no-flow test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=2                             # number of nodes
#SBATCH --ntasks=64                           # max number of MPI processes
#SBATCH --cpus-per-task=1                     # number of cores per process
#SBATCH --mem=100000                          # memory per node in MB

#SBATCH --time=00-07:55:00                    # dd-hh:mm:ss
#SBATCH --job-name=no_flow_2d                 # sensible job name
#SBATCH --output=no_flow_2d%j.out             # output file

# for performance measurements:
##SBATCH --cpu-freq=highm1                     # CPU frequency in kHz
##SBATCH --exclusive                           # don't share nodes with others

# ------------------------------------------------------------------------------
# Execution

srun ../ins_timeintegrator_3d_test ins_timeintegrator_3d_test_2d
