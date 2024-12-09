#!/bin/bash

# ==============================================================================
# SLURM script for 3D no-flow test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=64                           # max number of MPI processes
#SBATCH --cpus-per-task=1                     # number of cores per process
#SBATCH --mem=100000                          # memory per node in MB

#SBATCH --time=00-23:55:00                    # dd-hh:mm:ss
#SBATCH --job-name=no_flow_3d                 # sensible job name
#SBATCH --output=no_flow_3d.out               # output file

# for performance measurements:
##SBATCH --cpu-freq=highm1                     # CPU frequency in kHz
##SBATCH --exclusive                           # don't share nodes with others

# ------------------------------------------------------------------------------
# Execution

date

srun ../ins_timeintegrator_3d_test ins_timeintegrator_3d_test_3d

date

