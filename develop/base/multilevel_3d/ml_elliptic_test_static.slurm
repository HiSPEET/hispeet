#!/bin/bash

# ==============================================================================
# SLURM script for ML elliptic static test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                        # project
#SBATCH --nodes=1                                  # number of nodes
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1                          # number of cores per process
#SBATCH --mem=100000                               # memory per node in MB

#SBATCH --time=00-00:55:00                         # dd-hh:mm:ss
#SBATCH --job-name=ml_elliptic_test_static         # sensible job name
#SBATCH --output=ml_elliptic_test_static_%j.out    # output file

# for performance measurements:
##SBATCH --cpu-freq=highm1                         # CPU frequency in kHz
##SBATCH --exclusive                               # don't share nodes

# ------------------------------------------------------------------------------

date

srun ml_elliptic_test_static

date

