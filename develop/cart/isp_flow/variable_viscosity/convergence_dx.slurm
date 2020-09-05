#!/bin/bash

# ==============================================================================
# SLURM script for dt-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=triton                      # project (triton-ism|p_lvnsm)
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=1                            # max number of MPI processes
#SBATCH --cpu-freq=highm1                     # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
#SBATCH --cpus-per-task=12                    # number of cores per process
#SBATCH --mem=60000                           # memory per node in MB

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
#SBATCH --exclusive                           # don't share nodes with others
#SBATCH --job-name=convergence_dt             # sensible job name
#SBATCH --output=convergence_dt_%j.out        # output file

# ------------------------------------------------------------------------------

export EXEC="srun"

./convergence_dx.sh
