#!/bin/bash

# ==============================================================================
# SLURM script for h-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=8                            # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-32:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=p06_n02                    # sensible job name
#SBATCH --output=p06_n02_%j.out               # output file

# ------------------------------------------------------------------------------
# Environment

module load @ENV_MODULE_COMPILER@

PO_U=6
PO_P=5
PO_Q=10

NP1=2
NP2=2
NP3=2
EP1=1
EP2=1
EP3=1

CASE="p06_n02"

source run_case.sh
