#!/bin/bash

# ==============================================================================
# SLURM script for h-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=2                             # number of nodes
#SBATCH --ntasks=45                           # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-32:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=p06_n15                    # sensible job name
#SBATCH --output=p06_n15_%j.out               # output file

# ------------------------------------------------------------------------------
# Environment

PO_U=6
PO_P=5
PO_Q=10

NP1=5
NP2=3
NP3=3
EP1=3
EP2=5
EP3=5

CASE="p06_n15"

source run_case.sh
