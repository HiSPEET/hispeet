#!/bin/bash

# ==============================================================================
# SLURM script for h-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=11                           # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-32:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=p11_n11                    # sensible job name
#SBATCH --output=p11_n11_%j.out               # output file

# ------------------------------------------------------------------------------
# Environment

module load @ENV_MODULE_COMPILER@

PO_U=11
PO_P=10
PO_Q=17

NP1=1
NP2=1
NP3=11
EP1=11
EP2=11
EP3=1

CASE="p11_n11"

source run_case.sh
