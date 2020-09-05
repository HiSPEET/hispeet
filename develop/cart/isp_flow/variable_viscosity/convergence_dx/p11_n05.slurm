#!/bin/bash

# ==============================================================================
# SLURM script for h-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=2                             # number of nodes
#SBATCH --ntasks=25                           # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-32:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=p11_n05                    # sensible job name
#SBATCH --output=p11_n05_%j.out               # output file

# ------------------------------------------------------------------------------
# Environment

PO_U=11
PO_P=10
PO_Q=17

NP1=1
NP2=5
NP3=5
EP1=5
EP2=1
EP3=1

CASE="p11_n05"

source run_case.sh
