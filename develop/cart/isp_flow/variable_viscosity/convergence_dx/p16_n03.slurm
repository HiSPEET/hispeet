#!/bin/bash

# ==============================================================================
# SLURM script for h-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=2                             # number of nodes
#SBATCH --ntasks=27                           # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-32:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=p16_n03                    # sensible job name
#SBATCH --output=p16_n03_%j.out               # output file

# ------------------------------------------------------------------------------
# Environment

module load @ENV_MODULE_COMPILER@

PO_U=16
PO_P=15
PO_Q=24

NP1=3
NP2=3
NP3=3
EP1=1
EP2=1
EP3=1

CASE="p16_n03"

source run_case.sh
