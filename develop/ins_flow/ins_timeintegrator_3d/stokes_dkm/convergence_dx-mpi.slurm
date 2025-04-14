#!/bin/bash

# ==============================================================================
# SLURM script for dx-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=64                           # max number of MPI processes
#SBATCH --cpus-per-task=1                     # number of cores per process
#SBATCH --mem=100000                          # memory per node in MB

#SBATCH --time=00-23:55:00                    # dd-hh:mm:ss
#SBATCH --job-name=convergence_dx             # sensible job name
#SBATCH --output=convergence_dx_%j.out        # output file

# for performance measurements:
##SBATCH --cpu-freq=highm1                     # CPU frequency in kHz
##SBATCH --exclusive                           # don't share nodes with others

# ------------------------------------------------------------------------------
# Settings (overriding those in convergence_dt.sh)

export EXEC="srun"
export NPROC_MAX="$SLURM_NTASKS"
export PO_U="4" 
export PO_P="$((${PO_U} - 1))"
export PE_MU="1.0"

# ------------------------------------------------------------------------------
# Execution

./convergence_dx.sh
