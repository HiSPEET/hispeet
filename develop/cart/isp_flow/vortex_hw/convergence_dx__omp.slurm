#!/bin/bash

# ==============================================================================
# SLURM script for dx-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=1                            # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
##SBATCH --cpus-per-task=12                   # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-02:00:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=convergence_dt             # sensible job name
#SBATCH --output=convergence_dt_%j.out        # output file

# ------------------------------------------------------------------------------
# Settings (overriding those in convergence_dx.sh)

export EXEC="srun"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# number of partitions in directions 1-2 such that NP1*NP2 = ntasks
export NP1="2"
export NP2="2"

# investigated range of elements per partition in directions 1-2
export RANGE_EP="1 2 4"

# ------------------------------------------------------------------------------
# Execution

./convergence_dx.sh
