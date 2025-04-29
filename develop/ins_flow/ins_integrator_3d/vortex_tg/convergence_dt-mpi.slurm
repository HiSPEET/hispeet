#!/bin/bash

# ==============================================================================
# SLURM script for dt-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=16                           # max number of MPI processes
#SBATCH --cpus-per-task=1                     # number of cores per process
#SBATCH --mem=100000                          # memory per node in MB

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
#SBATCH --job-name=convergence_dt             # sensible job name
#SBATCH --output=convergence_dt_%j.out        # output file

# for performance measurements:
##SBATCH --cpu-freq=highm1                     # CPU frequency in kHz
##SBATCH --exclusive                           # don't share nodes with others

# ------------------------------------------------------------------------------
# Settings (overriding those in convergence_dt.sh)

export MPIRUN="srun"

# number of partitions in directions 1-2 such that NP*NP equals ntasks
export NP="4"

# number of elements per partition in directions 1-2
export EP="2"

# max number of time step subdivisions per series
ST_MAX="10"

# ------------------------------------------------------------------------------
# Execution

./convergence_dt.sh
