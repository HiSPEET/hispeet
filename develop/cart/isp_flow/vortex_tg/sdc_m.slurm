#!/bin/bash

# ==============================================================================
# SLURM script for SDC(M,K=3*M)-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_lvnsm                     # project (p_lvnsm)
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=4                            # max number of MPI processes
#SBATCH --cpu-freq=highm1                     # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
##SBATCH --cpus-per-task=1                    # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=sdc_m                      # sensible job name
#SBATCH --output=sdc_m_%j.out                 # output file

# ------------------------------------------------------------------------------
# Environment

module load ${ENV_MODULE_COMPILER}
export EXEC="srun"

# ------------------------------------------------------------------------------
# Settings (overriding those in sdc_m.sh)

# number of partitions in directions 1-2 such that NP1*NP2 = ntasks
export NP1="2"
export NP2="2"

# number of elements per partition in directions 1-2
export EP="1"

# final time
export T_END="1" 

# ------------------------------------------------------------------------------
# Execution

./sdc_m.sh
