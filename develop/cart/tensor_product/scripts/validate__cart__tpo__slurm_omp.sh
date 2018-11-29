#!/bin/bash

# ==============================================================================
# SLURM script for testing tensor-product operators with OpenMP
# ------------------------------------------------------------------------------

#SBATCH --account=p_lvnsm                     # project (triton-ism|p_lvnsm)
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=1                            # max number of MPI processes
#SBATCH --ntasks-per-socket=1                 # only one process per socket
#SBATCH --cpu-freq=highm1                     # highest non-turbo frequency

#SBATCH --partition=haswell                   # partition
#SBATCH --cpus-per-task=12                    # number of cores per process
#SBATCH --mem=60000                           # memory per node in MB

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
#SBATCH --exclusive                           # don't share nodes with others
#SBATCH --job-name=validate-omp               # sensible job name
#SBATCH --output=validate-omp_%j.out          # output file

# ------------------------------------------------------------------------------

module load ${ENV_MODULE_COMPILER}
export EXEC="srun"

./validate__cart__tpo.sh
