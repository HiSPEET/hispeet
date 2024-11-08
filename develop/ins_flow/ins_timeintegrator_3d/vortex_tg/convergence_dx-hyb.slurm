#!/bin/bash

# ==============================================================================
# SLURM script for dt-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=4                             # number of nodes
#SBATCH --ntasks=16                           # max number of MPI processes
#SBATCH --cpus-per-task=12                    # number of cores per process
#SBATCH --mem=100000                          # memory per node in MB

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
#SBATCH --job-name=convergence_dx             # sensible job name
#SBATCH --output=convergence_dx_%j.out        # output file

# for performance measurements:
##SBATCH --cpu-freq=highm1                     # CPU frequency in kHz
##SBATCH --exclusive                           # don't share nodes with others

# ------------------------------------------------------------------------------

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
  echo ""
  echo "SLURM_JOB_ID        = $SLURM_JOB_ID"
  echo "SLURM_JOB_PARTITION = $SLURM_JOB_PARTITION"
  echo "SLURM_JOB_NODELIST  = $SLURM_JOB_NODELIST"
  echo ""
fi

# ------------------------------------------------------------------------------

export MPIRUN="srun"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# time step width
export DT="1e-5"

# test configurations
export NP=( 2 4 4 )  # number of partitions in directions 1-2
export EP=( 1 1 2 )  # number of elements per partition in directions 1-2 

# ------------------------------------------------------------------------------
# Execution

./convergence_dx.sh
