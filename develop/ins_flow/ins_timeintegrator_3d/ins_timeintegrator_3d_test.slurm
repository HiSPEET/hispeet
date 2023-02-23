#!/bin/bash

# ==============================================================================
# SLURM script for INS time integrator 3D test
# ------------------------------------------------------------------------------

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=4                            # max number of MPI processes
##SBATCH --ntasks-per-socket=1                # max num MPI proc per socket
##SBATCH --cpus-per-task=1                    # number of cores per process
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
#SBATCH --mem=0                               # memory per node in MB (0 = all)

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=ins_timeintegrator_3d_test # sensible job name
#SBATCH --output=ins_timeintegrator_3d_test_%j.out

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
# Execution

export KMP_AFFINITY=compact,1,0

# 4 processes using 1 cpu
OMP_NUM_THREADS=1  srun -n 4 ./ins_timeintegrator_3d_test

