#!/bin/bash

# ==============================================================================
# SLURM script for elliptic_test__ip_ci
# ------------------------------------------------------------------------------

# Take care that the number of tasks matches the number of partitions defined by 
# parameter `np` in `elliptic_test__ip_ci.prm`, i.e., ntasks = np(1)*np(2)*np(3)

#SBATCH --account=p_hispeet                      # project
#SBATCH --nodes=1                                # number of nodes
#SBATCH --ntasks=1                               # max number of MPI processes
#SBATCH --cpu-freq=highm1                        # CPU frequency in kHz

#SBATCH --partition=haswell                      # partition
#SBATCH --cpus-per-task=12                       # number of cores per process
##SBATCH --mem=60000                             # memory per node in MB

#SBATCH --time=00-00:10:00                       # dd-hh:mm:ss
#SBATCH --exclusive                              # don't share nodes with others
#SBATCH --job-name=elliptic_test__ip_ci          # sensible job name
#SBATCH --output=elliptic_test__ip_ci_%j.out     # output file

# ------------------------------------------------------------------------------

export PROGRAM="elliptic_test__ip_ci"

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
  echo ""
  echo "SLURM_JOB_ID        = $SLURM_JOB_ID"
  echo "SLURM_JOB_PARTITION = $SLURM_JOB_PARTITION"
  echo "SLURM_JOB_NODELIST  = $SLURM_JOB_NODELIST"
  echo "RROGRAM             = $PROGRAM"
  echo ""
fi

# ------------------------------------------------------------------------------
# Execution

export OMP_NUM_THREADS=1
srun ./$PROGRAM

export OMP_NUM_THREADS=2
srun ./$PROGRAM

export OMP_NUM_THREADS=4
srun ./$PROGRAM

export OMP_NUM_THREADS=8
srun ./$PROGRAM

export OMP_NUM_THREADS=12
srun ./$PROGRAM

