#!/bin/bash

# ==============================================================================
# SLURM script for DG difusion test -- hybrid
# ------------------------------------------------------------------------------

# Take care that the number of tasks matches the number of partitions defined in 
# `dg_diffusion.prm`

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=2                            # max number of MPI processes
#SBATCH --ntasks-per-socket=1                 # max num MPI proc per socket
#SBATCH --cpus-per-task=12                    # number of cores per process
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
#SBATCH --mem=60000                           # memory per node in MB

#SBATCH --time=00-00:30:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=dg_diffusion3d_test        # sensible job name
#SBATCH --output=dg_diffusion3d_test_%j.out

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
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

OMP_NUM_THREADS=1  srun ./dg_diffusion3d_test
OMP_NUM_THREADS=2  srun ./dg_diffusion3d_test
OMP_NUM_THREADS=4  srun ./dg_diffusion3d_test
OMP_NUM_THREADS=8  srun ./dg_diffusion3d_test
OMP_NUM_THREADS=12 srun ./dg_diffusion3d_test

