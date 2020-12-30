#!/bin/bash

# ==============================================================================
# SLURM script for sequential testing of tensor-product operators
# ------------------------------------------------------------------------------

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=1                            # max number of MPI processes
#SBATCH --cpu-freq=highm1                     # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
#SBATCH --cpus-per-task=12                    # number of cores per process
#SBATCH --mem=60000                           # memory per node in MB

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
#SBATCH --exclusive                           # don't share nodes with others
#SBATCH --job-name=validate-seq               # sensible job name
#SBATCH --output=validate-seq_%j.out          # output file

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

export EXEC="srun"

./validate__tpo.sh
