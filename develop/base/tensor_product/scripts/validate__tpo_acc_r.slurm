#!/bin/bash

# ==============================================================================
# SLURM script for running sequential test with validate__cart__tpo_diffusion
# ------------------------------------------------------------------------------

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=1                            # max number of MPI processes

#SBATCH --partition=gpu2                      # partition
#SBATCH --gres=gpu:1                          # use one GPU per node

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
#SBATCH --exclusive                           # don't share nodes with others
#SBATCH --job-name=validate-acc               # sensible job name
#SBATCH --output=validate-acc_%j.out          # stdout file, %j is the job id

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

./validate__tpo_r.sh
