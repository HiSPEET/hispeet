#!/bin/bash

# ==============================================================================
# SLURM script for SDC(M,K)-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=8                            # max number of MPI processes
#SBATCH --cpu-freq=highm1                     # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-00:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=sdc_k                      # sensible job name
#SBATCH --output=sdc_k_%j.out                 # output file

# ------------------------------------------------------------------------------

if [ X"$SLURM_STEP_ID" = "X" -a X"$SLURM_PROCID" = "X"0 ]
then
  echo ""
  echo "SLURM_JOB_ID        = $SLURM_JOB_ID"
  echo "SLURM_JOB_PARTITION = $SLURM_JOB_PARTITION"
  echo "SLURM_JOB_NODELIST  = $SLURM_JOB_NODELIST"
  echo ""
fi


export EXEC="srun"

# ------------------------------------------------------------------------------
# Settings (overriding those in sdc_k.sh)

# number of partitions in directions 1-3 such that NP1*NP2*NP3 = ntasks
export NP1="2"
export NP2="2"
export NP3="2"

# number of elements per partition in directions 1-3
export EP1="1"
export EP2="1"
export EP3="1"

# range of number of sweeps
export RANGE_N_SWEEP="0 1 2 3 4 5 6"

# ------------------------------------------------------------------------------
# Execution

./sdc_k.sh
