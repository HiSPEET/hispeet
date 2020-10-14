#!/bin/bash

# ==============================================================================
# SLURM script for h-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=2                             # number of nodes
#SBATCH --ntasks=48                           # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-32:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=p11_n12                    # sensible job name
#SBATCH --output=p11_n12_%j.out               # output file

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
# Environment

PO_U=11
PO_P=10
PO_Q=17

NP1=4
NP2=4
NP3=3
EP1=3
EP2=3
EP3=4

CASE="p11_n12"

source run_case.sh
