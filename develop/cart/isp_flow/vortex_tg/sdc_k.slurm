#!/bin/bash

# ==============================================================================
# SLURM script for SDC(M,K)-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=4                            # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
##SBATCH --cpus-per-task=1                    # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-04.00:00                    # dd-hh:mm:ss
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

# number of partitions in directions 1-2 such that NP1*NP2 = ntasks
export NP1="2"
export NP2="2"

# number of elements per partition in directions 1-2
export EP="1"

# final time
export T_END="1" 

# number of subintervals for SDC
export N_SUB="3"

# max time step size, may be increased to 0.0625 for periodic case
export DT_MAX="0.015625"

# max number of time step reductions per series 
export ST_MAX="8"

# ------------------------------------------------------------------------------
# Execution

./sdc_k.sh
