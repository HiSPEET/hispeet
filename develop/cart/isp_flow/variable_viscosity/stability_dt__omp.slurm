#!/bin/bash

# ==============================================================================
# SLURM script for dt-stability test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_lvnsm                     # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=1                            # max number of MPI processes
#SBATCH --cpu-freq=highm1                     # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
#SBATCH --cpus-per-task=12                    # number of cores per process
#SBATCH --mem=60000                           # memory per node in MB

#SBATCH --time=32:00:00                       # dd-hh:mm:ss
##SBATCH --exclusive                           # don't share nodes with others
#SBATCH --job-name=stability_dt               # sensible job name
#SBATCH --output=stability_dt_%j.out          # output file

# ------------------------------------------------------------------------------

# run job with unbuffered output
export EXEC="srun -u"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ------------------------------------------------------------------------------
# Settings (overriding those in sdc_k.sh)

# number of partitions in directions 1-3 such that NP1*NP2*NP3 = ntasks
export NP1="2"
export NP2="2"
export NP3="2"

# polynomial orders for u and p and for integrating the nonlinear terms
export PO_U="16"
export PO_P="15"
export PO_Q="24"

# SDC parameters
export N_SUB="1"
export N_SWEEP="0"

# time stepping
export DT0="0"
export DT1="1"
export NBI="15"

# ------------------------------------------------------------------------------
# Execution

./stability_dt.sh
