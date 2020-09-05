#!/bin/bash

#!!!!! MUSS NOCH FUER OMP ANGEPASST WERDEN !!!!!

# ==============================================================================
# SLURM script for SDC(M,K=3*M)-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=4                            # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=12                   # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-04:00:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=sdc_m                      # sensible job name
#SBATCH --output=sdc_m_%j.out                 # output file

# ------------------------------------------------------------------------------
# Settings (overriding those in sdc_m.sh)

export EXEC="srun"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# number of partitions in directions 1-2 such that NP1*NP2 = ntasks
export NP1="2"
export NP2="2"

# number of elements per partition in directions 1-2
export EP="1"

# final time
export T_END="1" 

# max number of subintervals for SDC
export N_SUB_MAX="3"

# max time step size, may be increased to 0.0625 for periodic case
export DT_MAX="0.015625"

# max number of time step reductions per series 
export ST_MAX="8"

# ------------------------------------------------------------------------------
# Execution

./sdc_m.sh
