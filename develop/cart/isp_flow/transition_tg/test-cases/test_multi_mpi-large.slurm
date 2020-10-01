#!/bin/bash

# ==============================================================================
# SLURM script for transition TG example -- hybrid
# ------------------------------------------------------------------------------

# Take care that the number of tasks matches the number of partitions defined by 
# parameter `np` in `isp_flow__sdc_test.prm`, i.e., ntasks = np(1)*np(2)*np(3)

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=15                            # number of nodes
#SBATCH --ntasks=343                          # max number of MPI processes
##SBATCH --ntasks-per-socket=1                # max num MPI proc per socket
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-02:00:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=tg3d                       # sensible job name
#SBATCH --output=test_multi_mpi-large_%j.out

# ------------------------------------------------------------------------------
# Execution

export KMP_AFFINITY=compact,1,0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun ../isp_flow__sdc_test test_multi_mpi-large

