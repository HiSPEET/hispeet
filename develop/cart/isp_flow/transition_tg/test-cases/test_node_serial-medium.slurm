#!/bin/bash

# ==============================================================================
# SLURM script for transition TG example -- serial
# ------------------------------------------------------------------------------

# Take care that the number of tasks matches the number of partitions defined by 
# parameter `np` in `isp_flow__sdc_test.prm`, i.e., ntasks = np(1)*np(2)*np(3)

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=1                            # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-02:00:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=tg3d                       # sensible job name
#SBATCH --output=test_node_serial-medium_%j.out

# ------------------------------------------------------------------------------
# Execution

export OMP_NUM_THREADS=1
srun ../isp_flow__sdc_test test_node_serial-medium

