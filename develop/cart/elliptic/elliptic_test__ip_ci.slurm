#!/bin/bash

# ==============================================================================
# SLURM script for elliptic_test__ip_ci
# ------------------------------------------------------------------------------

# Take care that the number of tasks matches the number of partitions defined by 
# parameter `np` in `isp_flow__sdc_test.prm`, i.e., ntasks = np(1)*np(2)*np(3)

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=1                            # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
##SBATCH --cpus-per-task=1                    # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-00:10:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=elliptic_test__ip_ci       # sensible job name
#SBATCH --output=elliptic_test__ip_ci_%j.out  # output file

# ------------------------------------------------------------------------------
# Execution

srun ../elliptic_test__ip_ci

