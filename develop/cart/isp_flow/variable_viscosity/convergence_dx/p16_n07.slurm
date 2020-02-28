#!/bin/bash

# ==============================================================================
# SLURM script for h-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=3                             # number of nodes
#SBATCH --ntasks=49                           # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-32:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=p16_n07                    # sensible job name
#SBATCH --output=p16_n07_%j.out               # output file

# ------------------------------------------------------------------------------
# Environment

module load @ENV_MODULE_COMPILER@

PO_U=16
PO_P=15
PO_Q=24

NP1=1
NP2=7
NP3=7
EP1=7
EP2=1
EP3=1

CASE="p16_n07"
LOG_FILE=.log
PRM_FILE=.prm
DAT_FILE=.dat

sed -e "s/<np1>/$NP1/g" \
    -e "s/<np2>/$NP2/g" \
    -e "s/<np3>/$NP3/g" \
    -e "s/<ep1>/$EP1/g" \
    -e "s/<ep2>/$EP2/g" \
    -e "s/<ep3>/$EP3/g" \
    -e "s/<po_u>/$PO_U/g" \
    -e "s/<po_p>/$PO_P/g" \
    -e "s/<po_q>/$PO_Q/g" \
    isp_flow__sdc_test.tmpl > 

date > 

srun ../../../isp_flow__sdc_test  2>&1 | tee -a 

date >> 

sleep 10

grep -e "#      t" -m 1  >  
grep -e "#last#$"        >> 
