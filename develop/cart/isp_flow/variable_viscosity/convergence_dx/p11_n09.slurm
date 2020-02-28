#!/bin/bash

# ==============================================================================
# SLURM script for h-convergence test using MPI
# ------------------------------------------------------------------------------

# Adjust the following arguments as appropriate

#SBATCH --account=p_hispeet                   # project
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=9                            # max number of MPI processes
##SBATCH --cpu-freq=highm1                    # CPU frequency in kHz

##SBATCH --partition=haswell                  # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
##SBATCH --mem=60000                          # memory per node in MB

#SBATCH --time=00-32:55:00                    # dd-hh:mm:ss
##SBATCH --exclusive                          # don't share nodes with others
#SBATCH --job-name=p11_n09                    # sensible job name
#SBATCH --output=p11_n09_%j.out               # output file

# ------------------------------------------------------------------------------
# Environment

module load @ENV_MODULE_COMPILER@

PO_U=11
PO_P=10
PO_Q=17

NP1=3
NP2=3
NP3=1
EP1=3
EP2=3
EP3=9

CASE="p11_n09"
LOG_FILE=${CASE}.log
PRM_FILE=${CASE}.prm
DAT_FILE=${CASE}.dat

sed -e "s/<np1>/$NP1/g" \
    -e "s/<np2>/$NP2/g" \
    -e "s/<np3>/$NP3/g" \
    -e "s/<ep1>/$EP1/g" \
    -e "s/<ep2>/$EP2/g" \
    -e "s/<ep3>/$EP3/g" \
    -e "s/<po_u>/$PO_U/g" \
    -e "s/<po_p>/$PO_P/g" \
    -e "s/<po_q>/$PO_Q/g" \
    isp_flow__sdc_test.tmpl > ${PRM_FILE}

date > ${LOG_FILE}

srun ../../../isp_flow__sdc_test ${CASE} 2>&1 | tee -a 

date >> ${LOG_FILE}

sleep 10

grep -e "#      t" -m 1  >  ${LOG_FILE}
grep -e "#last#$"        >> ${LOG_FILE}
