#!/usr/bin/env bash

# run with
#  sbatch isp_flow__sdc_test.slurm
#
# ------------------------------------------------------------------------------
# Control

#SBATCH --account=p_lvnsm                     # project (triton-ism|p_lvnsm)
#SBATCH --nodes=1                             # number of nodes
#SBATCH --ntasks=16                           # max number of MPI processes
#SBATCH --cpu-freq=highm1                     # CPU frequency in kHz

#SBATCH --partition=haswell                   # partition
#SBATCH --cpus-per-task=1                     # number of cores per process
#SBATCH --mem=60000                           # memory per node in MB

#SBATCH --time=00-01:55:00                    # dd-hh:mm:ss
#SBATCH --exclusive
#SBATCH --job-name=isp_flow__sdc_test
#SBATCH --output=isp_flow__sdc_test_%j.out

date

# ------------------------------------------------------------------------------
# Environment

module load compiler/intel

# ------------------------------------------------------------------------------
# Execution


N_PART=16
N_SUB=3
N_CPI=0
RANGE_N_SWEEP="6"
RANGE_DT="0.0025"

BASE="isp_flow__sdc_test"

for N_SWEEP in $RANGE_N_SWEEP; do
for DT in $RANGE_DT; do

    echo ""
    echo "========================================================================="
    echo "n_sweep = " $N_SWEEP ", dt =" $DT

    CASE=${BASE}_${N_SUB}-${N_SWEEP}-${I_PRC}-dt${DT}

    sed -e "s/<n_sub>/$N_SUB/g" \
        -e "s/<n_sweep>/$N_SWEEP/g" \
        -e "s/<n_cpi>/$N_CPI/g" \
        -e "s/<dt>/$DT/g" \
        ${BASE}.tmpl > ${CASE}.prm

    srun -n ${N_PART} isp_flow__sdc_test ${CASE}

done
done

date

