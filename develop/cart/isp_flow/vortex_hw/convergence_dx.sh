#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

# number of partitions in directions 1-2
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP=$((${NP1} * ${NP2}))

# use preset execution command, or mpirun, if not set
MPIRUN=${MPIRUN:-"mpirun"}
EXEC=${EXEC:-"${MPIRUN} -n ${NP}"}

# path to program
PROGRAM="../isp_flow__sdc_test"

# test case
CASE="convergence_dx"

# investigated range of elements per partition in directions 1-2
RANGE_EP=${RANGE_EP:-"2 4 8"}

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U="5"
PO_P="4"
PO_Q="9"

# time step size
DT="5e-5"

# SDC parameters
N_SUB=${N_SUB:-"1"}
N_SWEEP=${N_SWEEP:-"1"}

date > ${CASE}.log

for EP in $RANGE_EP; do

    echo "========================================================================="
    echo "ep =" $EP
    echo
    sed -e "s/<np1>/$NP1/g" \
        -e "s/<np2>/$NP2/g" \
        -e "s/<ep1>/$EP/g" \
        -e "s/<ep2>/$EP/g" \
        -e "s/<po_u>/$PO_U/g" \
        -e "s/<po_p>/$PO_P/g" \
        -e "s/<po_q>/$PO_Q/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<n_sub>/$N_SUB/g" \
        -e "s/<n_sweep>/$N_SWEEP/g" \
        isp_flow__sdc_test.tmpl > \
        isp_flow__sdc_test.prm

    ${EXEC} ${PROGRAM} 2>&1 | tee -a ${CASE}.log

done

grep -e "#      t" -m 1 ${CASE}.log >  ${CASE}.dat
grep -e "#last#$"       ${CASE}.log >> ${CASE}.dat

date >> ${CASE}.log
