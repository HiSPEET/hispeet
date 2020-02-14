#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l stability_inf-sup.sh`

# number of partitions in directions 1-2
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP=$((${NP1} * ${NP2}))

# use preset execution command, or mpirun, if not set
EXEC=${EXEC:-"mpirun -n ${NP}"}

# path to program
PROGRAM="../isp_flow__sdc_test"

# test case
CASE="stability_inf-sup"

# investigated range of elements per partition in directions 1-2
RANGE_EP=${RANGE_EP:-"2 4 8"}

# polynomial orders for u and p and size of time step
PO_U="3"
PO_P="2"
DT="1e-5"

# SDC parameters
N_SUB="1"
N_SWEEP="0 1"

date > ${CASE}.log

for EP in $RANGE_EP; do

    echo "========================================================================="
    echo " ep =" $EP
    echo
    sed -e "s/<np1>/$NP1/g" \
        -e "s/<np2>/$NP2/g" \
        -e "s/<ep1>/$EP/g" \
        -e "s/<ep2>/$EP/g" \
        -e "s/<po_u>/$PO_U/g" \
        -e "s/<po_p>/$PO_P/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<n_sub>/$N_SUB/g" \
        -e "s/<n_sweep>/$N_SWEEP/g" \
        isp_flow__sdc_test.tmpl > \
        isp_flow__sdc_test.prm

    ${EXEC} ${PROGRAM} 2>&1 | tee -a ${CASE}.log

done

grep -e "#     t" -m 1  ${CASE}.log >  ${CASE}.dat
grep -e ^" 1.00000E-01" ${CASE}.log >> ${CASE}.dat

date >> ${CASE}.log
