#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

# number of partitions in directions 1-2
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP3=${NP3:-"1"}
NP=$((${NP1} * ${NP2} * ${NP3}))

# use preset execution command, or mpirun, if not set
MPIRUN=${MPIRUN:-"mpirun"}
EXEC=${EXEC:-"${MPIRUN} -n ${NP}"}

# path to program
PROGRAM="../isp_flow__sdc_test"

# test case
CASE="convergence_dt_sdc_2_4"

# elements per partition in directions 1-2
EP=${EP:-"4"}
EP3=${EP3:-"1"}
# polynomial orders for u and p and for integrating the nonlinear terms
PO_U="10"
PO_P="9"
PO_Q="16"
# final time
T_END=${T_END:-"0.25"}

# investigated range of time step sizes
RANGE_DT=\
" .05 .025 .0125 6.25d-3 3.125d-3 1.5625d-3 7.8125d-4 3.90625d-4" # 1.953125d-4"

# SDC parameters
N_SUB=${N_SUB:-"2"}
N_SWEEP=${N_SWEEP:-"4"}
# RK parameters 
NS=${NS:-"3"}
METHOD=${METHOD:-"1"}
date > ${CASE}.log

for DT in $RANGE_DT; do

    echo "========================================================================="
    echo "dt =" $DT
    echo
    sed -e "s/<np1>/$NP1/g" \
        -e "s/<np2>/$NP2/g" \
        -e "s/<np3>/$NP3/g" \
        -e "s/<ep1>/$EP/g" \
        -e "s/<ep2>/$EP/g" \
        -e "s/<ep3>/$EP3/g" \
        -e "s/<po_u>/$PO_U/g" \
        -e "s/<po_p>/$PO_P/g" \
        -e "s/<po_q>/$PO_Q/g" \
        -e "s/<t_end>/$T_END/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<n_sub>/$N_SUB/g" \
        -e "s/<n_sweep>/$N_SWEEP/g" \
        -e "s/<ns>/$NS/g" \
        -e "s/<method>/$METHOD/g" \
        isp_flow__sdc_test.tmpl > \
        isp_flow__sdc_test.prm


    ${EXEC} ${PROGRAM} 2>&1 | tee -a ${CASE}.log

done

grep -e "#      t" -m 1 ${CASE}.log >  ${CASE}.dat
grep -e "#last#$"       ${CASE}.log >> ${CASE}.dat

date >> ${CASE}.log

