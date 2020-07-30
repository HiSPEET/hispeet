#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

# test case
CASE="convergence_dt"

# path to program
PROGRAM="../isp_flow__sdc_test"

# number of partitions in directions 1-2
NP1=${NP1:-"1"}
NP2=${NP2:-"1"}
NP=$((${NP1} * ${NP2}))

# use preset execution command, or mpirun, if not set
MPIRUN=${MPIRUN:-"mpirun"}
EXEC=${EXEC:-"${MPIRUN} -n ${NP}"}

# elements per partition in directions 1-2
EP=${EP:-"4"}

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U="8"
PO_P="7"
PO_Q="12"

# time integration
T_END=${T_END:-"0.1"}
TIME_METHOD=${TIME_METHOD:-"5"}
N_SUB=${N_SUB:-"1"}
N_SWEEP=${N_SWEEP:-"1"}

# max time step size and max number of time step subdivisions per series
DT_MAX=${DT_MAX:-"0.2"}
ST_MAX=${ST_MAX:-"8"}

date > ${CASE}.log

for((s=0; s<=ST_MAX; s++)); do

    DT=$(bc -l <<< "${DT_MAX}/sqrt(2^${s})")

    echo "========================================================================="
    echo "dt =" $DT ", s = "${s}"/"${ST_MAX}
    echo
    sed -e "s/<np1>/$NP1/g" \
        -e "s/<np2>/$NP2/g" \
        -e "s/<ep1>/$EP/g" \
        -e "s/<ep2>/$EP/g" \
        -e "s/<po_u>/$PO_U/g" \
        -e "s/<po_p>/$PO_P/g" \
        -e "s/<po_q>/$PO_Q/g" \
        -e "s/<t_end>/$T_END/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<time_method>/$TIME_METHOD/g" \
        -e "s/<n_sub>/$N_SUB/g" \
        -e "s/<n_sweep>/$N_SWEEP/g" \
        isp_flow__sdc_test.tmpl > \
        isp_flow__sdc_test.prm

    ${EXEC} ${PROGRAM} 2>&1 | tee -a ${CASE}.log

done

grep -e "#      t" -m 1 ${CASE}.log >  ${CASE}.dat
grep -e "#last#$"       ${CASE}.log >> ${CASE}.dat

date >> ${CASE}.log
