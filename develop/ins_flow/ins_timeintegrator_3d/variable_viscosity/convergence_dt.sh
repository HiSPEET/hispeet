#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

# path to program
PROG="../ins_timeintegrator_3d_test"

# test case
CASE="convergence_dt"

# number of partitions in directions 1-3
NP1="1"
NP2="1"
NP3="1"
NP=$((${NP1} * ${NP2} * ${NP3}))

# use preset execution command, or mpirun, if not set
EXEC=${EXEC:-"mpirun -n ${NP}"}

# elements per partition in directions 1-2
EP1=${EP1:-"2"}
EP2=${EP2:-"2"}
EP3=${EP3:-"2"}

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U="16"
PO_P="15"

# time integration
T_END=${T_END:-"0.25"}
TIME_METHOD=${TIME_METHOD:-"1"}

# investigated range of time step sizes, dt = dt_0 / 2^(k/2), k = k_min, ... k_max
DT_MAX=${DT_MAX:-"0.025"}
ST_MAX=${ST_MAX:-"8"}

date > ${CASE}.log

  for((s=0; s<=ST_MAX; s++)); do

    DT=$(bc -l <<< "${DT_MAX}/sqrt(2^${s})")

    echo "========================================================================="
    echo "dt =" $DT ", s = "${s}"/"${ST_MAX}
    echo
    sed -e "s/<np1>/$NP1/g" \
        -e "s/<np2>/$NP2/g" \
        -e "s/<np3>/$NP3/g" \
        -e "s/<ep1>/$EP1/g" \
        -e "s/<ep2>/$EP2/g" \
        -e "s/<ep3>/$EP3/g" \
        -e "s/<po_u>/$PO_U/g" \
        -e "s/<po_p>/$PO_P/g" \
        -e "s/<t_end>/$T_END/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<time_method>/$TIME_METHOD/g" \
        test_convergence.tmpl > \
        test_convergence_dt.prm

    ${EXEC} ${PROG} 2>&1 test_convergence_dt | tee -a ${CASE}.log

done

grep -e "#      t" -m 1 ${CASE}.log >  ${CASE}.dat
grep -e "#last#$"       ${CASE}.log >> ${CASE}.dat

date >> ${CASE}.log
