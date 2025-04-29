#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

# test case
CASE="convergence_dt"

# mpirun script
MPIRUN=${MPIRUN:-"mpirun"}

# path to program
PROG="../ins_integrator_3d_test"

# kinematic shear viscosity, used for scaling the stabilizing bulk viscosity
NU=0.01

# number of partitions in directions 1-2
NP=${NP:-"1"}
NPROC=$((${NP} * ${NP}))

# elements per partition in directions 1-2
EP=${EP:-"4"}

# domain extension in direction 3
L3=$(bc -l <<< "1/(${NP} * ${EP})")

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U=${PO_U:-"5"}
PO_P=$((${PO_U} - 1))
PO_Q=$((3*${PO_U}/2 + 1))

# time integration
TIME_METHOD=${TIME_METHOD:-"2"}
T_END=${T_END:-"0.25"}

# stabilizing bulk viscosity
MU=$(bc -l <<< "sqrt(${NU}^2 + 4*${L3}^2)")

# max time step size and max number of time step subdivisions per series
DT_MAX=${DT_MAX:-"0.015625"}
ST_MAX=${ST_MAX:-"8"}

date > ${CASE}.log

for((s=0; s<=ST_MAX; s++)); do

    DT=$(bc -l <<< "${DT_MAX}/sqrt(2^${s})")

    echo "========================================================================="
    echo "dt =" $DT ", s = "${s}"/"${ST_MAX}
    echo
    sed -e "s/<l3>/$L3/g" \
        -e "s/<mu>/$MU/g" \
        -e "s/<np1>/$NP/g" \
        -e "s/<np2>/$NP/g" \
        -e "s/<ep1>/$EP/g" \
        -e "s/<ep2>/$EP/g" \
        -e "s/<po_u>/$PO_U/g" \
        -e "s/<po_p>/$PO_P/g" \
        -e "s/<po_q>/$PO_Q/g" \
        -e "s/<t_end>/$T_END/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<time_method>/$TIME_METHOD/g" \
        test_convergence.tmpl > \
        test_convergence.prm

    ${MPIRUN} -n ${NPROC} ${PROG} test_convergence 2>&1 | tee -a ${CASE}.log

done

grep -e "#\s\s*t" -m 1 ${CASE}.log >  ${CASE}.dat
grep -e "#last#"       ${CASE}.log >> ${CASE}.dat

date >> ${CASE}.log
