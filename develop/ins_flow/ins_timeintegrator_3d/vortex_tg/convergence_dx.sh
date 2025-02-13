#!/usr/bin/env bash

# execute within slurm script or stand-alone using
# `bash -l convergence_dt.sh`

# test case
CASE="convergence_dx"

# mpirun script
MPIRUN=${MPIRUN:-"mpirun"}

# path to program
PROGRAM="../ins_timeintegrator_3d_test"

# kinematic shear viscosity, used for scaling the stabilizing bulk viscosity
NU=0.01

# test configurations
if [[ $HPC ]]
then
    # for HPC system
    NP=( 2 4 4 )  # number of partitions in directions 1-2 
    EP=( 1 1 2 )  # number of elements per partition in directions 1-2 
else
    # for personal computers or laptops
    NP=( 1 1 1 )  # number of partitions in directions 1-2 
    EP=( 2 4 8 )  # number of elements per partition in directions 1-2 
fi
NC=${#NP[@]}  # number of configurations used

# polynomial orders for u and p and for integrating the nonlinear terms
PO_U=${PO_U:-"5"}
PO_P=$((${PO_U} - 1))
PO_Q=$((3*${PO_U}/2 + 1))

# time integration
TIME_METHOD=${TIME_METHOD:-"2"}
T_END=${T_END:-"0.25"}

# time step, must be small enough to to neglect temporal error
DT=${DT:-"5e-5"}

date > ${CASE}.log

# test configurations 0 .. NC-1 
for i in $(seq 0 $(($NC - 1))) ; do

    # domain length in direction 3
    L3=$(bc -l <<< "1/(${NP[$i]} * ${EP[$i]})")
    MU=$(bc -l <<< "sqrt(${NU}^2 + 4*${L3}^2)")
    NPROC=$((${NP[$i]} * ${NP[$i]}))

    echo "========================================================================="
    echo "i =" ${i}
    echo
    sed -e "s/<l3>/$L3/g" \
        -e "s/<mu>/$MU/g" \
        -e "s/<np1>/${NP[$i]}/g" \
        -e "s/<np2>/${NP[$i]}/g" \
        -e "s/<ep1>/${EP[$i]}/g" \
        -e "s/<ep2>/${EP[$i]}/g" \
        -e "s/<po_u>/$PO_U/g" \
        -e "s/<po_p>/$PO_P/g" \
        -e "s/<po_q>/$PO_Q/g" \
        -e "s/<t_end>/$T_END/g" \
        -e "s/<dt>/$DT/g" \
        -e "s/<time_method>/$TIME_METHOD/g" \
        test_convergence.tmpl > \
        test_convergence.prm

    ${MPIRUN} -n ${NPROC} ${PROGRAM} test_convergence 2>&1 | tee -a ${CASE}.log

done

grep -e "#\s\s*t" -m 1 ${CASE}.log >  ${CASE}.dat
grep -e "#last#"       ${CASE}.log >> ${CASE}.dat

date >> ${CASE}.log
